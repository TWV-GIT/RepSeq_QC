# Versioning, libraries, and instructions---------------------------------------

script_version = '1.2'
# 1.1 <- added support for targets
# 1.2 <- added standard curve, absolute quantification and other improvements

library(reshape2)
library(dplyr)
library(ggplot2)
library(pracma)

# instructions:
# the library wells need to be named according to the picomolarity (e.g. 100 pM = sample name '100')

# Set parameters----------------------------------------------------------------
win_dir <- 'C:\\Users\\tverdonckt\\OneDrive - ITG\\Desktop\\DenMark_VLAIO\\Data\\qPCR\\J9\\24-06-28 - 205_208 - libs'
Amp_results <- "208-205 -  Quantification Amplification Results_SYBR.csv"
Quant_summary <- "208-205 -  Quantification Summary_0.csv"

Cycle_cutoff <- 30
lower_Cycle_cutoff <- 0
Sample_detect_threshold <- 500
wavelength <- 'SYBR' # e.g. '483-533' or '465-510'
Abs_threshold <- 188
baseline <- c(lower_Cycle_cutoff: 3)

Well_filter <- c()

plot_title <- "NEBNext IS library quantification plot"
plot_sizes <- c(8,6) # width, height

# samples that should not be plotted
filternames <- c('NTC', '')
# Name of standard samples (must be equal to standard concentration)
laddernames <- c('100', '10', '1', '0.1')
ladder_amplicon_size <- 399
IGH_amplicon_size <- 700
BCR_amplicon_size <- 635
TCR_amplicon_size <- 660
sample_dilution <- 2000

# Samples to plot uniquely
highlightnames <- c('A')

# save_tables()

# save_plots()

# define functions -------------------------------------------------------------
save_tables <- function(){
  write.table(average_sd, file = paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_Ct_values_", Abs_threshold, "_treshold.tsv"))
  write.table(df_summary, file = paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_Library-concentrations_", Abs_threshold, "_treshold.tsv"))
  
}
save_plots <- function(){
  pdf(file = paste0(format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_library_quantification_plots.PDF"), width = plot_sizes[1], height = plot_sizes[2])
  plot(1, type="n", xlab="", ylab="", xlim=c(0,10), ylim=c(0,10), axes=FALSE) # Create an empty plot with no axes
  lines <- c(paste0("PDF generated with the NEBNext_IS_qPCR_Absolute_Quantification_", script_version, ".R script."),
             paste0("Amplification input file: ", Amp_results),
             paste0("Quantification summary file: ", Quant_summary),
             paste0("The following samples were filtered from plotting: ", paste(filternames, collapse = ", ")),
             paste0("The following wells were filtered from plotting and calculations: ", paste(Well_filter, collapse = ", ")),
             paste0("Threshold value: ", Abs_threshold),
             paste0("Baseline range: ", lower_Cycle_cutoff, " to ", max(baseline), " cycles" ),
             paste0("PDF generated on ", format(Sys.time(), "%Y-%m-%d at %H:%M:%S")))
  text(x=5, y=seq(10, 1, by=-1.2), labels=lines, cex=1) # Increase cex to make the text larger
  print(p_log)
  print(p_linear)
  print(p_s_linear)
  print(reg)
  print(con)
  dev.off()}
# Load datasets ----------------------------------------------------------------
win_lin <- gsub("\\\\", "/", win_dir)
setwd(win_lin)
RFU <- read.csv(Amp_results, header = TRUE, check.names = FALSE)
RFU <- subset(RFU, select = -1)
RFU$Cycle <- as.numeric(rownames(RFU))
RFU_long <- melt(RFU, id.vars = "Cycle", variable.name = "Well", value.name = "fluorescence")
RFU_long$Well <- sprintf("%s%02d", substr(RFU_long$Well, 1, 1), as.integer(substr(RFU_long$Well, 2, nchar(as.character(RFU_long$Well)))))

Cq <- read.csv(Quant_summary, header = TRUE, check.names = FALSE)
Cq <- Cq[, !(names(Cq) == "")]


merged_table <- merge(RFU_long, Cq[,c("Well", "Sample", "Target")], by = "Well") # Include "Target" column
merged_table <- merged_table[order(as.numeric(merged_table$Sample), merged_table$Cycle), ]

# Calculations -----------------------------------------------------------------
#Filter out wells where the maximum fluorescence does not exceed the threshold
data_filtered <- merged_table %>%
  group_by(Sample, Target) %>% # Include "Target" column
  filter(max(fluorescence) > Sample_detect_threshold)

# Filter out the specified wells
data_filtered <- data_filtered[!data_filtered$Well %in% Well_filter,]

# Filter out the specified samples
data_filtered <- data_filtered[!data_filtered$Sample %in% filternames, ]

# Cutoff cycle numbers
data_filtered <- subset(data_filtered, Cycle < Cycle_cutoff)
data_filtered <- subset(data_filtered, Cycle > lower_Cycle_cutoff)

df_baseline <- data_filtered[data_filtered$Cycle %in% baseline,]
baseline_fluorescence <- aggregate(df_baseline$fluorescence, by=list(Well=df_baseline$Well), FUN=mean)
colnames(baseline_fluorescence) <- c("Well", "BaselineFluorescence")
data_filtered <- merge(data_filtered, baseline_fluorescence, by="Well")
data_filtered$AdjustedFluorescence <- data_filtered$fluorescence - data_filtered$BaselineFluorescence

# Calculate Ct values
data_filtered <- data_filtered %>%
  group_by(Well, Target) %>% # Include "Target" column
  mutate(Ct = approx(AdjustedFluorescence, `Cycle`, xout = Abs_threshold)$y)

# Calculate average Ct for each sample
average_ct_values <- data_filtered %>%
  group_by(Sample, Target) %>% # Include "Target" column
  summarize(Avg_Ct = mean(Ct, na.rm = TRUE))

# Calculate standard deviation of Ct values
sd_ct_values <- data_filtered %>%
  group_by(Sample, Target) %>% # Include "Target" column
  summarize(sd = sd(Ct))

average_sd <- merge(average_ct_values, sd_ct_values[,c("Sample", "Target", "sd")], by = c("Sample", "Target")) # Include "Target" column

# Create a new column for modified Sample names with average and sd
average_sd <- average_sd %>%
  mutate(ModifiedSampleName = paste0(Sample, " ", Target, " (", round(Avg_Ct, 2), "±", round(sd, 3), ")"))

# Join data with average_ct_values
data_filtered <- left_join(data_filtered, average_sd, by = c("Sample", "Target")) # Include "Target" column

# filter data for selected samples
data_filtered_select <- data_filtered[data_filtered$Sample %in% highlightnames, ]
average_ct_values_select <- average_ct_values[average_ct_values$Sample %in% highlightnames, ]

# Plot the amplification curves ------------------------------------------------
p_log <- ggplot(data_filtered, aes(x = `Cycle`, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
  geom_line() +
  geom_hline(yintercept = Abs_threshold, linetype = "dashed", color = "black") +
  geom_vline(data = average_ct_values, aes(xintercept = Avg_Ct), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = seq(min(data_filtered$`Cycle`), max(data_filtered$`Cycle`), by = 1)) +
  scale_y_log10(breaks = round(logspace(log10(1), log10(max(data_filtered$fluorescence)), 20), 0)) +
  labs(x = "Cycle", y = paste0("Fluorescence (", wavelength, ")"), title = plot_title,
       subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Log-linear plot. Excluding wells ", 
                         paste(Well_filter, collapse = ","), ". Script version ", script_version)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  annotate("text", x = 1, y = Abs_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
  scale_color_discrete(name = "Samples (Ct val.±SD)")

p_linear <- ggplot(data_filtered, aes(x = `Cycle`, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
  geom_line() +
  geom_hline(yintercept = Abs_threshold, linetype = "dashed", color = "black") +
  geom_vline(data = average_ct_values, aes(xintercept = Avg_Ct), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = seq(min(data_filtered$`Cycle`), max(data_filtered$`Cycle`), by = 1)) +
  scale_y_continuous(breaks = seq(0, max(data_filtered$fluorescence), by = ceiling(max(data_filtered$fluorescence)/10))) +
  labs(x = "Cycle", y = paste0("Fluorescence (", wavelength, ")"), title = plot_title,
       subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Linear-linear plot. Excluding wells ", 
                         paste(Well_filter, collapse = ","), ". Script version ", script_version)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +
  annotate("text", x = 1, y = Abs_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
  scale_color_discrete(name = "Samples (Ct val.±SD)")

p_s_linear <- ggplot(data_filtered_select, aes(x = `Cycle`, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
  geom_line() +
  geom_hline(yintercept = Abs_threshold, linetype = "dashed", color = "black") +
  geom_vline(data = average_ct_values_select, aes(xintercept = Avg_Ct), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = seq(min(data_filtered$`Cycle`), max(data_filtered$`Cycle`), by = 1)) +
  scale_y_continuous(breaks = seq(0, max(data_filtered$fluorescence), by = ceiling(max(data_filtered$fluorescence)/10))) +
  labs(x = "Cycle", y = paste0("Fluorescence (", wavelength, ")"), title = paste(plot_title, 'of selected samples'),
       subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Linear-linear plot. Excluding wells ", 
                         paste(Well_filter, collapse = ","), ". Script version ", script_version)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +
  annotate("text", x = 1, y = Abs_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
  scale_color_discrete(name = "Samples (Ct val.±SD)")

# Plot
p_log
p_linear
p_s_linear

# Determine standard regression-------------------------------------------------
# Cut out the ladder data
data_filtered_meta <- data_filtered %>% 
  select(-fluorescence, -Cycle, -AdjustedFluorescence) %>% distinct()

data_ladder <- data_filtered_meta[data_filtered_meta$Sample %in% laddernames,]
data_ladder$Sample <- as.numeric(data_ladder$Sample)

model <- lm(Ct ~ log(Sample, 10), data = data_ladder)
r_squared <- summary(model)$r.squared
slope <- summary(model)$coefficients[2, 1]
intercept <- summary(model)$coefficients[1,1]
efficiency <- (10^(-1/slope))-1

# calculate concentration of samples based on slope and Ct
data_filtered_meta$avg_ct_fit <- 10^((data_filtered_meta$Avg_Ct-intercept)/slope)

# plot values on linear regression
reg <- ggplot(data_ladder, aes(x=Sample, y=Ct)) +
  geom_smooth(method="lm", col="grey") +
  geom_point() +
  theme_minimal() +
  scale_x_log10() +
  labs(title="Linear regression of standards",
       subtitle = paste("R-squared: ", round(r_squared, 3), "; Slope: ", round(slope, 2),
                        "; Efficiency: ", round(efficiency, 3)*100,"%"),
       x="Concentration (pM)", y="Ct-value") +
  geom_point(data = data_filtered_meta %>% filter(Target != "ladder"), 
             aes(x=avg_ct_fit, y=Avg_Ct, color=Target), position = "identity")
print(reg)
# calculate undiluted concentration---------------------------------------------
# adjust for dilution
data_filtered_meta$undiluted <- (10^((data_filtered_meta$Ct-intercept)/slope))*sample_dilution
# adjust for amplicon size
data_filtered_meta <- data_filtered_meta %>% mutate(concentration = 
                        ifelse(Target =="IGH", undiluted*(!!ladder_amplicon_size/!!IGH_amplicon_size),
                          ifelse(Target =="BCR", undiluted*(!!ladder_amplicon_size/!!BCR_amplicon_size),
                            ifelse(Target =="TCR", undiluted*(!!ladder_amplicon_size/!!TCR_amplicon_size), undiluted))))

# plot results
# Load necessary libraries
# Calculate means and standard errors
df_summary <- data_filtered_meta[!data_filtered_meta$Sample %in% laddernames,] %>%
  group_by(Sample, Target) %>%
  summarise(
    mean_concentration = mean(concentration/1000, na.rm = TRUE),
    se_concentration = sd(concentration/1000, na.rm = TRUE) / sqrt(n())
  )

# Create the bar plot
con <- ggplot(df_summary, aes(x = Sample, y = mean_concentration, fill = Target)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    aes(ymin = mean_concentration - se_concentration, ymax = mean_concentration + se_concentration),
    width = 0.2,
    position = position_dodge(0.9)
  ) +
  labs(title = 'Absolute concentration of libraries ± SE',
       subtitle = paste0('Sample dilution = ', sample_dilution, "X;\nLibrary amplicon size of IGH / BCR / TCR: ", 
                         IGH_amplicon_size, ' / ', BCR_amplicon_size, ' / ', TCR_amplicon_size,
                         '\nStandard amplicon size: 399'),
       x = "Sample", y = "Concentration (nM)") +
  theme_minimal()

print(con)
