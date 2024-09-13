# Versioning and libraries------------------------------------------------------
script_version = '1.2'

library(reshape2)
library(dplyr)
library(ggplot2)
library(pracma)
library(scales)
library(gridExtra)

# Set parameters----------------------------------------------------------------
Curve_top <- 1
Cycle_cutoff <- 27
lower_Cycle_cutoff <- 0
Sample_detect_threshold <- 500
wavelength <- 'SYBR' # e.g. '483-533' or '465-510'
Ct_factor <- 0.666 #fraction of plateau fluorescence used for Ct threshold
Abs_threshold <- 128

Well_filter <- c('',
                 '')
plot_title <- "NEBNext IS library quantification plot"


filternames <- c('100', '50', '12.5', '25', '3.13', '6.25', '1.56','NTC') # samples that should not be plotted (e.g. ladders)

highlightnames <- c('A','B') # Samples to plot uniquely

#save_tables()
#save_plots()

# Load and format data ---------------------------------------------------------
# Select and load files, define output directory
if (rstudioapi::isAvailable()) {
  # Using rstudioapi (works in RStudio on all operating systems)
  RFU_file <- rstudioapi::selectFile(
    caption = "Select the Quantification Amplification Results",
    filter = "Results Files (*.csv)",
    existing = TRUE
  )
} else if (.Platform$OS.type == "windows") {
  # Using file.choose() on Windows
  selected_file <- file.choose()
} else {
  # Fallback for non-Windows systems without RStudio
  cat("Please enter the path to the file: ")
  selected_file <- readline()
}
if (rstudioapi::isAvailable()) {
  # Using rstudioapi (works in RStudio on all operating systems)
  Cq_file <- rstudioapi::selectFile(
    caption = "Select the Quantification Cq Results",
    filter = "Results Files (*.csv)",
    existing = TRUE
  )
} else if (.Platform$OS.type == "windows") {
  # Using file.choose() on Windows
  selected_file <- file.choose()
} else {
  # Fallback for non-Windows systems without RStudio
  cat("Please enter the path to the file: ")
  selected_file <- readline()
}
RFU <- read.csv(RFU_file, header = TRUE, check.names = FALSE)
Cq <- read.csv(Cq_file, header = TRUE, check.names = FALSE)
if (rstudioapi::isAvailable()) {
  # Using rstudioapi (works in RStudio on all operating systems)
  output_dir <- rstudioapi::selectDirectory(
    caption = "Select output directory",
    label = "Select",
    path = getwd()
  )
} else if (.Platform$OS.type == "windows") {
  # Using choose.dir() on Windows
  output_dir <- choose.dir(default = getwd(), caption = "Select output directory")
} else {
  # Fallback for non-Windows systems without RStudio
  cat("Please enter the path to the output directory: ")
  output_dir <- readline()
}

RFU <- subset(RFU, select = -1)
RFU$Cycle <- as.numeric(rownames(RFU))
RFU_long <- melt(RFU, id.vars = "Cycle", variable.name = "Well", value.name = "fluorescence")
RFU_long$Well <- sprintf("%s%02d", substr(RFU_long$Well, 1, 1), as.integer(substr(RFU_long$Well, 2, nchar(as.character(RFU_long$Well)))))

Cq <- subset(Cq, select = -1)
merged_table <- merge(RFU_long, Cq[,c("Well", "Sample", "Target")], by = "Well") # Include "Target" column
merged_table <- merged_table[order(as.numeric(merged_table$Sample), merged_table$Cycle), ]

# Calculations -----------------------------------------------------------------
baseline <- c(lower_Cycle_cutoff: 3)
#Filter out wells where the maximum fluorescence does not exceed the threshold
data_filtered <- merged_table %>%
  group_by(Sample, Target) %>% # Include "Target" column
  filter(max(fluorescence) > Sample_detect_threshold)

# Filter out the specified wells
data_filtered <- data_filtered[!data_filtered$Well %in% Well_filter,]

# Cutoff cycle numbers
data_filtered2 <- subset(data_filtered, Cycle < Cycle_cutoff)
data_filtered2 <- subset(data_filtered2, Cycle > lower_Cycle_cutoff)

# Normalize fluorescence values (Min-Max normalization)
normalize_fluorescence <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

data_normalized <- data_filtered2 %>%
  group_by(Well, Target) %>% # Include "Target" column
  mutate(`fluorescence_normalized` = normalize_fluorescence(fluorescence))

df_baseline <- data_filtered2[data_filtered2$Cycle %in% baseline,]
baseline_fluorescence <- aggregate(df_baseline$fluorescence, by=list(Well=df_baseline$Well), FUN=mean)
colnames(baseline_fluorescence) <- c("Well", "BaselineFluorescence")
data_filtered2 <- merge(data_filtered2, baseline_fluorescence, by="Well")
data_filtered2$AdjustedFluorescence <- data_filtered2$fluorescence - data_filtered2$BaselineFluorescence

# Calculate Ct values
Ct_threshold <- Curve_top*Ct_factor

data_normalized <- data_normalized %>%
  group_by(Well, Target) %>% # Include "Target" column
  mutate(Ct = approx(`fluorescence_normalized`, `Cycle`, xout = Ct_threshold)$y)

data_filtered3 <- data_filtered2 %>%
  group_by(Well, Target) %>% # Include "Target" column
  mutate(Ct = approx(AdjustedFluorescence, `Cycle`, xout = Abs_threshold)$y)

# Calculate average Ct for each sample
average_ct_values_n <- data_normalized %>%
  group_by(Sample, Target) %>% # Include "Target" column
  summarize(Avg_Ct = mean(Ct, na.rm = TRUE))

average_ct_values <- data_filtered3 %>%
  group_by(Sample, Target) %>% # Include "Target" column
  summarize(Avg_Ct = mean(Ct, na.rm = TRUE))

# Calculate standard deviation of Ct values
sd_ct_values_n <- data_normalized %>%
  group_by(Sample, Target) %>% # Include "Target" column
  summarize(sd = sd(Ct))

sd_ct_values <- data_filtered3 %>%
  group_by(Sample, Target) %>% # Include "Target" column
  summarize(sd = sd(Ct))

average_sd_n <- merge(average_ct_values_n, sd_ct_values_n[,c("Sample", "Target", "sd")], by = c("Sample", "Target")) # Include "Target" column
average_sd <- merge(average_ct_values, sd_ct_values[,c("Sample", "Target", "sd")], by = c("Sample", "Target")) # Include "Target" column

# Create a new column for modified Sample names with average and sd
average_sd_n <- average_sd_n %>%
  mutate(ModifiedSampleName = paste0(Sample, " ", Target, " (", round(Avg_Ct, 2), "±", round(sd, 3), ")"))

average_sd <- average_sd %>%
  mutate(ModifiedSampleName = paste0(Sample, " ", Target, " (", round(Avg_Ct, 2), "±", round(sd, 3), ")"))

# Group the samples based on the required Cycles and plot table
wrap_text <- function(text, width) {
  sapply(text, function(x) paste(strwrap(x, width=width), collapse="\n"))
}

Cycle_data_n <- average_sd_n %>%
  arrange(Target, Sample) %>%  # Order by Target, then by Sample
  mutate(Cycles = round(Avg_Ct)) %>%
  group_by(Cycles) %>%
  summarise(Samples = wrap_text(toString(paste(Sample, Target, sep = " ")), width = 42))

table_grob_n <- tableGrob(Cycle_data_n, rows = NULL)

average_sd_n <- average_sd_n %>%
  mutate(Cycles = round(Avg_Ct)) %>%
  left_join(Cycle_data_n, by = "Cycles")
  
# Join data with average_ct_values
data_normalized <- left_join(data_normalized, average_sd_n, by = c("Sample", "Target")) # Include "Target" column

data_filtered3 <- left_join(data_filtered3, average_sd, by = c("Sample", "Target")) # Include "Target" column

data_filtered4 <- data_filtered3[!data_filtered3$Sample %in% filternames, ]
data_normalized_filtered <- data_normalized[!data_normalized$Sample %in% filternames, ]

# Plot the amplification curves ------------------------------------------------
p_log <- ggplot(data_filtered4, aes(x = `Cycle`, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
  geom_line() +
  geom_hline(yintercept = Abs_threshold, linetype = "dashed", color = "black") +
  geom_vline(data = average_ct_values, aes(xintercept = Avg_Ct), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = seq(min(data_filtered3$`Cycle`), max(data_filtered3$`Cycle`), by = 1)) +
  scale_y_log10(breaks = round(logspace(log10(1), log10(max(data_filtered3$fluorescence)), 20), 0)) +
  labs(x = "Cycle", y = paste0("Fluorescence (", wavelength, ")"), title = plot_title,
       subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Log-linear plot. Excluding wells ", 
                         paste(Well_filter, collapse = ","), ". Script version ", script_version)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  annotate("text", x = 1, y = Abs_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
  scale_color_discrete(name = "Samples (Ct val.±SD)")

p_linear <- ggplot(data_filtered4, aes(x = `Cycle`, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
  geom_line() +
  geom_hline(yintercept = Abs_threshold, linetype = "dashed", color = "black") +
  geom_vline(data = average_ct_values, aes(xintercept = Avg_Ct), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = seq(min(data_filtered3$`Cycle`), max(data_filtered3$`Cycle`), by = 1)) +
  scale_y_continuous(breaks = seq(0, max(data_filtered3$fluorescence), by = ceiling(max(data_filtered3$fluorescence)/10))) +
  labs(x = "Cycle", y = paste0("Fluorescence (", wavelength, ")"), title = plot_title,
       subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Linear-linear plot. Excluding wells ", 
                         paste(Well_filter, collapse = ","), ". Script version ", script_version)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +
  annotate("text", x = 1, y = Abs_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
  scale_color_discrete(name = "Samples (Ct val.±SD)")

p_n_linear <- ggplot(data_normalized_filtered, aes(x = `Cycle`, y = `fluorescence_normalized`, group = Well, color = ModifiedSampleName)) +
  geom_line() +
  geom_hline(yintercept = Ct_threshold, linetype = "dashed", color = "black") +
  geom_hline(yintercept = Curve_top, linetype = "dashed", color = "red") +
  geom_vline(data = average_ct_values_n, aes(xintercept = Avg_Ct), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = seq(min(data_normalized$`Cycle`), max(data_normalized$`Cycle`), by = 1)) +
  scale_y_continuous(breaks = seq(0, max(data_normalized$`fluorescence_normalized`), by = 0.1)) +
  labs(x = "Cycle", y = paste0("Normalized Fluorescence (", wavelength, ")"), title = plot_title,
       subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Normalized linear-linear plot. Excluding wells ", 
                         paste(Well_filter, collapse = ","), ". Script version ", script_version)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +
  annotate("text", x = 1, y = Ct_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
  annotate("text", x = 1, y = Curve_top, label = "Plateau fluorescence", hjust = 0.1, vjust = -0.2, color = "red") +
  scale_color_discrete(name = "Samples (Ct val.±SD)")


p_n_linear_grouped <- ggplot(data_normalized_filtered, aes(x = `Cycle`, y = `fluorescence_normalized`, group = Well, color = as.character(Cycles))) +
  geom_line() +
  geom_hline(yintercept = Ct_threshold, linetype = "dashed", color = "black") +
  geom_hline(yintercept = Curve_top, linetype = "dashed", color = "red") +
  geom_vline(data = Cycle_data_n, aes(xintercept = Cycles), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = seq(min(data_normalized$`Cycle`), max(data_normalized$`Cycle`), by = 1)) +
  scale_y_continuous(breaks = seq(0, max(data_normalized$`fluorescence_normalized`), by = 0.1)) +
  labs(x = "Cycle", y = paste0("Normalized Fluorescence (", wavelength, ")"), title = plot_title,
       subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Normalized linear-linear plot. Excluding wells ", 
                         paste(Well_filter, collapse = ","), ". Script version ", script_version)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +
  annotate("text", x = 1, y = Ct_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
  annotate("text", x = 1, y = Curve_top, label = "Plateau fluorescence", hjust = 0.1, vjust = -0.2, color = "red") +
  scale_color_discrete(name = "Samples (Ct val.±SD)") +
  theme(legend.position = "none")

data_filtered3_select <- data_filtered3[data_filtered3$Sample %in% highlightnames, ]
data_normalized_select <- data_normalized[data_normalized$Sample %in% highlightnames, ]
average_ct_values_select <- average_ct_values[average_ct_values$Sample %in% highlightnames, ]
average_ct_values_n_select <- average_ct_values_n[average_ct_values_n$Sample %in% highlightnames, ]


p_s_linear <- ggplot(data_filtered3_select, aes(x = `Cycle`, y = AdjustedFluorescence, group = Well, color = ModifiedSampleName)) +
  geom_line() +
  geom_hline(yintercept = Abs_threshold, linetype = "dashed", color = "black") +
  geom_vline(data = average_ct_values_select, aes(xintercept = Avg_Ct), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = seq(min(data_filtered3$`Cycle`), max(data_filtered3$`Cycle`), by = 1)) +
  scale_y_continuous(breaks = seq(0, max(data_filtered3$fluorescence), by = ceiling(max(data_filtered3$fluorescence)/10))) +
  labs(x = "Cycle", y = paste0("Fluorescence (", wavelength, ")"), title = plot_title,
       subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Linear-linear plot. Excluding wells ", 
                         paste(Well_filter, collapse = ","), ". Script version ", script_version)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +
  annotate("text", x = 1, y = Abs_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
  scale_color_discrete(name = "Samples (Ct val.±SD)")

p_s_n_linear <- ggplot(data_normalized_select, aes(x = `Cycle`, y = `fluorescence_normalized`, group = Well, color = ModifiedSampleName)) +
  geom_line() +
  geom_hline(yintercept = Ct_threshold, linetype = "dashed", color = "black") +
  geom_hline(yintercept = Curve_top, linetype = "dashed", color = "red") +
  geom_vline(data = average_ct_values_n_select, aes(xintercept = Avg_Ct), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = seq(min(data_normalized$`Cycle`), max(data_normalized$`Cycle`), by = 1)) +
  scale_y_continuous(breaks = seq(0, max(data_normalized$`fluorescence_normalized`), by = 0.1)) +
  labs(x = "Cycle", y = paste0("Normalized Fluorescence (", wavelength, ")"), title = plot_title,
       subtitle = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " - Linear-linear plot. Excluding wells ", 
                         paste(Well_filter, collapse = ","), ". Script version ", script_version)) +
  theme_minimal() +
  theme(panel.grid.minor.x = element_blank()) +
  annotate("text", x = 1, y = Ct_threshold, label = "Ct threshold", hjust = 0.1, vjust = -0.2, color = "black") +
  annotate("text", x = 1, y = Curve_top, label = "Plateau fluorescence", hjust = 0.1, vjust = -0.2, color = "red") +
  scale_color_discrete(name = "Samples (Ct val.±SD)")

# Plot
p_log
p_linear
p_n_linear
p_s_linear
p_s_n_linear
grid.arrange(p_n_linear_grouped, table_grob_n, ncol = 2)


# Save output ------------------------------------------------------------------
save_tables <- function(){
  write.table(average_sd, file = paste0(output_dir, '/', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_Ct_values_", Abs_threshold, "_treshold.tsv"))
  write.table(average_sd_n, file = paste0(output_dir, '/', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_Ct_values_norm_", Ct_factor, "_treshold.tsv"))
}
save_plots <- function(){
  pdf(file = paste0(output_dir, '/', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_library_quantification_plots.PDF"), width = 10, height = 7)
  plot(1, type="n", xlab="", ylab="", xlim=c(0,10), ylim=c(0,10), axes=FALSE) # Create an empty plot with no axes
  lines <- c(paste0("PDF generated with the NEBNext_IS_PCR2_Ct_V", script_version, ".R script."),
             paste0("Amplification input file: ", sub("^.*/", "", RFU_file)),
             paste0("Quantification summary file: ", sub("^.*/", "", Cq_file)),
             paste0("The following samples were filtered from plotting: ", paste(filternames, collapse = ", ")),
             paste0("The following wells were filtered from plotting and calculations: ", paste(Well_filter, collapse = ", ")),
             paste0("Threshold value: ", Abs_threshold, "; for normalized plots: ", Ct_factor),
             paste0("Baseline range: ", lower_Cycle_cutoff, " to ", max(baseline), " cycles" ),
             paste0("PDF generated on ", format(Sys.time(), "%Y-%m-%d at %H:%M:%S")))
  text(x=5, y=seq(10, 1, by=-1.2), labels=lines, cex=1) # Increase cex to make the text larger
  print(p_log)
  print(p_linear)
  print(p_n_linear)
  print(p_s_linear)
  print(p_s_n_linear)
  grid.arrange(p_n_linear_grouped, table_grob_n, ncol = 2)
  dev.off()}
