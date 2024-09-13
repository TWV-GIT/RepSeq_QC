# Versioning and libraries------------------------------------------------------
script_version = '1.4'
# 1.1: added peak-detection reporting and stock molarity and concentration
# 1.2: added size-labels to vertical lines
# 1.3: changed export file name.

library(ggplot2)
library(bioanalyzeR)
library(gridExtra)
library(grid)
library(dplyr)

# set parameters: --------------------------------------------------------------
annotations_file <- rstudioapi::selectFile(
  caption = "Select the Annotations file",
  filter = "CSV Files (*.csv)",
  existing = TRUE)

XML_file <- rstudioapi::selectFile(
  caption = "Select the .xml file",
  filter = "CSV Files (*.XML)",
  existing = TRUE)

CSV_file <- rstudioapi::selectFile(
  caption = "Select the Electropherogram .csv file",
  filter = "CSV Files (*.csv)",
  existing = TRUE)

output_dir <- rstudioapi::selectDirectory(
  caption = "Select output directory",
  label = "Select",
  path = getwd())

standard_sample_order <- TRUE
exclude_wells <- c('')
dilution <- 10

BCR_refline <- 635
TCR_refline <- 660
IGH_refline <- 700
# save_plots()

# set standard annotations -----------------------------------------------------
standard_annotations <- df <- data.frame(
  well.number = c("B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2", "C2", "D2", "E2","F2","G2", "H2"),
  sample.day = c("D0","D0","D0","D8","D8","D8","D10","D10","D10","D14","D14","D14","D28" ,"D28" ,"D28"),
  sample.type = c("IGH" ,"BCR" ,"TCR" ,"IGH" ,"BCR" ,"TCR" ,"IGH" ,"BCR" ,"TCR" ,"IGH" ,"BCR" ,"TCR" ,"IGH" ,"BCR" ,"TCR")
)


custom_annotations <- df <- data.frame(
  well.number = c("B1", "C1", "D1", "E1", "F1", "G1", "H1", "A2", "B2", "C2", "D2", "E2","F2","G2", "H2"),
  sample.day = c("D14","D14","D14","D28","D28","D28","D28","D28","D28","D0","D8","D14","D28" ,"D28" ,"D28"),
  sample.type = c("IGH" ,"BCR" ,"TCR" ,"IGH" ,"BCR" ,"TCR" ,"IGH" ,"BCR" ,"TCR" ,"BCR" ,"BCR" ,"TCR" ,"IGH" ,"BCR" ,"TCR")
)
# define functions -------------------------------------------------------------
save_plots <- function(){
  pdf(file = paste0(output_dir, '/', run_id, "_library_size_plot.PDF"), width = 8, height = 6)
  plot(1, type="n", xlab="", ylab="", xlim=c(0,10), ylim=c(0,10), axes=FALSE) # Create an empty plot with no axes
  lines <- c(paste0("PDF generated with the NEBNext_TapeStation_QC_", script_version, ".R script."),
             paste0("Run ID: ", run_id),
             report_text,
             paste0("PDF generated on ", format(Sys.time(), "%Y-%m-%d at %H:%M:%S")))
  text(x=5, y=seq(5, 1, by=-1.2), labels=lines, cex=1) # Increase cex to make the text larger
  print(q)
  print(pl)
  print(pl2)
  dev.off()}
# read -------------------------------------------------------------------------
if (standard_sample_order == TRUE) {
  annotations <- standard_annotations
  report_text <- "Standard sample annotations were used."
} else { 
  annotations <- read.csv(annotations_file)
  report_text <- paste0("Sample annotations file: ", annotations_file)
  }
electrophoresis <- read.tapestation(
  XML_file,
  CSV_file,
  method = "hyman",
  extrapolate = FALSE
)
# process ----------------------------------------------------------------------
run_id <- gsub(".*/|\\.xml", "", XML_file)

an_electrophoresis <- annotate.electrophoresis(
  electrophoresis,
  annotations,
  header = TRUE,
  row.names = NULL,
  sep = "\t",
  stringsAsFactors = TRUE,
)

f_electrophoresis <- subset(an_electrophoresis, !(well.number %in% exclude_wells))

table_data <- subset(an_electrophoresis[["samples"]], select = c('well.number', 'sample.name', 'sample.type', 'sample.day'))

table_data$molarity <- dilution*round(integrate.custom(
  an_electrophoresis,
  lower.bound = 500,
  upper.bound = 800,
  bound.variable = "length",
  sum.variable = "molarity"
), 0)

table_data$concentration <- dilution*round(integrate.custom(
  an_electrophoresis,
  lower.bound = 500,
  upper.bound = 800,
  bound.variable = "length",
  sum.variable = "concentration"
), 0)

# plot -------------------------------------------------------------------------
# plot table
my_text <- paste0("Wells excluded from plotting: ", paste(exclude_wells, collapse = ", "),
                  '.\n Molarity and concentration values are for fragments of size 500 to 800 bp.',
                  '\nMolarity is expressed in pM. Concentration is expressed in pg/µl.', 
                  '\nThese are values for undiluted samples (considering ', dilution, 'X dilution before analysis).')

p <- ggplot() +
  theme_void() +  # Remove axes and labels
  xlim(0, 16) +
  ylim(0, 16)

tab <- tableGrob(table_data, rows = NULL, theme = ttheme_minimal(base_size = 10))

title_grob <- textGrob(my_text, gp = gpar(fontsize = 10, fontface = "bold"))

combined <- arrangeGrob(title_grob, tab, ncol = 1, heights = c(0.1, 0.4))

q <- p + annotation_custom(grob = combined, xmin = 0, xmax = 16, ymin = 0, ymax = 16)
print(q)

# plot charts
p <- qplot.electrophoresis(
  f_electrophoresis,
  x = "length",
  y = "concentration",
  color = sample.day,
  log = "",
  normalize = "none",
  facets = ~sample.type,
  margins = FALSE,
  scales = "fixed",
  geom = c("line", "area"),
  include.ladder = FALSE,
  include.markers = FALSE,
  lower.marker.spread = 10,
  xlim = c(200, 1000),
  ylim = c(NA, NA),
  show.peaks = "none",
  region.alpha = 0.2,
  area.alpha = 0.2,
  title = 'NEBNext IS (human) library size distribution for undiluted samples',
  xlab = NULL,
  ylab = NULL
) + theme_minimal() +
  theme(panel.spacing = unit(2, "lines")) +
  scale_color_discrete(name = "Timepoint")

# scale based on dilution
conc_df <- p[["data"]]
conc_df <- conc_df %>% mutate(y.scaled = y.scaled * dilution)
  
p[["data"]] <- conc_df

# Add the vertical lines and values
vlines <- data.frame(
  sample.type = c("IGH", "BCR", "TCR"),
  xintercept = c(IGH_refline, BCR_refline, TCR_refline),
  linetype = c("dashed", "dashed", "dashed"),
  color = c("black", "black",  "black"),
  alpha = c(0.3, 0.3, 0.3)
)
vlines$text_label <- c(IGH_refline, BCR_refline, TCR_refline)

pl <- p + geom_vline(data = vlines, aes(xintercept = xintercept), linetype = vlines$linetype, 
                     color = vlines$color, alpha = vlines$alpha) +
          geom_text(data = vlines, aes(x = xintercept+20, y = 0, label = text_label),
                    vjust = 1.5, hjust = 0, size = 3.5)  

print(pl)


p2 <- qplot.electrophoresis(
  f_electrophoresis,
  x = "length",
  y = "molarity",
  color = sample.day,
  log = "",
  normalize = "none",
  facets = ~sample.type,
  margins = FALSE,
  scales = "fixed",
  geom = c("line", "area"),
  include.ladder = FALSE,
  include.markers = FALSE,
  lower.marker.spread = 10,
  xlim = c(200, 1000),
  ylim = c(NA, NA),
  show.peaks = "none",
  region.alpha = 0.2,
  area.alpha = 0.2,
  title = 'NEBNext IS (human) library size distribution for undiluted samples',
  xlab = NULL,
  ylab = NULL
) + theme_minimal() +
  theme(panel.spacing = unit(2, "lines")) +
  scale_color_discrete(name = "Timepoint")

# scale based on dilution
mol_df <- p2[["data"]]
mol_df <- mol_df %>% mutate(y.scaled = y.scaled * dilution)

p2[["data"]] <- mol_df


pl2 <- p2 + 
  geom_vline(data = vlines, aes(xintercept = xintercept), linetype = vlines$linetype, 
             color = vlines$color, alpha = vlines$alpha) + 
  geom_text(data = vlines, aes(x = xintercept+20, y = 0, label = text_label),
            vjust = 1.5, hjust = 0, size = 3.5)

print(pl2)



filtered_peaks <- an_electrophoresis$peaks %>%
  filter(length >= 500 & length <= 800) %>%
  group_by(sample.index)

# Find the highest "length" for each "sample.index"
result <- filtered_peaks %>%
  summarise(max_length = max(length))

# Print the result
print(result)

