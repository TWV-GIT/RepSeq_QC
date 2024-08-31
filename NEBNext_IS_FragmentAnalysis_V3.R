# Load necessary libraries
library(ggplot2)

# set working directory
win_dir <- 'C:\\Users\\tverdonckt\\OneDrive - ITG\\Desktop\\DenMark_VLAIO\\Data\\Fragment Analyzer\\24-03-01_RUN 189\\Library_BCR_11-38-53'
#win_dir <- 'C:\\Users\\tverdonckt\\OneDrive - ITG\\Desktop\\DenMark_VLAIO\\Data\\Fragment Analyzer\\24-03-01_RUN 189\\Library TCR 12-45-26'
#win_dir <- 'C:\\Users\\tverdonckt\\OneDrive - ITG\\Desktop\\DenMark_VLAIO\\Data\\Fragment Analyzer\\24-03-01_RUN 189\\PCR1 BCR 13-52-00'
#win_dir <- "C:\\Users\\tverdonckt\\OneDrive - ITG\\Desktop\\DenMark_VLAIO\\Data\\Fragment Analyzer\\24-03-01_RUN 189\\PCR1 TCR 14-58-32"
#win_dir <- 'C:\\Users\\tverdonckt\\OneDrive - ITG\\Desktop\\DenMark_VLAIO\\Data\\Fragment Analyzer\\23-12-19 - RUN 147 12-56-21'

win_lin <- gsub("\\\\", "/", win_dir)
setwd(win_lin)

# standard csv: if col sep , and decimal .
data <- read.csv("2023 12 19 12H 56M Electropherogram Libraries.csv", header = TRUE, stringsAsFactors = FALSE)

# European csv: if col sep ; and decimal ,
#data <- read.csv2("2023 12 19 12H 56M Electropherogram PCR.csv", header = TRUE, stringsAsFactors = FALSE)

name_filter <- c('B', 
                 'Blank',
                 'Da B6',
                 'G8',
                 'H8',
                 'KD D8'
                 )
ladder_name <- 'ladder'

ladder_peaks_norm <- c(100, 200, 300, 400, 600, 700, 800, 900, 1200, 1500, 2000, 3000)

# define smear limits
smear <- c(500, 675)

# Replace dots with spaces and colons with underscores in the column names
colnames(data) <- gsub("\\.\\.", ": ", colnames(data))
colnames(data) <- gsub("\\.", " ", colnames(data))
colnames(data) <- gsub("Size: bp ", "Size", colnames(data))

# Reshape data for ggplot
data_long <- tidyr::pivot_longer(data, -'Size', names_to = "Sample", values_to = "Value")

# Extract well position and sample name from Sample column
data_long$Well <- gsub(":.*", "", data_long$Sample)
data_long$Sample_Name <- gsub(".*: ", "", data_long$Sample)

# Filter out samples based on name_filter
filtered_data <- data_long[!data_long$Sample_Name %in% name_filter, ]

# Calculate the average fluorescence of the smaller ladder peaks
average_value <- mean(filtered_data[filtered_data$Sample_Name == ladder_name & filtered_data$Size %in% ladder_peaks_norm,]$Value)

# Normalize the fluorescence towards the ladder
filtered_data$Normalized_values <- filtered_data$Value/average_value

# Calculate highest value for each sample
highest_values <- aggregate(Value ~ Sample_Name, data = filtered_data, FUN = max)
highest_values <- subset(highest_values, Sample_Name != ladder_name)

highest_values_norm <- aggregate(Normalized_values ~ Sample_Name, data = filtered_data, FUN = max)
highest_values_norm <- subset(highest_values_norm, Sample_Name != ladder_name)

filtered_data <- filtered_data %>%
  group_by(Sample_Name) %>%
  mutate(Max_Normalized = max(Normalized_values))
filtered_data <- filtered_data %>%
  mutate(Sample_Name_max = paste(Sample_Name, " (", round(Max_Normalized, 1), ")", sep=""))

# calculate average fragment size for each sample
library(pracma)

filtered_data_smear <- filtered_data %>%
  filter(Size >= min(smear) & Size <= max(smear))

    # get the unique samples
samples <- unique(filtered_data_smear$Sample)

    # for each sample, calculate the value of x that divides the area under the curve into two equal parts
for (sample in samples) {
  # subset the dataframe for the current sample
  df_sample <- filtered_data_smear[filtered_data_smear$Sample == sample, ]
  
  # calculate total area under the curve
  total_area <- pracma::trapz(df_sample$Size, df_sample$Normalized_values)
  
  # define a function for the area from 0 to x
  area_func <- function(x) {
    x_values <- df_sample$Size[df_sample$Size <= x]
    y_values <- df_sample$Normalized_values[df_sample$Size <= x]
    pracma::trapz(x_values, y_values)
  }
  
  # define a function for the difference between the total area / 2 and the area from 0 to x
  diff_func <- function(x) {
    total_area / 2 - area_func(x)
  }
  
  # find the root of the diff_func
  result <- uniroot(diff_func, interval = c(min(df_sample$Size), max(df_sample$Size)))
  
  # print the value of x that divides the area under the curve into two equal parts for the current sample
  print(paste("Sample:", sample, "Size:", result$root))
}


# Plotting using ggplot2 with filtered data and horizontal dotted lines
ggplot(filtered_data, aes(x=Size, y=Normalized_values, color=Sample_Name_max)) +
  geom_line() +
  geom_hline(data = highest_values_norm, aes(yintercept = Normalized_values, linetype = Sample_Name), 
             color = "black", linetype = "dotted") + # Add dotted lines
  labs(x="Fragment Size (bp)", y="Normalized fluorescence", title="NEBNext IS library prep", color="Sample names (max val.)") +
  theme_minimal() +
  scale_color_manual(values = rainbow(length(unique(filtered_data$Sample_Name)))) +
  scale_x_continuous(trans = "log10", limits = c(400, 1200), breaks = c(50, 100, 200,
                                                                        300, 400, 500,
                                                                        600, 700, 800,
                                                                        900, 1000, 1200,
                                                                        1500, 2000, 3000)) +
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, by = 1))



# Optimized for PCR products

filtered_data2 <- filtered_data[filtered_data$Sample_Name != ladder_name,]

ggplot(filtered_data2, aes(x=Size, y=Normalized_values, color=Sample_Name_max)) +
  geom_line() +
  labs(x="Fragment Size (bp)", y="Normalized fluorescence", title="NEBNext IS library prep", color="Sample names (max val.)") +
  theme_minimal() +
  scale_color_manual(values = rainbow(length(unique(filtered_data$Sample_Name)))) +
  scale_x_continuous(trans = "log10", limits = c(400, 1200), breaks = c(50, 100, 200,
                                                                        300, 400, 500,
                                                                        600, 700, 800,
                                                                        900, 1000, 1200,
                                                                        1500, 2000, 3000)) +
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.5, by = 0.1))
