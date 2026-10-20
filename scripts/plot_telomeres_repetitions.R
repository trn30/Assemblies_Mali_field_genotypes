

# Run from the repository root. Generate these inputs before running this script:
# samtools faidx data/reference/PlasmoDB-63_Pfalciparum3D7_Genome.fasta
# samtools faidx data/assemblies/KK2504_assembly_ILRA.fasta
# tidk search -s TTTAGGG -w 10000 -d results/telomeres -o 3D7 data/reference/PlasmoDB-63_Pfalciparum3D7_Genome.fasta
# tidk search -s TTTAGGG -w 10000 -d results/telomeres -o KK2504 data/assemblies/KK2504_assembly_ILRA.fasta


# Load required libraries
library(ggplot2)
library(scales)

# =========================================================================
# 1. LOAD CHROMOSOME BACKGROUNDS (.fai files)
# =========================================================================
# Function to automatically read and format .fai files
prepare_fai <- function(file_path, strain_name) {
  df <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE)[, 1:2]
  colnames(df) <- c("Chr", "Length")

  # Clean chromosome names (removes "_v3" and "Pf3D7_") to get numeric values
  df$Chromosome <- as.numeric(gsub("_v3", "", gsub("Pf3D7_", "", df$Chr)))
  df <- df[!is.na(df$Chromosome) & df$Chromosome >= 1 & df$Chromosome <= 14, ]

  df$Strain <- strain_name
  df$Start <- 0
  return(df)
}

# Add your samples here (Update file paths accordingly)
bg_3D7     <- prepare_fai("data/reference/PlasmoDB-63_Pfalciparum3D7_Genome.fasta.fai", "3D7")
bg_KK2504 <- prepare_fai("data/assemblies/KK2504_assembly_ILRA.fasta.fai", "KK2504")

# Combine all backgrounds into a single data frame
chromosome_backgrounds <- rbind(bg_3D7, bg_KK2504)

# Force the order so chromosome 1 is at the top and 14 at the bottom
chromosome_backgrounds$Chromosome <- factor(chromosome_backgrounds$Chromosome, levels = 14:1)

# =========================================================================
# 2. PROCESS TIDK REPEAT DATA (.tsv files)
# =========================================================================
process_tidk_repeats <- function(tsv_file, strain_name, window_size = 10000) {
  df <- read.delim(tsv_file, header = TRUE, stringsAsFactors = FALSE)

  # Rename standard tidk output columns (Column 2 is the actual coordinate)
  colnames(df)[1:4] <- c("Chr", "Real_Start", "Fwd", "Rev")

  # Calculate total repeats
  df$Total_Repeats <- df$Fwd + df$Rev

  # Filter out noise (e.g., <= 10 repeats)
  df <- subset(df, Total_Repeats > 10)
  if(nrow(df) == 0) return(NULL)

  # Set genomic coordinates based on window size
  df$Start <- df$Real_Start
  df$End <- df$Real_Start + window_size

  # Clean chromosome names to match the background data
  df$Chromosome <- as.numeric(gsub("_v3", "", gsub("Pf3D7_", "", df$Chr)))
  df <- subset(df, !is.na(Chromosome) & Chromosome >= 1 & Chromosome <= 14)

  # Categorize repeats by frequency
  df$Category <- ifelse(df$Total_Repeats >= 100, ">100 repeats",
                 ifelse(df$Total_Repeats >= 30, "30-100 repeats", "10-30 repeats"))
  df$Strain <- strain_name

  # Return cleaned dataframe
  return(df[, c("Strain", "Chromosome", "Start", "End", "Category", "Total_Repeats", "Fwd", "Rev")])
}

# Process the corresponding .tsv files for the samples
telomere_data <- rbind(
  process_tidk_repeats("results/telomeres/3D7_telomeric_repeat_windows.tsv", "3D7", 10000),
  process_tidk_repeats("results/telomeres/KK2504_telomeric_repeat_windows.tsv", "KK2504", 10000)
)

# Set factor levels for plotting order
telomere_data$Chromosome <- factor(telomere_data$Chromosome, levels = 14:1)
telomere_data$Category <- factor(telomere_data$Category, levels = c(">100 repeats", "30-100 repeats", "10-30 repeats"))

# Order data so the highest repeat values are plotted on top
telomere_data <- telomere_data[order(telomere_data$Category, decreasing = TRUE), ]

# =========================================================================
# 3. GENERATE THE FINAL PLOT
# =========================================================================
# Define custom colors
telo_colors <- c(">100 repeats" = "#2ECC71",
                 "30-100 repeats" = "#FFC101",
                 "10-30 repeats" = "#F57C00")

ideogram_plot <- ggplot() +
  # LAYER 1: Chromosome backgrounds (grey bars)
  geom_segment(data = chromosome_backgrounds,
               aes(x = Start, xend = Length, y = Chromosome, yend = Chromosome),
               color = "grey95", linewidth = 4, lineend = "round") +

  # LAYER 2: Telomere repeat marks
  geom_segment(data = telomere_data,
               aes(x = Start, xend = End, y = Chromosome, yend = Chromosome, color = Category),
               linewidth = 3, lineend = "round") +

  # Facet by strain (adjust ncol as needed for more samples)
  facet_wrap(~Strain, ncol = 2, scales = "free_x") +
  scale_color_manual(values = telo_colors, name = "Telomeric repeat\nsequence region", drop = FALSE) +

  # Format X axis to Mbp and Y axis to "chrXX" format
  scale_x_continuous(labels = function(x) paste0(x / 1000000, "Mbp")) +
  scale_y_discrete(labels = function(x) sprintf("chr%02d", as.numeric(as.character(x)))) +

  labs(x = "Length (Mbp)", y = "") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(face = "bold", size = 9, color = "black"),
        strip.text = element_text(face = "bold", size = 13))

# Save outputs
if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)
pdf("results/figures/telomeric_repeats.pdf", width = 10, height = 10)
print(ideogram_plot)
dev.off()

ggsave("results/figures/telomeric_repeats.png", plot = ideogram_plot, width = 12, height = 10)
