# =============================================================================
# VISIUM HD ANALYSIS CONFIGURATION FILE
# =============================================================================
# 
# This file contains all the parameters needed to customize the Visium HD
# analysis workflow for your specific dataset. Modify the values below
# and source this file in your R analysis script.
#
# Usage: source("config_visium_hd.R")
#
# =============================================================================

# =============================================================================
# DATASET INFORMATION
# =============================================================================

# Project name and description
PROJECT_NAME <- "My_VisiumHD_Project"
PROJECT_DESCRIPTION <- "Description of your Visium HD experiment"

# Sample information
SAMPLE_NAMES <- c("SAMPLE1", "SAMPLE2")  # Update with your actual sample names
DATA_BASE_DIR <- "data/raw/"             # Path to your raw data directory

# Image names for each sample (if multiple images per sample)
# Format: list("SAMPLE1" = c("IMAGE1", "IMAGE2"), "SAMPLE2" = c("IMAGE3"))
IMAGE_NAMES <- list(
  "SAMPLE1" = c("IMAGE1", "IMAGE2"),
  "SAMPLE2" = c("IMAGE3", "IMAGE4")
)

# Bin size for analysis (typically 2, 8, or 16 micrometers)
BIN_SIZE <- c(2)

# =============================================================================
# QUALITY CONTROL PARAMETERS
# =============================================================================

# Minimum UMI counts per spot (spots below this will be filtered out)
MIN_COUNTS <- 50

# Maximum mitochondrial gene percentage (spots above this will be filtered out)
MAX_MT_PERCENT <- 20

# Minimum number of features per spot (optional additional filter)
MIN_FEATURES <- 10

# Maximum number of features per spot (optional additional filter)
MAX_FEATURES <- Inf

# =============================================================================
# ANALYSIS PARAMETERS
# =============================================================================

# Number of cells to use for sketching (for large datasets)
SKETCH_NCELLS <- 50000

# Clustering resolutions to test
CLUSTER_RESOLUTIONS <- c(0.2, 0.5, 1.0, 2.0)

# Optimal resolution for detailed analysis (will be selected from above)
OPTIMAL_RESOLUTION <- 0.5

# Number of variable features to identify
N_VARIABLE_FEATURES <- 5000

# PCA dimensions to use for clustering and UMAP
PCA_DIMS <- 1:50

# =============================================================================
# BANKSY SPATIAL ANALYSIS PARAMETERS
# =============================================================================

# BANKSY lambda values (balance between expression and spatial information)
# 0 = expression only, 1 = spatial only, 0.2 = recommended balance
BANKSY_LAMBDA <- c(0.2, 0.8)

# BANKSY k-neighbors (number of neighboring spots to consider)
BANKSY_K <- c(15, 20)

# =============================================================================
# VISUALIZATION PARAMETERS
# =============================================================================

# Point size factor for spatial plots
SPATIAL_PT_SIZE <- 2.5

# Label size for cluster labels
LABEL_SIZE <- 4

# Color palette for clusters (optional)
CLUSTER_COLORS <- NULL  # Set to NULL for automatic colors

# Figure dimensions for saving plots
FIGURE_WIDTH <- 12
FIGURE_HEIGHT <- 10

# =============================================================================
# COMPUTATIONAL PARAMETERS
# =============================================================================

# Number of CPU cores to use for parallel processing
N_CORES <- future::availableCores()

# Maximum memory for future operations (in bytes)
MAX_MEMORY <- 2*10^12

# =============================================================================
# OUTPUT DIRECTORIES
# =============================================================================

# Base output directory
OUTPUT_BASE_DIR <- "output"

# Subdirectories (will be created automatically)
PLOT_DIR <- file.path(OUTPUT_BASE_DIR, "plots")
TABLE_DIR <- file.path(OUTPUT_BASE_DIR, "tables")
RDS_DIR <- file.path(OUTPUT_BASE_DIR, "RDS-files")
REPORT_DIR <- file.path(OUTPUT_BASE_DIR, "reports")

# =============================================================================
# GENE SETS AND ANNOTATIONS (OPTIONAL)
# =============================================================================

# Mitochondrial gene pattern (species-specific)
MT_PATTERN <- "^mt-|^MT-"  # Use "^MT-" for human, "^mt-" for mouse

# Ribosomal gene pattern
RIBO_PATTERN <- "^Rpl|^Rps|^RPL|^RPS"

# Hemoglobin gene pattern
HB_PATTERN <- "^Hb|^HB"

# Custom gene sets for pathway analysis (optional)
CUSTOM_GENE_SETS <- list(
  "Inflammatory" = c("Il1b", "Tnf", "Il6", "Ccl2"),
  "Proliferation" = c("Mki67", "Pcna", "Top2a", "Ccnb1"),
  "ECM" = c("Col1a1", "Col3a1", "Fn1", "Lam1")
)

# =============================================================================
# SPECIES-SPECIFIC SETTINGS
# =============================================================================

SPECIES <- "mouse"  # Options: "mouse", "human"

# Gene nomenclature adjustments based on species
if (SPECIES == "human") {
  MT_PATTERN <- "^MT-"
  RIBO_PATTERN <- "^RPL|^RPS"
  HB_PATTERN <- "^HB"
} else if (SPECIES == "mouse") {
  MT_PATTERN <- "^mt-"
  RIBO_PATTERN <- "^Rpl|^Rps"
  HB_PATTERN <- "^Hb"
}

# =============================================================================
# ADVANCED PARAMETERS
# =============================================================================

# SCTransform parameters
SCT_VFCELLS <- NULL  # Number of cells to use for variance calculation
SCT_NFEATURES <- 3000  # Number of features to return

# Integration parameters (if combining multiple samples)
INTEGRATION_METHOD <- "CCA"  # Options: "CCA", "RPCA", "harmony"
INTEGRATION_DIMS <- 1:30

# Differential expression parameters
DE_TEST <- "wilcox"  # Statistical test for marker detection
DE_MIN_PCT <- 0.1    # Minimum percentage of cells expressing the gene
DE_LOGFC_THRESHOLD <- 0.25  # Minimum log fold change

# =============================================================================
# FUNCTION TO VALIDATE CONFIGURATION
# =============================================================================

validate_config <- function() {
  cat("Validating configuration...\n")
  
  # Check required parameters
  if (length(SAMPLE_NAMES) == 0) {
    stop("SAMPLE_NAMES must be specified")
  }
  
  if (!dir.exists(dirname(DATA_BASE_DIR))) {
    warning("DATA_BASE_DIR parent directory does not exist: ", dirname(DATA_BASE_DIR))
  }
  
  if (MIN_COUNTS <= 0) {
    stop("MIN_COUNTS must be positive")
  }
  
  if (MAX_MT_PERCENT <= 0 || MAX_MT_PERCENT > 100) {
    stop("MAX_MT_PERCENT must be between 0 and 100")
  }
  
  if (SKETCH_NCELLS <= 0) {
    stop("SKETCH_NCELLS must be positive")
  }
  
  # Check BANKSY parameters
  if (any(BANKSY_LAMBDA < 0) || any(BANKSY_LAMBDA > 1)) {
    stop("BANKSY_LAMBDA values must be between 0 and 1")
  }
  
  if (any(BANKSY_K <= 0)) {
    stop("BANKSY_K values must be positive")
  }
  
  cat("Configuration validation completed successfully!\n")
  
  # Print summary
  cat("\n=== CONFIGURATION SUMMARY ===\n")
  cat("Project:", PROJECT_NAME, "\n")
  cat("Samples:", length(SAMPLE_NAMES), "\n")
  cat("Species:", SPECIES, "\n")
  cat("Min counts:", MIN_COUNTS, "\n")
  cat("Max MT%:", MAX_MT_PERCENT, "\n")
  cat("Sketch cells:", SKETCH_NCELLS, "\n")
  cat("Clustering resolutions:", paste(CLUSTER_RESOLUTIONS, collapse = ", "), "\n")
  cat("BANKSY lambda:", paste(BANKSY_LAMBDA, collapse = ", "), "\n")
  cat("BANKSY k:", paste(BANKSY_K, collapse = ", "), "\n")
  cat("==============================\n\n")
}

# =============================================================================
# FUNCTION TO CREATE OUTPUT DIRECTORIES
# =============================================================================

create_output_dirs <- function() {
  dirs_to_create <- c(OUTPUT_BASE_DIR, PLOT_DIR, TABLE_DIR, RDS_DIR, REPORT_DIR)
  
  for (dir in dirs_to_create) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      cat("Created directory:", dir, "\n")
    }
  }
}

# =============================================================================
# AUTO-EXECUTION
# =============================================================================

# Automatically validate configuration and create directories when sourced
if (interactive()) {
  validate_config()
  create_output_dirs()
}

# =============================================================================
# TEMPLATE FUNCTIONS FOR COMMON CUSTOMIZATIONS
# =============================================================================

# Function to update sample information
update_samples <- function(sample_names, data_dir, image_names = NULL) {
  SAMPLE_NAMES <<- sample_names
  DATA_BASE_DIR <<- data_dir
  if (!is.null(image_names)) {
    IMAGE_NAMES <<- image_names
  }
  cat("Updated sample configuration\n")
}

# Function to update QC thresholds
update_qc_thresholds <- function(min_counts = NULL, max_mt = NULL, min_features = NULL) {
  if (!is.null(min_counts)) MIN_COUNTS <<- min_counts
  if (!is.null(max_mt)) MAX_MT_PERCENT <<- max_mt
  if (!is.null(min_features)) MIN_FEATURES <<- min_features
  cat("Updated QC thresholds\n")
}

# Function to update clustering parameters
update_clustering <- function(resolutions = NULL, optimal_res = NULL, sketch_cells = NULL) {
  if (!is.null(resolutions)) CLUSTER_RESOLUTIONS <<- resolutions
  if (!is.null(optimal_res)) OPTIMAL_RESOLUTION <<- optimal_res
  if (!is.null(sketch_cells)) SKETCH_NCELLS <<- sketch_cells
  cat("Updated clustering parameters\n")
}

cat("Visium HD configuration loaded successfully!\n")
cat("Modify the parameters above and re-source this file to update settings.\n")