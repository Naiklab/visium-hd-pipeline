#!/usr/bin/env Rscript

# =============================================================================
# VISIUM HD ANALYSIS PIPELINE LAUNCHER
# =============================================================================
#
# This script demonstrates how to use the general Visium HD workflow
# with different dataset configurations.
#
# Usage:
#   Rscript run_visium_hd_analysis.R [config_file]
#
# If no config file is specified, it will use the default config_visium_hd.R
#
# =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
config_file <- if (length(args) > 0) args[1] else "config_visium_hd.R"

cat("=============================================================================\n")
cat("VISIUM HD ANALYSIS PIPELINE\n")
cat("=============================================================================\n")
cat("Using configuration file:", config_file, "\n")
cat("Start time:", Sys.time(), "\n")
cat("=============================================================================\n\n")

# Check if configuration file exists
if (!file.exists(config_file)) {
  stop("Configuration file not found: ", config_file)
}

# Source the configuration
source(config_file)

# Load required libraries with informative messages
cat("Loading required libraries...\n")
required_packages <- c(
  "Seurat", "ggplot2", "patchwork", "dplyr", "tidyr", 
  "RColorBrewer", "ggpubr", "jsonlite", "knitr"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing missing package:", pkg, "\n")
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Optional packages (install if available)
optional_packages <- c("Banksy", "SeuratWrappers")
for (pkg in optional_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    library(pkg, character.only = TRUE)
    cat("Loaded optional package:", pkg, "\n")
  } else {
    cat("Optional package not available:", pkg, "\n")
  }
}

# Set up parallel processing
future::plan("multicore", workers = N_CORES)
options('future.globals.maxSize' = MAX_MEMORY)

cat("Configuration loaded and validated.\n")
cat("Parallel processing set up with", N_CORES, "cores.\n\n")

# =============================================================================
# ANALYSIS EXECUTION FUNCTIONS
# =============================================================================

run_analysis <- function(render_report = TRUE) {
  start_time <- Sys.time()
  
  cat("Starting Visium HD analysis...\n")
  cat("Project:", PROJECT_NAME, "\n")
  cat("Samples:", paste(SAMPLE_NAMES, collapse = ", "), "\n\n")
  
  # Create output directories
  create_output_dirs()
  
  if (render_report) {
    # Render the R Markdown report
    cat("Rendering analysis report...\n")
    
    output_file <- file.path(REPORT_DIR, paste0(PROJECT_NAME, "_analysis_report.html"))
    
    # Create a temporary config chunk to inject parameters
    config_chunk <- paste0(
      "```{r config, include=FALSE}\n",
      "# Configuration loaded from: ", config_file, "\n",
      "source('", config_file, "')\n\n",
      "# Override parameters in workflow\n",
      "sample.names <- SAMPLE_NAMES\n",
      "data.base.dir <- DATA_BASE_DIR\n",
      "image.names <- IMAGE_NAMES\n",
      "bin.size <- BIN_SIZE\n",
      "min.counts <- MIN_COUNTS\n",
      "max.mt.percent <- MAX_MT_PERCENT\n",
      "sketch.ncells <- SKETCH_NCELLS\n",
      "cluster.resolutions <- CLUSTER_RESOLUTIONS\n",
      "banksy.lambda <- BANKSY_LAMBDA\n",
      "banksy.k <- BANKSY_K\n",
      "```\n"
    )
    
    # Read the workflow template
    workflow_path <- "General_VisiumHD_Analysis_Workflow.Rmd"
    if (!file.exists(workflow_path)) {
      stop("Workflow template not found: ", workflow_path)
    }
    
    workflow_content <- readLines(workflow_path)
    
    # Find the setup chunk and insert config chunk after it
    setup_idx <- grep("```{r setup", workflow_content)
    if (length(setup_idx) > 0) {
      # Find the end of the setup chunk
      setup_end <- grep("```", workflow_content[-(1:setup_idx[1])])[1] + setup_idx[1]
      
      # Insert config chunk
      new_content <- c(
        workflow_content[1:setup_end],
        "",
        strsplit(config_chunk, "\n")[[1]],
        "",
        workflow_content[(setup_end+1):length(workflow_content)]
      )
      
      # Create temporary workflow file
      temp_workflow <- tempfile(fileext = ".Rmd")
      writeLines(new_content, temp_workflow)
      
      # Render the report
      rmarkdown::render(
        temp_workflow,
        output_file = output_file,
        params = list(
          project_name = PROJECT_NAME,
          config_file = config_file
        )
      )
      
      # Clean up
      unlink(temp_workflow)
      
      cat("Analysis report saved to:", output_file, "\n")
    } else {
      warning("Could not find setup chunk in workflow template")
    }
  }
  
  end_time <- Sys.time()
  duration <- end_time - start_time
  
  cat("\n=============================================================================\n")
  cat("ANALYSIS COMPLETED\n")
  cat("=============================================================================\n")
  cat("Total runtime:", format(duration), "\n")
  cat("End time:", Sys.time(), "\n")
  
  if (render_report) {
    cat("Report location:", file.path(REPORT_DIR, paste0(PROJECT_NAME, "_analysis_report.html")), "\n")
  }
  
  cat("Results saved in:", OUTPUT_BASE_DIR, "\n")
  cat("=============================================================================\n")
}

# =============================================================================
# EXAMPLE CONFIGURATIONS FOR DIFFERENT DATASETS
# =============================================================================

# Function to create example configurations
create_example_configs <- function() {
  
  # Example 1: Mouse wound healing time course
  cat("# Example configuration for mouse wound healing experiment\n", 
      file = "config_mouse_wound_healing.R")
  cat('
# Mouse wound healing Visium HD analysis configuration
PROJECT_NAME <- "Mouse_Wound_Healing_VisiumHD"
PROJECT_DESCRIPTION <- "Time course analysis of wound healing in mouse skin"

SAMPLE_NAMES <- c("CONTROL", "DAY1", "DAY3", "DAY7")
DATA_BASE_DIR <- "data/raw/wound_healing/"
IMAGE_NAMES <- list(
  "CONTROL" = c("unwounded_region"),
  "DAY1" = c("wound_day1"),
  "DAY3" = c("wound_day3"), 
  "DAY7" = c("wound_day7")
)

BIN_SIZE <- c(2)
MIN_COUNTS <- 50
MAX_MT_PERCENT <- 25
SKETCH_NCELLS <- 75000
CLUSTER_RESOLUTIONS <- c(0.1, 0.2, 0.5, 1.0, 1.5)
OPTIMAL_RESOLUTION <- 0.5

SPECIES <- "mouse"
CUSTOM_GENE_SETS <- list(
  "Wound_healing" = c("Tgfb1", "Vegfa", "Pdgfa", "Fgf2"),
  "Inflammation" = c("Il1b", "Tnf", "Il6", "Ccl2", "Cxcl1"),
  "Proliferation" = c("Mki67", "Pcna", "Top2a", "Ccnb1"),
  "ECM_remodeling" = c("Col1a1", "Col3a1", "Fn1", "Mmp2", "Mmp9"),
  "Angiogenesis" = c("Vegfa", "Angpt1", "Tie1", "Kdr")
)

source("config_visium_hd.R", local = FALSE)
', file = "config_mouse_wound_healing.R", append = TRUE)

  # Example 2: Human tumor analysis
  cat('
# Human tumor Visium HD analysis configuration  
PROJECT_NAME <- "Human_Tumor_VisiumHD"
PROJECT_DESCRIPTION <- "Spatial analysis of human tumor heterogeneity"

SAMPLE_NAMES <- c("TUMOR_CORE", "TUMOR_EDGE", "NORMAL_TISSUE")
DATA_BASE_DIR <- "data/raw/tumor_analysis/"
IMAGE_NAMES <- list(
  "TUMOR_CORE" = c("core_region1", "core_region2"),
  "TUMOR_EDGE" = c("invasive_front"),
  "NORMAL_TISSUE" = c("adjacent_normal")
)

BIN_SIZE <- c(2)
MIN_COUNTS <- 100
MAX_MT_PERCENT <- 15
SKETCH_NCELLS <- 100000
CLUSTER_RESOLUTIONS <- c(0.2, 0.5, 1.0, 2.0)
OPTIMAL_RESOLUTION <- 1.0

SPECIES <- "human"
CUSTOM_GENE_SETS <- list(
  "Tumor_markers" = c("MKI67", "TP53", "EGFR", "HER2"),
  "Immune_infiltration" = c("CD3E", "CD8A", "CD4", "FOXP3", "CD68"),
  "Angiogenesis" = c("VEGFA", "ANGPT1", "TIE1", "KDR"),
  "EMT" = c("VIM", "CDH1", "CDH2", "SNAI1", "TWIST1"),
  "Hypoxia" = c("HIF1A", "VEGFA", "LDHA", "PKM")
)

source("config_visium_hd.R", local = FALSE)
', file = "config_human_tumor.R")

  # Example 3: Development/organogenesis
  cat('
# Developmental biology Visium HD analysis configuration
PROJECT_NAME <- "Mouse_Development_VisiumHD" 
PROJECT_DESCRIPTION <- "Spatial transcriptomics of organ development"

SAMPLE_NAMES <- c("E12_5", "E14_5", "E16_5", "E18_5")
DATA_BASE_DIR <- "data/raw/development/"
IMAGE_NAMES <- list(
  "E12_5" = c("embryo_E12_5"),
  "E14_5" = c("embryo_E14_5"),
  "E16_5" = c("embryo_E16_5"),
  "E18_5" = c("embryo_E18_5")
)

BIN_SIZE <- c(2)
MIN_COUNTS <- 30
MAX_MT_PERCENT <- 30
SKETCH_NCELLS <- 60000
CLUSTER_RESOLUTIONS <- c(0.1, 0.3, 0.5, 1.0)
OPTIMAL_RESOLUTION <- 0.5

SPECIES <- "mouse"
CUSTOM_GENE_SETS <- list(
  "Neural_development" = c("Sox2", "Nestin", "Pax6", "Ngn2"),
  "Mesenchymal" = c("Twist1", "Snai1", "Vim", "Fn1"),  
  "Epithelial" = c("Cdh1", "Krt8", "Krt18", "Epcam"),
  "Proliferation" = c("Mki67", "Pcna", "Ccnb1", "Top2a"),
  "Apoptosis" = c("Casp3", "Bax", "Bcl2", "P53")
)

source("config_visium_hd.R", local = FALSE)
', file = "config_development.R")

  cat("Created example configuration files:\n")
  cat("- config_mouse_wound_healing.R\n")
  cat("- config_human_tumor.R\n") 
  cat("- config_development.R\n\n")
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# If running interactively, provide helpful information
if (interactive()) {
  cat("Interactive mode detected.\n")
  cat("Available functions:\n")
  cat("- run_analysis(render_report = TRUE): Run the full analysis\n")
  cat("- create_example_configs(): Create example configuration files\n")
  cat("- validate_config(): Validate current configuration\n")
  cat("- create_output_dirs(): Create output directories\n\n")
  
  cat("Current configuration:\n")
  validate_config()
  
} else {
  # Running as script - execute analysis
  run_analysis(render_report = TRUE)
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Function to check dataset integrity
check_dataset <- function() {
  cat("Checking dataset integrity...\n")
  
  issues <- c()
  
  for (sample in SAMPLE_NAMES) {
    sample_dir <- file.path(DATA_BASE_DIR, sample)
    
    if (!dir.exists(sample_dir)) {
      issues <- c(issues, paste("Sample directory not found:", sample_dir))
    } else {
      # Check for required files
      required_files <- c("filtered_feature_bc_matrix.h5", "spatial/tissue_positions.csv")
      
      for (req_file in required_files) {
        full_path <- file.path(sample_dir, req_file)
        if (!file.exists(full_path)) {
          issues <- c(issues, paste("Required file not found:", full_path))
        }
      }
    }
  }
  
  if (length(issues) == 0) {
    cat("Dataset integrity check passed!\n")
    return(TRUE)
  } else {
    cat("Dataset integrity issues found:\n")
    for (issue in issues) {
      cat("-", issue, "\n")
    }
    return(FALSE)
  }
}

# Function to estimate memory requirements
estimate_memory <- function() {
  cat("Estimating memory requirements...\n")
  
  # Rough estimates based on typical Visium HD datasets
  spots_per_sample <- 1000000  # Typical for 2μm bins
  genes_per_sample <- 30000
  
  total_spots <- length(SAMPLE_NAMES) * spots_per_sample
  
  # Memory estimates (rough)
  raw_data_gb <- (total_spots * genes_per_sample * 4) / (1024^3)  # 4 bytes per value
  analysis_overhead <- raw_data_gb * 3  # 3x overhead for analysis objects
  
  total_memory_gb <- raw_data_gb + analysis_overhead
  
  cat("Estimated memory requirements:\n")
  cat("- Raw data:", round(raw_data_gb, 1), "GB\n")
  cat("- Analysis overhead:", round(analysis_overhead, 1), "GB\n") 
  cat("- Total estimated:", round(total_memory_gb, 1), "GB\n")
  cat("- Recommended RAM:", round(total_memory_gb * 1.5, 1), "GB\n\n")
  
  # Check against current settings
  current_limit_gb <- MAX_MEMORY / (1024^3)
  cat("Current memory limit:", round(current_limit_gb, 1), "GB\n")
  
  if (current_limit_gb < total_memory_gb) {
    cat("WARNING: Current memory limit may be insufficient!\n")
    cat("Consider increasing MAX_MEMORY in your config file.\n")
  }
}

cat("Visium HD Analysis Pipeline loaded successfully!\n")
cat("Run run_analysis() to start the analysis or create_example_configs() for examples.\n")