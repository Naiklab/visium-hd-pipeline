# Visium HD Analysis Pipeline

A comprehensive, generalizable workflow for analyzing Visium HD spatial transcriptomics data, derived from the wound healing analysis template.

## Overview

This pipeline provides a standardized approach to analyze Visium HD data with the following features:

- **Modular Design**: Easy to customize for different experimental designs
- **Quality Control**: Comprehensive QC metrics and filtering
- **Multiple Clustering Approaches**: Standard clustering + BANKSY spatial clustering
- **Scalable**: Handles large datasets using sketching approaches
- **Reproducible**: Parameterized configuration system
- **Comprehensive Output**: Plots, tables, and analysis reports

## Files in this Pipeline

### Core Files
- `General_VisiumHD_Analysis_Workflow.Rmd` - Main analysis workflow (R Markdown)
- `config_visium_hd.R` - Configuration file with all parameters
- `run_visium_hd_analysis.R` - Pipeline launcher and utility functions
- `submit_template.sh` - LSF/`bsub` job template (Minerva) with the required module + BLAS-thread environment (see FAQ)
- `renv.lock` - Pinned R package versions (R 4.4.1, Bioconductor 3.20) for a reproducible `renv::restore()`
- `README.md` - This file

### Template Files
- `pk_wound_healing_visiumHD_analysis_first_pass.Rmd` - Original template script

> **New here? Read the [FAQ](FAQ.md) first.** It captures the non-obvious HPC, deconvolution, and statistical gotchas that cost real debugging time on actual Visium HD projects.

## Quick Start

### 1. Basic Setup

```bash
# Clone or download the pipeline files to your working directory
# Ensure your data is organized as follows:
data/
├── raw/
│   ├── SAMPLE1/
│   │   ├── filtered_feature_bc_matrix.h5
│   │   ├── spatial/
│   │   │   ├── tissue_positions.csv
│   │   │   ├── scalefactors_json.json
│   │   │   └── tissue_hires_image.png
│   │   └── cloupe_csvs/ (optional, for splitting samples)
│   └── SAMPLE2/
│       └── ... (same structure)
```

### 2. Configure Your Analysis

Edit `config_visium_hd.R` to match your dataset:

```r
# Update these key parameters:
PROJECT_NAME <- "Your_Project_Name"
SAMPLE_NAMES <- c("SAMPLE1", "SAMPLE2", "SAMPLE3")
DATA_BASE_DIR <- "data/raw/"

# Adjust QC thresholds as needed:
MIN_COUNTS <- 50        # Minimum UMI counts per spot
MAX_MT_PERCENT <- 20    # Maximum mitochondrial gene percentage

# Set species:
SPECIES <- "mouse"      # or "human"
```

### 3. Run the Analysis

#### Option A: Interactive R Session
```r
source("run_visium_hd_analysis.R")
run_analysis()
```

#### Option B: Command Line
```bash
Rscript run_visium_hd_analysis.R
```

#### Option C: Custom Configuration
```bash
Rscript run_visium_hd_analysis.R config_mouse_wound_healing.R
```

## Detailed Usage

### Configuration Options

The `config_visium_hd.R` file contains all customizable parameters:

#### Dataset Parameters
```r
SAMPLE_NAMES <- c("CONTROL", "TREATMENT")     # Your sample names
DATA_BASE_DIR <- "data/raw/"                  # Path to data
BIN_SIZE <- c(2)                             # Analysis resolution (μm)
```

#### Quality Control
```r
MIN_COUNTS <- 50          # Minimum UMI counts per spot
MAX_MT_PERCENT <- 20      # Maximum mitochondrial percentage
MIN_FEATURES <- 10        # Minimum features per spot
```

#### Analysis Parameters
```r
SKETCH_NCELLS <- 50000                    # Cells for sketching
CLUSTER_RESOLUTIONS <- c(0.2, 0.5, 1, 2) # Clustering resolutions
N_VARIABLE_FEATURES <- 5000               # Variable features to identify
```

#### Spatial Analysis (BANKSY)
```r
BANKSY_LAMBDA <- c(0.2, 0.8)  # Balance expression vs spatial info
BANKSY_K <- c(15, 20)         # Number of spatial neighbors
```

### Example Configurations

Generate example configs for common use cases:

```r
source("run_visium_hd_analysis.R")
create_example_configs()
```

This creates:
- `config_mouse_wound_healing.R` - Time course wound healing
- `config_human_tumor.R` - Tumor heterogeneity analysis  
- `config_development.R` - Developmental biology

### Output Structure

```
output/
├── plots/
│   ├── spatial_qc_plots.pdf
│   ├── umap_multiple_resolutions.pdf
│   ├── spatial_clusters_multiple_resolutions.pdf
│   ├── banksy_clustering_results.pdf
│   └── top_marker_features.pdf
├── tables/
│   ├── markers_resolution_0.5.csv
│   ├── cluster_composition.csv
│   └── analysis_summary.json
├── RDS-files/
│   ├── merged_object_clustered.RDS
│   └── final_analysis_object.RDS
└── reports/
    └── [PROJECT_NAME]_analysis_report.html
```

## Advanced Usage

### Custom Gene Sets

Add pathway-specific gene sets to your config:

```r
CUSTOM_GENE_SETS <- list(
  "Inflammatory" = c("Il1b", "Tnf", "Il6", "Ccl2"),
  "Proliferation" = c("Mki67", "Pcna", "Top2a", "Ccnb1"),
  "ECM" = c("Col1a1", "Col3a1", "Fn1", "Lam1")
)
```

### Multi-Image Samples

For samples with multiple capture areas:

```r
IMAGE_NAMES <- list(
  "SAMPLE1" = c("REGION_A", "REGION_B"),
  "SAMPLE2" = c("REGION_C", "REGION_D")
)
```

### Memory Management

For large datasets:

```r
# Increase memory limits
MAX_MEMORY <- 4*10^12  # 4TB limit

# Reduce sketching size
SKETCH_NCELLS <- 25000

# Use fewer cores to reduce memory per process
N_CORES <- 4
```

### Batch Processing

Process multiple projects:

```bash
# Process different configurations
Rscript run_visium_hd_analysis.R config_wound_healing.R
Rscript run_visium_hd_analysis.R config_tumor.R
Rscript run_visium_hd_analysis.R config_development.R
```

## Workflow Steps

The analysis workflow includes:

1. **Data Loading**: Load Visium HD data from 10X format
2. **Quality Control**: Calculate QC metrics and filter low-quality spots
3. **Normalization**: Normalize expression data
4. **Integration**: Merge multiple samples if applicable
5. **Sketching**: Subsample for efficient analysis of large datasets
6. **Clustering**: Standard and spatial-aware clustering (BANKSY)
7. **Visualization**: Generate comprehensive plots and reports
8. **Marker Analysis**: Identify cluster-specific genes
9. **Export**: Save results and generate summary report

## Dependencies

### Required R Packages
- `Seurat` (≥4.0)
- `ggplot2`
- `patchwork` 
- `dplyr`
- `tidyr`
- `RColorBrewer`
- `ggpubr`
- `jsonlite`
- `knitr`
- `rmarkdown`

### Optional Packages
- `Banksy` (for spatial clustering)
- `SeuratWrappers`

### Installation

```r
# Install required packages
install.packages(c("Seurat", "ggplot2", "patchwork", "dplyr", 
                   "tidyr", "RColorBrewer", "ggpubr", "jsonlite",
                   "knitr", "rmarkdown"))

# Install BANKSY for spatial analysis
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Banksy")

# Install SeuratWrappers
remotes::install_github("satijalab/seurat-wrappers")
```

## Troubleshooting

### Common Issues

1. **Memory Errors**
   ```r
   # Reduce sketch size or increase memory limit
   SKETCH_NCELLS <- 25000
   MAX_MEMORY <- 4*10^12
   ```

2. **Missing Data Files**
   ```r
   # Check dataset integrity
   source("run_visium_hd_analysis.R")
   check_dataset()
   ```

3. **Species-Specific Issues**
   ```r
   # Ensure correct species setting
   SPECIES <- "mouse"  # or "human"
   ```

4. **BANKSY Installation Issues**
   - BANKSY is optional - analysis will continue without it
   - Install using BiocManager as shown above

### Performance Tips

1. **Large Datasets**: Use sketching and reduce resolution count
2. **Memory Optimization**: Process samples individually then merge
3. **Speed**: Reduce number of variable features and PCA dimensions
4. **Storage**: Use compressed RDS files and limit intermediate saves

## FAQ / Lessons Learned

Non-obvious HPC, deconvolution, and statistical gotchas from real Visium HD runs are
collected in **[FAQ.md](FAQ.md)** — e.g. single-core RCTD segfaults (BLAS threads),
tissue-specific reference limitations, why an n=2 vs n=2 comparison can't reach
significance, annotation-label hygiene, and Seurat v5 pitfalls.

## Customization

### Adding New Analysis Steps

1. Modify `General_VisiumHD_Analysis_Workflow.Rmd`
2. Add parameters to `config_visium_hd.R`
3. Update validation in `run_visium_hd_analysis.R`

### Custom Visualization

Add custom plotting functions to the workflow or create separate scripts that load the saved RDS objects.

### Integration with Other Tools

The pipeline saves standard Seurat objects that can be used with:
- Cell2location
- SPOTlight  
- RCTD
- BayesSpace
- Other spatial analysis tools

## Citation

If you use this pipeline, please cite:
- The original Seurat paper
- BANKSY if used for spatial analysis
- Any other methods incorporated in your analysis

## Contributing

This pipeline was developed for the Naik Lab. For questions or contributions:
1. Check existing issues and documentation
2. Submit detailed bug reports or feature requests
3. Follow coding standards for any contributions

## License

This pipeline is provided as-is for research use. Please respect the licenses of all incorporated tools and methods.