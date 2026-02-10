# SpatialDryArtifacts

**Adaptive edge artifact detection for 10x Visium and Visium HD.**

## Overview

**SpatialDryArtifacts** implements enhanced methods for detecting edge artifacts and dry spots—technical artifacts caused by incomplete reagent coverage or tissue handling issues. This package extends existing spatial transcriptomics quality control methods by incorporating **spatial neighborhood information** and **morphological vision**.

## Key Features

- **Dual-Platform Support**: Works on both standard **10x Visium** and high-resolution **Visium HD**.
- **Morphological Detection**: Uses raster-based focal transformations (fill, outline, star-pattern) to intelligently identify artifact clusters.
- **Hierarchical Classification**: Categorizes artifacts into actionable groups (e.g., *Large Edge Artifact*, *Small Interior Artifact*).
- **Fast & Efficient**: Optimized with `terra` for handling large datasets.

## Installation

You can install the development version from [GitHub](https://github.com/CambridgeCat13/SpatialArtifacts_v2) with:

```r
# install.packages("devtools")
devtools::install_github("CambridgeCat13/SpatialArtifacts_v2")
```

## Quick Start

```r
library(SpatialDryArtifacts)
library(SpatialExperiment)

# 1. Detect artifacts
# Option A: Standard Visium (Hexagonal Grid)
# Just specify platform = "visium"
spe <- detectEdgeArtifacts(spe, platform = "visium", qc_metric = "sum_umi")

# Option B: Visium HD (Square Grid)
# Specify platform = "visiumhd" AND the required resolution ("8um" or "16um")
# spe <- detectEdgeArtifacts(spe, platform = "visiumhd", resolution = "8um")

# 2. Classify results (Platform independent)
# Note: For Visium HD, remember to increase min_spots (e.g., min_spots = 400)
spe <- classifyEdgeArtifacts(spe)

# 3. View classification
table(spe$edge_artifact_classification)
```

## Documentation

Visit the package website: **[https://CambridgeCat13.github.io/SpatialArtifacts_v2/](https://CambridgeCat13.github.io/SpatialArtifacts_v2/)**

## Contributors

- **Harriet Jiali He**
- **Jacqueline Thompson**
- **Michael Totty**
- **Stephanie Hicks**
