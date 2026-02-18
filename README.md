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

## Important: `shifted` Parameter Warning

**For standard Visium data loaded via `read10xVisium()` or `SpatialExperiment`: ALWAYS use `shifted = FALSE` (default).**

### Why this matters:
Standard Visium array coordinates (`array_col`: 0–127, `array_row`: 0–77) already form a **regular integer grid** where every spot has a unique (array_row, array_col) combination. The hexagonal offset pattern is already encoded in these coordinate values—**no adjustment is needed**.

### When to use `shifted = TRUE`:
`shifted = TRUE` is **only** relevant for non-standard coordinate systems where odd columns have a physical half-unit offset that is **not** already corrected in the coordinate values themselves. This is extremely rare and typically only applies to custom data formats.

### How to verify your data:
You can check whether your coordinates are already regular (and thus require `shifted = FALSE`) by running:
```r
# If this returns TRUE → use shifted = FALSE (coordinates are already regular)
nrow(unique(colData(spe)[, c("array_row", "array_col")])) == ncol(spe)
```

If this check returns `TRUE`, your coordinates are **already properly formatted** for standard Visium data, and you should **NOT** use `shifted = TRUE`.

### What happens if you incorrectly use `shifted = TRUE`:
Using `shifted = TRUE` on standard Visium data will cause **coordinate transposition errors** that corrupt the spatial topology:
- Artifact detection will fail to identify true edge regions
- Spots will be incorrectly mapped in raster space
- Classification results will be unreliable

**Bottom line**: Unless you have a custom, non-standard coordinate system (very rare), always use the default `shifted = FALSE`.

## Documentation
Visit the package website: **[https://CambridgeCat13.github.io/SpatialArtifacts_v2/](https://CambridgeCat13.github.io/SpatialArtifacts_v2/)**

## Contributors
- **Harriet Jiali He**
- **Jacqueline Thompson**
- **Michael Totty**
- **Stephanie Hicks**
