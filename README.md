# GenoEnrich

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![R](https://img.shields.io/badge/R-≥4.0-blue.svg)](https://cran.r-project.org/)

**Genomic Region Enrichment Analysis Toolkit**

`GenoEnrich` is an R package for statistical enrichment analysis of genomic regions (e.g., NADs, LADs, enhancers, promoters) against functional annotations (e.g., chromatin states, chromatin compartmentalization). It quantifies enrichment/depletion using statistical tests and provides visualization tools.

---

## Installation

```r
if (!require("devtools")) install.packages("devtools")
devtools::install_github("yourusername/GenoEnrich")
```

> Automatically installs dependencies: `GenomicRanges`, `IRanges`, `dplyr`, `ggplot2`, `tidyr`, `scales`.

---

## Methodology

- **Statistical test**: Hypergeometric test (equivalent to one-sided Fisher’s exact test)
- **Enrichment metrics**:
  - **Enrichment Ratio** = Observed overlap / Expected overlap
  - **Odds Ratio (OR)** = (Overlap × Background non-overlap) / (Target non-overlap × State non-overlap)
- **Significance**: Raw p-values with **FDR correction** (Benjamini-Hochberg, default)
- **Background**: Customizable (default = total length of all annotated regions)
- **Input**: 
  - Target regions: `BED` file(s)
  - Feature annotations: `BED` file with a `state` column

---

## Functions

| Function | Purpose |
|----------|---------|
| `read_bed_as_gr()` | Read 3-column BED files into `GRanges` |
| `regionEnrich()` | Single-region enrichment analysis |
| `batch_regionEnrich()` | Batch analysis for multiple target regions |
| `plotEnrichmentBubble()` | Bubble plot: enrichment ratio (size), OR (color), significance (stars) |
| `plotEnrichmentBar()` | Bar plot: total state size vs. target overlap (grouped/stacked) |

---

## Quick Example

```r
library(GenoEnrich)

# Load data
target_gr <- read_bed_as_gr("targetRegion.bed")
state_gr <- read_bed_as_gr("stateRegion.bed")

# Enrichment analysis
res <- regionEnrich(target_gr, state_gr)

# Visualization
p1 <- plotEnrichmentBubble(res, title = "Target Genomic Region Enrichment")
p2 <- plotEnrichmentBar(res, title = "Target Genomic Region Coverage")
```

---

## Citation

If you use `GenoEnrich` in your research, please cite:
```
Haozhe Zhu (2025). GenoEnrich: Genomic Region Enrichment Analysis Toolkit.  
GitHub repository. https://github.com/BennyHz/GenoEnrich
```

---
