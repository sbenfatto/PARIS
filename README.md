---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->



# PARIS

<!-- badges: start -->
<!-- badges: end -->

PAn-canceR Inferred Synthetic lethalities (PARIS), a machine learning-based approach to predict cancer vulnerabilities.

## Installation

You can install the released version of PARIS with:

```r
devtools::install_github()
```

## Example

How to run PARIS

```r
library(PARIS)
run_PARIS(imp.score = "raw", omic = "expression", outdir = "/PARIS_test", genelist_dep = yourgenelist, genelist_feat = yourgenelist2)
```

