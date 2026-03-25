required_cran <- c(
  "data.table",
  "qqman",
  "ggplot2",
  "dplyr",
  "stringr",
  "coloc",
  "tidyr"
)

required_bioc <- c(
  "DESeq2",
  "EnhancedVolcano"
)

install_if_missing <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
  }
}

install_bioc_if_missing <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg, ask = FALSE, update = FALSE)
    }
  }
}

install_if_missing(required_cran)
install_bioc_if_missing(required_bioc)

message("All required packages are installed.")
