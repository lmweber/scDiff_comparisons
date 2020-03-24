######################################
# Example script for 'scDiff' analyses
# Lukas Weber, Mar 2020
######################################


# ------------
# Installation
# ------------

# install from source (note: requires compiler)
#install.packages("scDiff_0.99.tar.gz", repos = NULL, type = "source")


# --------------
# Simple example
# --------------

# simple example showing how to perform differential testing using scDiff, with
# dataset from HDCytoData package

library(scDiff)
library(HDCytoData)
library(SummarizedExperiment)


# load dataset: Bodenmiller_BCR_XL

d_SE <- Bodenmiller_BCR_XL_SE()


# subset a small number of cells per sample for faster runtime

n_sub <- 100

sample_names <- levels(rowData(d_SE)$sample_id)
sample_names

# subset and put back together into original object format
rd <- DataFrame()
cd <- colData(d_SE)
md <- metadata(d_SE)
ad <- c()

set.seed(123)
for (i in seq_along(sample_names)) {
  ix_i <- seq_len(nrow(d_SE))[rowData(d_SE)$sample_id == sample_names[i]]
  subset_ix_i <- sample(ix_i, n_sub)
  rd <- rbind(rd, rowData(d_SE)[subset_ix_i, , drop = FALSE])
  ad <- rbind(ad, assay(d_SE)[subset_ix_i, , drop = FALSE])
}

d_SE_sub <- SummarizedExperiment(
  assays = list(exprs = ad), 
  rowData = rd, 
  colData = cd, 
  metadata = md
)

d_SE_sub

dim(d_SE_sub)


# transform data: using standard transform asinh(x/5) for CyTOF data (similar to
# log, but linear near zero and allows negative values)

cols_transf <- colData(d_SE_sub)$marker_class != "none"
cofactor <- 5
ad_transf <- assay(d_SE_sub)
ad_transf[, cols_transf] <- asinh(ad_transf[, cols_transf, drop = FALSE] / cofactor)
assay(d_SE_sub) <- ad_transf

summary(assay(d_SE_sub))


# run scDiff

res_scDiff <- test_CYTOF(
  d_SE_sub, 
  perm_cells = FALSE, 
  logarithm = FALSE, 
  P = 10^3, 
  N_breaks = 10, 
  min_non_zero_cells_per_group = 0
)


# output
dim(res_scDiff)
head(res_scDiff)
#View(res_scDiff)


# --------------
# Interpretation
# --------------

# Since we are testing for differential states (DS) within cell populations, for
# CyTOF data, it only makes sense to test 'state' markers. Therefore, we subset
# the final results to show 'state' markers only.

# However, this may break the p-value adjustments for multiple testing.

markers_type <- colData(d_SE)$marker_name[colData(d_SE)$marker_class == "type"]
markers_state <- colData(d_SE)$marker_name[colData(d_SE)$marker_class == "state"]

markers_type
markers_state

res_scDiff <- res_scDiff[res_scDiff$gene %in% markers_state, , drop = FALSE]
rownames(res_scDiff) <- NULL

dim(res_scDiff)
head(res_scDiff)
#View(res_scDiff)


