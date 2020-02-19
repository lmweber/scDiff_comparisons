#######################################
# Example script demonstrating 'scDiff'
# Lukas Weber, Feb 2020
#######################################


# ------------
# installation
# ------------

# install from source (note: requires compiler)
#install.packages("scDiff_0.99.tar.gz", repos = NULL, type = "source")


# -----------
# run example
# -----------

# example demonstrating differential testing using scDiff with dataset from
# HDCytoData package

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
rd <- data.frame()
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

ad_transf <- asinh(assay(d_SE_sub) / 5)
assay(d_SE_sub) <- ad_transf

summary(assay(d_SE_sub))


# run scDiff

permute_samples <- test_CYTOF(
  d_SE_sub, 
  perm_cells = FALSE, 
  logarithm = FALSE, 
  P = 10^3, 
  N_breaks = 10, 
  min_non_zero_cells_per_group = 0
)


# output
dim(permute_samples)
head(permute_samples)


