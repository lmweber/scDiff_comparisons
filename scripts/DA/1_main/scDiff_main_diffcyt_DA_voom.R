######################################
# scDiff comparisons: main simulations
# Lukas Weber, Apr 2020
######################################

# Script to run diffcyt methods on main simulations using datasets from
# HDCytoData package (Weber_AML_sim)


library(diffcyt)
library(SummarizedExperiment)
library(HDCytoData)



###############
# Load datasets
###############

# load data from HDCytoData package
data_main <- Weber_AML_sim_main_1pc_SE()

dim(data_main)


# subset object to keep CN and healthy conditions
data_sub <- data_main[rowData(data_main)$group_id %in% c("healthy", "CN"), ]

# drop empty levels in row data
rowData(data_sub) <- droplevels(rowData(data_sub))
# drop empty levels in experiment_info
experiment_info <- metadata(data_sub)$experiment_info
experiment_info <- droplevels(experiment_info[experiment_info$group_id %in% c("healthy", "CN"), ])
metadata(data_sub)$experiment_info <- experiment_info
# drop empty levels in n_cells
keep <- grepl("healthy|CN", names(metadata(data_sub)$n_cells))
metadata(data_sub)$n_cells <- metadata(data_sub)$n_cells[keep]

dim(data_sub)



#####################
# Run diffcyt-DA-voom
#####################


##################
# diffcyt pipeline
##################

# --------------------
# pre-processing steps
# --------------------

runtime_preprocessing <- system.time({
  
  # get dataset
  d_se <- data_sub
  
  colnames(d_se)[colData(d_se)$marker_class == "type"]
  colnames(d_se)[colData(d_se)$marker_class == "state"]
  
  # transform data
  d_se <- transformData(d_se, cofactor = 5)
  
  # clustering
  # (runtime: ~15 sec with xdim = 20, ydim = 20)
  seed <- 123
  d_se <- generateClusters(d_se, xdim = 20, ydim = 20, seed_clustering = seed)
  
  length(table(rowData(d_se)$cluster_id))  # number of clusters
  nrow(rowData(d_se))                      # number of cells
  sum(table(rowData(d_se)$cluster_id))
  min(table(rowData(d_se)$cluster_id))     # size of smallest cluster
  max(table(rowData(d_se)$cluster_id))     # size of largest cluster
  
  # calculate cluster cell counts
  d_counts <- calcCounts(d_se)
  
  dim(d_counts)
  rowData(d_counts)
  length(assays(d_counts))
  
  # calculate cluster medians
  # (note: requires updating column names due to change in formatting of
  # column names in HDCytoData package vs. diffcyt paper)
  colnames(d_se) <- gsub("\\(.*$", "", colnames(d_se))
  
  d_medians <- calcMedians(d_se)
  
  dim(d_medians)
  rowData(d_medians)
  length(assays(d_medians))
  names(assays(d_medians))
  
  # calculate medians by cluster and marker
  d_medians_by_cluster_marker <- calcMediansByClusterMarker(d_se)
  
  dim(d_medians_by_cluster_marker)
  length(assays(d_medians_by_cluster_marker))
  
  # calculate medians by sample and marker
  d_medians_by_sample_marker <- calcMediansBySampleMarker(d_se)
  
  dim(d_medians_by_sample_marker)
  length(assays(d_medians_by_sample_marker))
  
})


# ---------------------------------
# store data objects (for plotting)
# ---------------------------------

out_objects_diffcyt_DA_voom_main <- list(
  d_se = d_se, 
  d_counts = d_counts, 
  d_medians = d_medians, 
  d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
  d_medians_by_sample_marker = d_medians_by_sample_marker
)


# -----------------------------------------
# test for differentially abundant clusters
# -----------------------------------------

# contrast (to test 'CN' vs. 'healthy')
# note: include fixed effects for 'patient_id'
contrast_vec <- c(0, 1, 0, 0, 0, 0)

runtime_tests <- system.time({
  
  # set up design matrix
  # note: include fixed effects for 'patient_id'
  design <- createDesignMatrix(metadata(d_se)$experiment_info, 
                               cols_design = c("group_id", "patient_id"))
  design
  
  # set up contrast matrix
  contrast <- createContrast(contrast_vec)
  contrast
  
  # run tests
  res <- testDA_voom(d_counts, design, contrast, plot = FALSE)
  
})

# show results
rowData(res)

# sort to show top (most highly significant) clusters first
res_sorted <- rowData(res)[order(rowData(res)$p_adj), ]
print(head(res_sorted, 10))

# number of significant tests (note: one test per cluster)
print(table(res_sorted$p_adj <= 0.1))

# formatted summary table
topTable(res)

# runtime
runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_tests[["elapsed"]]
print(runtime_total)

runtime_diffcyt_DA_voom_main <- runtime_total


# ---------------------------------------------
# store results at cluster level (for plotting)
# ---------------------------------------------

res_clusters <- as.data.frame(rowData(res))

out_clusters_diffcyt_DA_voom_main <- res_clusters



##############################
# Return results at cell level
##############################

# Note: diffcyt methods return results for each cluster. To enable performance
# comparisons between methods at the cell level, we assign the cluster-level p-values
# to all cells within each cluster.


# spike-in status for each cell
is_spikein <- rowData(data_sub)$spikein


# match cluster-level p-values to individual cells

stopifnot(nrow(rowData(res)) == length(levels(rowData(d_se)$cluster_id)))
stopifnot(all(rowData(res)$cluster_id == levels(rowData(d_se)$cluster_id)))

rowData(res)$cluster_id <- factor(rowData(res)$cluster_id, levels = levels(rowData(d_se)$cluster_id))

# match cells to clusters
ix_match <- match(rowData(d_se)$cluster_id, rowData(res)$cluster_id)

p_vals_clusters <- rowData(res)$p_val
p_adj_clusters <- rowData(res)$p_adj

p_vals_cells <- p_vals_clusters[ix_match]
p_adj_cells <- p_adj_clusters[ix_match]


# set up data frame with results and true spike-in status at cell level

stopifnot(length(p_vals_cells) == length(is_spikein))
stopifnot(length(p_adj_cells) == length(is_spikein))

res_p_vals <- p_vals_cells
res_p_adj <- p_adj_cells

# replace NAs (due to filtering) to ensure same cells are returned for all methods
res_p_vals[is.na(res_p_vals)] <- 1
res_p_adj[is.na(res_p_adj)] <- 1

# return values for this condition and healthy
res <- data.frame(p_val = res_p_vals, 
                  p_adj = res_p_adj, 
                  spikein = is_spikein)

# store results
out_diffcyt_DA_voom_main <- res



################
# Output objects
################

# names of output objects from this script (to save or use in next script)

#out_diffcyt_DA_voom_main
#runtime_diffcyt_DA_voom_main
#out_clusters_diffcyt_DA_voom_main
#out_objects_diffcyt_DA_voom_main



