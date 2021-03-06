######################################
# scDiff comparisons: null simulations
# Lukas Weber, Mar 2020
######################################

# Script to run diffcyt methods on null simulations using datasets from
# HDCytoData package (Weber_BCR_XL_sim)


library(diffcyt)
library(SummarizedExperiment)
library(HDCytoData)



###############
# Load datasets
###############

# Load data from HDCytoData package

data_null_sims <- list(
  rep1 = Weber_BCR_XL_sim_null_rep1_SE(), 
  rep2 = Weber_BCR_XL_sim_null_rep2_SE(), 
  rep3 = Weber_BCR_XL_sim_null_rep3_SE()
)



######################
# Run diffcyt-DS-limma
######################

# names of replicates
rep_names <- names(data_null_sims)

# lists to store objects
out_diffcyt_DS_limma_null  <- 
  out_clusters_diffcyt_DS_limma_null <- 
  out_objects_diffcyt_DS_limma_null <- 
  runtime_diffcyt_DS_limma_null <- vector("list", length(rep_names))
names(out_diffcyt_DS_limma_null) <- 
  names(out_clusters_diffcyt_DS_limma_null) <- 
  names(out_objects_diffcyt_DS_limma_null) <- 
  names(runtime_diffcyt_DS_limma_null) <- rep_names


# run for each replicate

for (s in 1:length(rep_names)) {
  
  
  ##################
  # diffcyt pipeline
  ##################
  
  # --------------------
  # pre-processing steps
  # --------------------
  
  runtime_preprocessing <- system.time({
    
    # get dataset
    d_se <- data_null_sims[[s]]
    
    colnames(d_se)[colData(d_se)$marker_class == "type"]
    colnames(d_se)[colData(d_se)$marker_class == "state"]
    
    # transform data
    d_se <- transformData(d_se, cofactor = 5)
    
    # clustering
    # (runtime: ~5 sec with xdim = 10, ydim = 10)
    seed <- 1234
    d_se <- generateClusters(d_se, xdim = 10, ydim = 10, seed_clustering = seed)
    
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
  
  out_objects_diffcyt_DS_limma_null[[s]] <- list(
    d_se = d_se, 
    d_counts = d_counts, 
    d_medians = d_medians, 
    d_medians_by_cluster_marker = d_medians_by_cluster_marker, 
    d_medians_by_sample_marker = d_medians_by_sample_marker
  )
  
  
  # --------------------------------------------
  # test for differential states within clusters
  # --------------------------------------------
  
  # contrast (to compare 'null2' vs. 'null1')
  # note: include fixed effects for 'patient_id'
  contrast_vec <- c(0, 1, 0, 0, 0, 0, 0, 0, 0)
  
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
    res <- testDS_limma(d_counts, d_medians, design, contrast, plot = FALSE)
    
  })
  
  # show results
  rowData(res)
  
  # check p-value distribution for null simulation
  # (note: one test per cluster-marker combination)
  hist(rowData(res)$p_val)
  
  # runtime
  runtime_total <- runtime_preprocessing[["elapsed"]] + runtime_tests[["elapsed"]]
  print(runtime_total)
  
  runtime_diffcyt_DS_limma_null[[s]] <- runtime_total
  
  
  # ---------------------------------------------
  # store results at cluster level (for plotting)
  # ---------------------------------------------
  
  res_clusters <- as.data.frame(rowData(res))
  
  out_clusters_diffcyt_DS_limma_null[[s]] <- res_clusters
  
  
  
  ##############################
  # Return results at cell level
  ##############################
  
  # Note: diffcyt methods return results for each cluster-marker combination. To enable
  # performance comparisons between methods at the cell level, we assign the same p-values
  # to all cells within a given cluster-marker combination.
  
  # Note: return cell-level results for marker pS6 only, since the comparative evaluations
  # (ROC curves etc) are based on pS6 in B cells only.
  
  
  # identify B cells (these contain the true differential signal)
  is_B_cell <- rowData(d_se)$B_cell
  
  
  # match cluster-level p-values for marker pS6 to individual cells
  
  n_state_markers <- length(colnames(d_se)[colData(d_se)$marker_class == "state"])
  stopifnot(nrow(rowData(res)) == nlevels(rowData(d_se)$cluster_id) * n_state_markers)
  stopifnot(all(levels(rowData(res)$cluster_id) == levels(rowData(d_se)$cluster_id)))
  stopifnot(all(levels(rowData(res)$cluster_id) %in% rowData(res)$cluster_id))
  
  # select results for pS6
  res_pS6 <- res[rowData(res)$marker_id == "pS6", ]
  
  # match cells to clusters
  ix_match <- match(rowData(d_se)$cluster_id, rowData(res_pS6)$cluster_id)
  
  p_vals_clusters <- rowData(res_pS6)$p_val
  p_adj_clusters <- rowData(res_pS6)$p_adj
  
  p_vals_cells <- p_vals_clusters[ix_match]
  p_adj_cells <- p_adj_clusters[ix_match]
  
  
  # set up data frame with results (for marker pS6) and true B-cell status at cell level
  
  res_p_vals <- p_vals_cells
  res_p_adj <- p_adj_cells
  
  # replace NAs (due to filtering) to ensure same cells are returned for all methods
  res_p_vals[is.na(res_p_vals)] <- 1
  res_p_adj[is.na(res_p_adj)] <- 1
  
  stopifnot(length(res_p_vals) == length(res_p_adj))
  stopifnot(length(res_p_vals) == length(is_B_cell))
  
  res <- data.frame(p_val = res_p_vals, 
                    p_adj = res_p_adj, 
                    B_cell = is_B_cell)
  
  # store results
  out_diffcyt_DS_limma_null[[s]] <- res
  
}



################
# Output objects
################

# names of output objects from this script (to save or use in next script)

#out_diffcyt_DS_limma_null
#runtime_diffcyt_DS_limma_null
#out_clusters_diffcyt_DS_limma_null
#out_objects_diffcyt_DS_limma_null



