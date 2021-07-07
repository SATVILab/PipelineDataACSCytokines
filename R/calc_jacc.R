#' @title Calculate cluster stability for each cluster
#' 
#' @param fs flowSet. Assumed appropriately transformed (but not necessary scaled). 
#' @param fsom_orig output from FlowSOM::FlowSOM, as applied to original data. 
#' @param boot integer. Number of bootstrap samples to take. 
#' @param stim_and_prog_equal logical. If \code{TRUE}, then for each stimulation
#' and progressor combination an equal number of cells are sampled as in the original dataset. 
#' @param sample_by 'cell and sample'. If \code{'cell and sample'}, we randomly select a sample with replacement, 
#' and then randomly select cells from each sample with replacement, until we have
#' the same number of cells as in the original dataset. 
#' @param scale logical. If \code{TRUE}, the variables are centred and scaled. 
#' @param seed integer. If not \code{NULL}, then the seed is set to this. 
#' @param colsToUse,nClus same meaning as in FlowSOM::FlowSOM.
#' 
#' @return A tibble::tibble with a column for each cluster, where the i-th
#' entry for a particular column is the jaccard coefficient for that cluster for the i-th bootstrapped run. 

calc_jacc <- function(fs, fsom_orig, scale,
                      seed = NULL, boot, 
                      colsToUse, nClus, 
                      stim_and_prog_equal = TRUE, 
                      sample_by = 'cell and sample'){
  
  fs_input_orig <-  FlowSOM::ReadInput(input = fs, pattern = ".fcs", 
                                       compensate = FALSE, 
                                       spillover = NULL, 
                                       transform = FALSE, 
                                       toTransform = NULL, 
                                       transformFunction = flowCore::logicleTransform(),
                                       scale = scale, 
                                       silent = TRUE)
  
  jacc_tbl_out <- purrr::map_df(1:boot, function(boot_ind){
    
    if(!is.null(seed)) set.seed(seed)
    
    if(boot_ind %% 10 == 0) print(boot_ind)
    
    # get bootstrap sample
    fs_input_boot <- fs_input_orig
    if(stim_and_prog_equal){
      stim_prog_tbl_by_ind <- purrr::map_df(seq_along(fsom_orig$FlowSOM$metaData), function(i){
        fn <- names(fsom_orig$FlowSOM$metaData)[i]
        fcs_loc_end <- str_locate(fn, "/fcs/")[1,"end"]
        fn_short <- str_sub(fn, fcs_loc_end + 1)
        last_us_loc_tbl <- str_locate_all(fn_short, "_")[[1]]
        last_us_loc <- last_us_loc_tbl[nrow(last_us_loc_tbl), "end"]
        stim <- str_sub(fn_short, last_us_loc + 1, -5)
        prog <- ifelse(str_detect(fn_short, "Progressor"), "Progressor", 
                       "Control")
        tibble::tibble(stim = stim, prog = prog, 
               start_ind = fsom_orig$FlowSOM$metaData[[i]][1], 
               end_ind = fsom_orig$FlowSOM$metaData[[i]][2], 
               n_cell = end_ind - start_ind + 1)
      })
      if(sample_by == 'cell and sample'){
        prog_vec <- unique(stim_prog_tbl_by_ind$prog)
        stim_vec <- unique(stim_prog_tbl_by_ind$stim)
        ind_vec <- purrr::map(prog_vec, function(prog){
          purrr::map(stim_vec, function(stim_sample){
            sample_tbl <- stim_prog_tbl_by_ind %>%
              dplyr::filter(prog == .env$prog, 
                     stim == stim_sample)
            n_cell_total <- sum(sample_tbl$n_cell)
            n_sample_curr <- nrow(sample_tbl)
            sample_vec_curr <- NULL
            while(length(sample_vec_curr) < n_cell_total){
              sample_ind <- sample(1:n_sample_curr, 1)
              ind_vec_sample_orig <- sample_tbl[sample_ind,][["start_ind"]]:sample_tbl[sample_ind,][["end_ind"]]
              n_max_sample <- min(n_cell_total - length(sample_vec_curr), 
                                  diff(range(ind_vec_sample_orig)) + 1)
              sample_vec_curr <- c(sample_vec_curr, 
                                   sample(ind_vec_sample_orig, 
                                          n_max_sample, 
                                          replace = TRUE))
              
            }
            sample_vec_curr
          })
        }) %>%
          unlist()
      } else if(is.null(sample_by)){
        ind_vec <- sample(1:nrow(fs_input_orig$data),
                          nrow(fs_input_orig$data), 
                          replace = TRUE)
      }
    }
    
    fs_input_boot$data <- fs_input_orig$data[ind_vec,]
    
    # build SOM
    fsom <- FlowSOM::BuildSOM(fsom = fs_input_boot, 
                              colsToUse = colsToUse, 
                              silent = TRUE, 
                              importance = NULL)
    
    # build mst
    fsom <- FlowSOM::BuildMST(fsom, silent = TRUE)
    
    # metaclustering
    cluster_lab_vec <- FlowSOM::metaClustering_consensus(fsom$map$codes,
                                                         nClus, 
                                                         seed = 2 * seed)
    
    # clusters for each observations
    cluster_vec <- cluster_lab_vec[fsom$map$mapping[,1]]
    
    # clusters for each original data point
    cluster_vec_orig <- fsom_orig$metaclustering[fsom_orig$FlowSOM$map$mapping[,1]]
    
    # clusters for each original data point that is also now clustered
    ind_vec_uni <- unique(ind_vec)
    # this says what each original cell (whose index is in ind_vec_uni)
    # was mapped to in the origin clustering
    cluster_vec_orig_sample <- as.character(cluster_vec_orig[ind_vec_uni])
    # this says what each original cell (whose index is in ind_vec_uni)
    # was mapped to in the new clustering
    cluster_vec_matched_to_orig <- rep(NA_character_, length(ind_vec_uni))
    for(i in seq_along(ind_vec_uni)){
      cluster_vec_matched_to_orig[i] <- cluster_vec[ind_vec == ind_vec_uni[i]][1]
    }
    
    # jaccard coefficients for clusters in bootstrapped dataset 
    clusters_in_boot_sample_vec <- unique(cluster_vec_orig_sample)
    jc_vec <- purrr::map(clusters_in_boot_sample_vec, function(clust_orig){
      
      jacc_vec_boot <- purrr::map_dbl(1:nClus, function(clust_boot){
        ind_boot <- which(cluster_vec_matched_to_orig == clust_boot)
        ind_orig <- which(cluster_vec_orig_sample == clust_orig)
        length(intersect(ind_boot, ind_orig))/length(union(ind_boot, ind_orig))
      })
      #print(clusters_in_boot_sample_vec[which(jacc_vec_boot == max(jacc_vec_boot))[1]])
      max(jacc_vec_boot)
    }) %>%
      setNames(clusters_in_boot_sample_vec)
    
    # jaccard coefficient set to -1 for clusters not in bootstrapped
    # dataset
    clusters_not_in_boot_sample_vec <- setdiff(as.character(1:nClus), 
                                               clusters_in_boot_sample_vec)
    
    for(i in seq_along(clusters_not_in_boot_sample_vec)){
      jc_vec <- c(jc_vec, 
                  setNames(-1, clusters_not_in_boot_sample_vec[i]))
    }
    
    seed <<- seed + 1
    
    out_tbl <- data.frame(jc_vec)
    colnames(out_tbl) <- names(jc_vec)
    out_tbl[,as.character(1:nClus)] %>%
      dplyr::mutate(boot = boot_ind) %>%
      dplyr::select(boot, everything())
    
  })
  
}