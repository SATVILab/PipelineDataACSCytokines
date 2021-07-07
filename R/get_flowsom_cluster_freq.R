#' @export
get_fs_cluster_freq <- function(flowsom_out, uns_chr = 'uns', dir_save_cluster_results = NULL, stats_combn_tbl, stim){
  
  # ------------------
  # Get total classified cells per sample and stim combination
  # ------------------
  
  
  # extract the total classified cells for each sample for stim (and unstim)
  # --------------------
  
  # subset to one row per sample and stim combn
  cell_count_tbl <- stats_combn_tbl %>%
    dplyr::group_by(SubjectID, VisitType, stim) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(n_cell = n_cell_stim)
  
  # extract uns
  cell_count_tbl_uns <- stats_combn_tbl %>% 
    dplyr::group_by(SubjectID, VisitType) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::rename(n_cell = n_cell_uns) %>%
    dplyr::mutate(stim = uns_chr)
  
  # bind stim and uns sample counts together
  cell_count_tbl %<>%
    dplyr::bind_rows(cell_count_tbl_uns) %>%
    dplyr::select(SubjectID, VisitType, stim, n_cell) %>%
    dplyr::mutate(sampleid_stim = paste0(SubjectID, "_", VisitType, "_", stim))
  
  # get vector id for sampleid and stim combination
  cell_count_lab_vec <- setNames(cell_count_tbl$n_cell, cell_count_tbl$sampleid_stim)
  
  # add SubjectID_VisitType_stim column to flowsom_out
  flowsom_out %<>%
    dplyr::mutate(SampleID_Stim = paste0(SubjectID, "_", VisitType, "_", stim))
  
  # ------------------
  # Calculate cluster counts and frequencies per sample
  # ------------------
  
  # calculate frequency for each sample and stim combination
  # -----------------------
  
  flowsom_out_freq <- flowsom_out %>%
    dplyr::group_by(fcs, SubjectID, VisitType, Progressor, timeToTB, stim) %>%
    dplyr::mutate(n_cell_cyt_stim = n(), 
           sampleid_stim = paste0(SubjectID, "_", VisitType, "_", stim),
           n_cell_tot_stim = cell_count_lab_vec[sampleid_stim]) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(fcs, SubjectID, VisitType, Progressor, timeToTB, stim, n_cell_tot_stim, n_cell_cyt_stim, clust, stability) %>%
    dplyr::summarise(count_stim = n(), .groups = 'drop') %>%
    dplyr::mutate(freq_cyt_stim = count_stim/n_cell_cyt_stim * 1e2, freq_tot_stim = count_stim/n_cell_tot_stim*1e2)
  
  # add zero frequency for samples that responded for that stimulation but 
  # did not actually have cells that belong to a specific cluster
  # -----------------------
  
  stability_tbl <- flowsom_out_freq %>% 
    dplyr::group_by(clust) %>%
    dplyr::slice(1) %>%
    dplyr::select(clust, stability) %>%
    dplyr::ungroup()
  stability_lab_vec <- setNames(stability_tbl$stability, stability_tbl$clust)
  
  # loop over stimulations in flowsom_out_freq
  # to ensure that only people who responded to a particular stimulation
  # get a zero frequency for that stimulation
  for(stim_ind in unique(flowsom_out_freq$stim)){
    flowsom_out_freq_stim <- flowsom_out_freq %>% dplyr::filter(stim == stim_ind)
    for(subjectid in unique(flowsom_out_freq$SubjectID)){
      flowsom_out_freq_subjectid <- flowsom_out_freq_stim %>% 
        dplyr::filter(SubjectID == subjectid)
      visittype_vec <- unique(flowsom_out_freq_subjectid$VisitType)
      for(visittype in visittype_vec){
        flowsom_out_freq_subjectid_visittype <- flowsom_out_freq_subjectid %>%
          dplyr::filter(VisitType == visittype)
        for(clust in unique(flowsom_out_freq$clust)){
          # if that cluster not found for this sample, 
          # then add a zero freq for it as clearly it's zero
          # as the sample was a responder for this stim (if filtering 
          # was done by responder status)
          if(!clust %in% flowsom_out_freq_subjectid_visittype$clust){
            flowsom_out_freq %<>% 
              dplyr::bind_rows(flowsom_out_freq_subjectid_visittype[1,] %>%
                          dplyr::mutate(count_stim = 0, freq_cyt_stim = 0, freq_tot_stim = 0, clust = .env$clust, 
                                 stability = stability_lab_vec[clust])) # surely you should then say clust =- clust?
          }
        }
      }
    }
  }
  
  # calculate freq_tot_bs and freq_cyt_bs
  if('uns' %in% unique(flowsom_out_freq$stim)){
    
    flowsom_out_freq_stim <- flowsom_out_freq 
    
    flowsom_out_freq_uns <- flowsom_out_freq %>% 
      dplyr::filter(stim == 'uns') %>%
      dplyr::rename(count_uns = count_stim, 
             n_cell_tot_uns = n_cell_tot_stim, 
             n_cell_cyt_uns = n_cell_cyt_stim,
             freq_cyt_uns = freq_cyt_stim, 
             freq_tot_uns = freq_tot_stim) %>%
      dplyr::select(SubjectID, VisitType, clust, n_cell_tot_uns, n_cell_cyt_uns, count_uns, freq_tot_uns, freq_cyt_uns)
    
    flowsom_out_freq <- flowsom_out_freq_stim %>%
      dplyr::left_join(flowsom_out_freq_uns, by = c("SubjectID", "VisitType", "clust")) %>%
      dplyr::mutate(freq_cyt_uns = ifelse(is.na(freq_cyt_uns), 0, freq_cyt_uns), 
             freq_tot_uns = ifelse(is.na(freq_tot_uns), 0, freq_tot_uns)) %>%
      dplyr::mutate(freq_cyt_bs = freq_cyt_stim - freq_cyt_uns, 
             freq_tot_bs = freq_tot_stim - freq_tot_uns) 
    
  } else{
    flowsom_out_freq %<>%
      dplyr::mutate(n_cell_tot_uns = NA, 
             n_cell_cyt_uns = NA, 
             count_uns = NA, 
             freq_cyt_uns = NA, 
             freq_tot_uns = NA,
             freq_cyt_bs = NA, 
             freq_tot_bs = NA)
  }
  
  # make cluster and stim variables factors
  # -----------------------
  
  # convert mtb to mtbaux if mtbaux is what's in flowsom_out
  if('mtbaux' %in% flowsom_out_freq$stim){
    flowsom_out_freq %<>% dplyr::mutate(stim = ifelse(stim == 'mtbaux', 'mtb', stim)) 
  }
  
  # get order vector (assuming unstim is available)
  stim_order_vec <- switch(stim, 
                           "all_u" = c("uns", "ebv", "mtb", "p1"), 
                           "all" = c("ebv", "mtb", "p1"), 
                           c('uns', stim))
  
  # remove unstim if it is available
  if(!'uns' %in% flowsom_out_freq$stim) stim_order_vec %<>% setdiff('uns')
  
  # make stim and cluster factors 
  stim_factor_vec <- factor(flowsom_out_freq$stim, stim_order_vec)
  
  flowsom_out_freq %<>%
    dplyr::mutate(stim = stim_factor_vec,
           clust = factor(as.character(clust), levels = as.character(sort(as.numeric(unique(clust))))))
  
  
  # arrange
  flowsom_out_freq %<>% dplyr::arrange(SubjectID, VisitType, stim, clust)
  
  if(!is.null(dir_save_cluster_results)){
    path_fn <- file.path(dir_save_cluster_results, "fs_freq.rds")
    if(file.exists(path_fn)) file.copy(path_fn, file.path(dir_save_cluster_results, "fs_freq_old.rds"))
    saveRDS(flowsom_out_freq, path_fn)
  }
  
  flowsom_out_freq 
    
}