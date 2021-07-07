#'
#' @description Loads flowSet for flowSOM, automatically filtering out non-responders (if desired), automatically 
#' excluding files that didn't read in correctly and 
#' 
#' @param stim 'all_u', 'all', 'p1', 'p4', 'ebv' or 'mtbaux'. Specifies stimulation(s) that should be read in together 
#' as one flowSet. 
#' @param dir_fcs character. Directory containing FCS files. 
#' @param stim_vec_u character vector. Names of all stimulations that could be read in together. 
#' @param uns_chr character. Name for unstim in \code{stim_vec_u}. Default is "uns". 
#' @param responders_list named list. Named list, where names are stimulations and each element is a character vector 
#' of samples that respond. 
#' @param conv_mtbaux_to_mtb logical. If \code{TRUE}, then when matching fcs files to responder lists, 
#' "mtbaux"  is converted to "mtb". Default is \code{FALSE}.
#' @param trans object of class 'transform'. Transformation to appy to data. If \code{NULL}, then no transformation is applied.  
#' Default is \code{flowCore::arcsinhTransform(b = 1/5,a = 0)}, 
#' which corresponds to asinh(x/5). 
#' @param chnl_sel character vector. Channel names to transform, if transformation is supplied. 
#' 
#' @details If \code{stim='all_u"} or \code{stim='all'}, then all samples are loaded from each stimulated sample that are responders (even if
#' not Mtb-specific). If \code{stim = 'all'), then unstim samples are only loaded from individuals that responded 
#' to at least one of the stimulations. 
#' @export
load_fs_for_flowsom <- function(dir_fcs, stim, stim_vec_u = NULL, responders_list, chnl_sel = NULL, 
                              uns_chr = "uns", conv_mtbaux_to_mtb = FALSE, 
                              trans = flowCore::arcsinhTransform(b = 1/5,a = 0)){
  
  # =========================
  # Preparation 
  # =========================
  
  # get list of stimulations to read in
  stim_vec_load <- switch(stim, 
                          "all" = setdiff(stim_vec_u, uns_chr), 
                          "all_u" = stim_vec_u, 
                          stim)
  
  
  # all fcs files that match stim condition(s)
  # ------------------------
  
  # get list of FCS files for each stim
  fcs_list <-  purrr::map(stim_vec_load, function(stim_ind){
    c(list.files(dir_fcs, pattern = paste0(stim_ind ,".fcs")), 
      list.files(dir_fcs, pattern = paste0(stim_ind ,"aux.fcs")))
  }) %>%
    setNames(stim_vec_load)
  
  # return an error if none are found that match the pattern
  if(length(unlist(fcs_list)) == 0) stop(paste0("no FCS files found that match pattern in ", dir_fcs))
  
  
  # remove file names corresponding to non-responders (if desired)
  # -------------------------
  
  # only run if responders_list is not NULL
  if(!is.null(responders_list)){
    
    # convert mtbaux in fcs_list to mtb if pre-specified or if mtbaux is in responders_list names
    if(conv_mtbaux_to_mtb || 'mtb' %in% names(responders_list)){
      names(fcs_list)[which(names(fcs_list) == 'mtbaux')] <- 'mtb'
    }
    
    # print a warning that the non-unstim stimulations do not match if so
    if(length(setdiff(setdiff(names(fcs_list), uns_chr), 
                      setdiff(names(responders_list), uns_chr))) != 0 || 
       length(setdiff(setdiff(names(responders_list), uns_chr), 
                      setdiff(names(fcs_list), uns_chr))) != 0){
      message(paste0("names from fcs_list (", paste0(setdiff(names(fcs_list), uns_chr), collapse = " "), 
                     ") do not match names from responders_list (", paste0(setdiff(names(responders_list), uns_chr), collapse = " "), ")"))

    }
    
    # loop over stimulations for which response probs have been calculated
    for(stim_resp in names(responders_list)){
      fcs_vec_stim <- fcs_list[[stim_resp]]
      for(fcs in fcs_vec_stim){
        fcs_resp <- any( purrr::map_lgl(responders_list[[stim_resp]], function(x) str_detect(fcs, x)))
        if(!fcs_resp) fcs_vec_stim <- fcs_vec_stim[-which(fcs_vec_stim == fcs)]
      }
      fcs_list[[stim_resp]] <- fcs_vec_stim
    }
    
    # only choose unstim samples that correspond to stimulated samples that responded
    if(stim == 'all_u'){
      resp_vec <- unlist(responders_list) %>% unique()
      fcs_vec_uns <- fcs_list[[uns_chr]]
      for(fcs in fcs_vec_uns){
        fcs_resp <- any( purrr::map_lgl(resp_vec, function(x) str_detect(fcs, x)))
        if(!fcs_resp) fcs_vec_uns <- fcs_vec_uns[-which(fcs_vec_uns == fcs)]
      }
      fcs_list[[uns_chr]] <- fcs_vec_uns
    }
  }
  
  
  # create character vector of all short file names to be read in
  fcs_vec <- unlist(fcs_list)


  # only choose fcs files that read in without error
  # -------------------------
  for(fcs in fcs_vec){
    ex <- try(flowCore::read.FCS(file.path(dir_fcs, fcs)), 
              silent = TRUE)
    
    if(class(ex) == 'try-error'){
      # print(fcs)
      fcs_vec <- fcs_vec[-which(fcs_vec == fcs)]
    } 
  }
  
  # =========================
  # Read in flowSet
  # =========================
  
  fs <- flowCore::read.flowSet(files = file.path(dir_fcs, fcs_vec))
  
  # =========================
  # Transform flowSet, if specified 
  # =========================
  
  if(!is.null(trans)){
    if(is.null(chnl_sel)) chnl_sel <- as.character(attr(exprs(fs[[1]]), "dimnames")[[2]])
    trans_list <- flowCore::transformList(chnl_sel, trans)
    fs <- flowCore::transform(fs, trans_list)
  }
  
  
  fs
}