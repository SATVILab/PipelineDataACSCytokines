#' @title Get path to FlowSOM results\\
#' 
#' @description 
#' Get the path to the directory containing 
#' FlowSOM results for a given dataset 
#' under a differernt parameterisation. 
#'
#' @export
get_dir_flowsom_cyt <- function(dir_gating,
                                pop,
                                chnl,
                                exc,
                                non_responders,
                                stim,
                                n_clust){
      
      # prep
      # -------------
      
      dir_gs <- switch(
        pop,
        "cd4" = 'gs_cytof_acs_cd4',
        "cd8" = 'gs_cytof_acs_cd8',
        "tcrgd" = 'gs_cytof_acs_tcrgd',
        "nk" = 'gs_cytof_acs_nk',
        "bcell" = 'gs_cytof_acs_bcell'
      )
      dir_chnl <- paste0(chnl, collapse = "_")
      dir_exc <- purrr::map_chr(exc, function(exc_curr){
        if(is.null(exc_curr)) return("exc-none")
        paste0("exc-", paste0(exc_curr, collapse = "_&_"))
      })
      
      dir_non_responders <- switch(
        as.character(is.character(non_responders)),
        "TRUE" = switch(
          non_responders,
          "exc" = ,
          "r_o" = "r_o",
          "inc" = ,
          "all" = "all",
          stop(paste0("non_responders value of ", non_responders, " not valid"))
        ),
        "FALSE" = switch(
          as.character(non_responders),
          "TRUE" = "all",
          "FALSE" = "r_o",
          stop(paste0("non_responders value of ", non_responders, " not valid"))
        ))
      
      dir_stim <- switch(
        stim,
        "all" = ,
        "all_u-all" = 'all_u-all',
        "p1" = ,
        "all_u-p1" = 'all_u-p1',
        "ebv" = ,
        "all_u-ebv" = 'all_u-ebv',
        "mtb" = ,
        "mtbaux" = ,
        "all_u-mtbaux" = ,
        "all_u-mtb" = 'all_u-mtb',
        stop(paste0("stim value of ", stim, " not valid"))
      )
      
      dir_n_clust <- switch(
        n_clust,
        'f' = ,
        'n_f' = 'n_f',
        'f_f' = ,
        'n_f_f' = 'n_f_f',
        stop(paste0("n_clust value of ", n_clust, " not valid"))
      )
      
      # out
      # -----------------
      
      file.path(
        dir_gating,
        dir_gs,
        "root",
        "ebvmtbp1p4uns",
        dir_chnl,
        dir_exc,
        "locb0.15_min_clust",
        "fs",
        dir_non_responders,
        'no_s',
        dir_n_clust, 
        "all_u",
        dir_stim
      )
      
    }