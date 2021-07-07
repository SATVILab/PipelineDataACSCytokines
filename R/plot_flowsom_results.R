#' @title Produce informative plots and results from a FlowSOM clustering for each stimulation combination
#'
#' 
#' @param responders_list list. Output from \code{get_post_probs} function for the current combination 
#' of excluded cytokine combinations, gate method and stimulation (all, all_u or a specific stim). 
#' @param chnl_lab named character vector. If not \code{NULL}, channels are renamed to markers using it. 
#' Default is \code{NULL}. 
#' @param stats_combn_tbl dataframe. Contains total classified cells for each sample (SubjectID and VisitType) for
#' each stimulation in a column called \code{n_cell_stim} and also has the total number of classified cells
#' in the corresponding unstim sample in a column called \code{n_cell_uns}. 
#' @param uns_chr character. Name for unstim in \code{stim_vec_u}. Default is "uns".
#' @param reuse,dir_save_cluster_results logical and character, respectively. If \code{reuse == TRUE} and \code{!is.null(dir_save_cluster_results)}, 
#' then if flowsom_freq is the same as flowsom_freq found 
#' @export 
plot_flowsom <- function(flowsom_out, flowsom_freq, flowsom_cluster_post_probs, stim, dir_save_base, chnl_lab, 
                         stats_combn_tbl, uns_chr = "uns", reuse = FALSE, dir_save_cluster_results = NULL, 
                         chnl_sel){
  
  if(reuse && !is.null(dir_save_cluster_results)){
    path_fs_freq_old <- file.path(dir_save_cluster_results, "fs_freq_old.rds")
    if(file.exists(path_fs_freq_old)){
      flowsom_freq_old <- try(readRDS(path_fs_freq_old))
      if(identical(flowsom_freq_old, flowsom_freq)){
        print("Skipping fs plots as reuse == TRUE and flowsom_freq has not changed.")
        return(invisible(TRUE))
      }
    }
  }
  
  flowsom_out %<>%
    dplyr::mutate(stim = stim %>%
             str_remove("all_u-") %>%
             str_remove("all-"))
  
  # TODO: Incorporate a way to automatically skip plotting if the clustering has not changed. 
  
  # ------------------
  # Preparation
  # ------------------
  
  # get stims to plot
  stim_vec_plot <- unique(c(as.character(flowsom_out$stim), stim)) %>%
    setdiff(uns_chr)
  
  # ------------------
  # create plots for each stimulation
  # ------------------

  purrr::walk(stim_vec_plot, function(stim_curr){
    
    # get directory to save to 
    # -----------
    
    # get how to recognise stim in path
    stim_add <- ifelse(stim_curr %in% c('all', 'all_u'), paste0(stim, "-all"), 
                       paste0(stim, "-", stim_curr))
    
    # get overall directory
    par_list <- attr(flowsom_out, "par")
    dir_save_stim <- file.path(dir_save_base, stim_add, 
                               par_list$scale, 
                               par_list$n_clust)
    if(!dir.exists(dir_save_stim)) dir.create(dir_save_stim, recursive = TRUE)
    
    # choose which stimulations to plot
    stim_plot <- switch(stim_curr, 
                        "all_u" = setdiff(stim_vec_plot, "all_u"),
                        "all" = setdiff(stim_vec_plot, "all"), 
                        stim_curr)
    
    # do dimension reduction plots only 
    # if the current stim is all/all_u or if no stim is all or all_u
    plot_dim_red <- stim_curr %in% c('all', 'all_u') || !any(stim %in% c('all', 'all_u'))
    plot_dim_red <- TRUE
    
    # add unstim if available
    
    # summarise and plot flowsom results for selected stim
    .plot_flowsom(flowsom_out = flowsom_out, 
                  flowsom_freq = flowsom_freq,
                  dir_save = dir_save_stim, 
                  chnl_lab = chnl_lab,
                  chnl_sel = chnl_sel,
                  stim_plot = stim_plot, 
                  dim_red = plot_dim_red)
  }) 
  
  invisible(TRUE)
}



#' @title Plot informative plots to assist cluster interpretation for an individual set of stimulations
#' @param cell_count_lab named numeric vector. Each name indicates sample id and stimulation (paste0(SubjectID, "_", VisitType, 
#' "_", stim)), and each element is the number of classified cells for that combination of sample id and stimulation. 
#' @param stim_plot character vector. Stimulations to plot. If 'uns' and other stimulations are also in flowsom_out$stim, then it will be plotted 
#' where appropriate as well. 
#' @param stim_clust character: 'all', 'all_u' or individual stim. Specifies stimulations that were used to cluster. 
#' @export
.plot_flowsom <- function(flowsom_out, flowsom_freq, dir_save, chnl_lab, chnl_sel, stim_plot, dim_red = TRUE){
  
  # identify specific stimulations needed
  stim_filter_vec <- switch(as.character(any(stim_plot %in% c("all", "all_u"))), 
                           'TRUE' = setdiff(unique(flowsom_out$stim), "uns"), 
                           'FALSE' = stim_plot)
  
  # ===========================
  # Plot frequencies
  # ===========================
  
  purrr::walk(list.files(dir_save, full.names = TRUE, pattern = ".png"), file.remove)
  
  # plot for unstim as well
  plot_fs_cluster_freq(flowsom_freq = flowsom_freq %>% 
                         dplyr::filter(stim %in% c(stim_filter_vec, 'uns')),
                       dir_save = dir_save)
  
  # =============================
  # Plot marker density by cluster
  # =============================

  # plot only for stim to plot
  plot_marker_density_by_cluster(flowsom_out = flowsom_out, chnl_lab = chnl_lab, 
                                 dir_save = dir_save, stim_plot = stim_plot)  
  
  # =============================
  # Plot dimensionality-reduction plots
  # =============================
  
  if(dim_red) plot_dim_red(flowsom_out = flowsom_out, dir_save = dir_save, 
                           stim_plot = stim_plot)
  
  # =============================
  # Plot heat maps
  # =============================
  
  # get empirical cdfs to calculate scaling if not supplied
  ecdf_list <- get_ecdf_list(
    flowsom_out = flowsom_out, 
    chnl_lab = chnl_lab, 
    chnl_sel = chnl_sel
    )
  
  # get order of markers along rows and clusters along columns
  order_list <- get_order_list(
    flowsom_out = flowsom_out, 
    ecdf_list = ecdf_list, 
    chnl_lab = chnl_lab, 
    chnl_sel = chnl_sel
    )
  
  # plot only stable clusters
  plot_flowsom_grid_quantile(
    flowsom_out = flowsom_out, 
    order_list = order_list, 
    ecdf_list = ecdf_list, 
    even_col = TRUE, 
    expand_mid_col = TRUE, 
    stability_min = 0.5, 
    stability_show = TRUE, 
    stim_plot = stim_plot, 
    dir_save = dir_save, 
    chnl_lab = chnl_lab
    )
  
  # plot all clusters
  plot_flowsom_grid_quantile(
    flowsom_out = flowsom_out, 
    order_list = order_list, 
    ecdf_list = ecdf_list, 
    even_col = TRUE, 
    expand_mid_col = TRUE, 
    stability_min = 0, 
    stability_show = TRUE, 
    stim_plot = stim_plot,
    dir_save = dir_save, 
    chnl_lab = chnl_lab
    )
  
  invisible(TRUE)
}

#' @title Plot marker density for each cluster
#' 
#' @param stim_plot character: 'all', 'all_u', or individual stim. If 'all' or 'all_u', then 
#' all cells in flowsom_out belonging to a particular cluster are plotted for that cluster, including unstim. 
#' If an individual stim, then only cells belonging to that stim are plotted for that cluster. 
#' When plotting non-cluster cells, all non-cluster cells are plotted regardless of the value for 
#' \code{stim_clust}. 
#' @export
plot_marker_density_by_cluster <- function(flowsom_out, chnl_lab, dir_save, stim_plot){
  chnl_sel <- attr(flowsom_out, "par_fs")$chnl_sel
  scale_fs <- attr(flowsom_out, "par_fs")$scale
  p_list <- purrr::map(unique(flowsom_out$clust), function(clust){
    flowsom_out_plot_clust_specific <- flowsom_out %>%
      tidyr::pivot_longer(Ba137Di:Time, names_to = "marker", values_to = "expr") %>%
      dplyr::filter(marker %in% chnl_lab[chnl_sel]) %>%
      dplyr::mutate(clust15 = ifelse(.data$clust == .env$clust, 
                              paste0("Cluster ", clust), 'Other'))
    
    
    # remove cells from cluster of interest if they do not belong to a specific stim 
    # of interest (if not plotting all stims used to cluster)
    if(!stim_plot %in% c('all', 'all_u')){
      flowsom_out_plot_clust_specific %<>%
        dplyr::filter(clust15 == 'Other' | stim %in% stim_plot)
    }
    
    if(scale_fs == 'scale'){
      scale_vec <- attr(flowsom_out, "scale")
      med_scale <- median(scale_vec[chnl_sel])
      med_scale_inv <- 1/med_scale
      marker_scale_lab_vec <- setNames(paste0(chnl_lab[chnl_sel], " (", round(1/scale_vec[chnl_sel]/med_scale_inv, 1), ")"), chnl_lab[chnl_sel])
    } else marker_scale_lab_vec <- setNames(chnl_lab[chnl_sel], chnl_lab[chnl_sel])
    
    ggplot(flowsom_out_plot_clust_specific, aes(x = expr)) +
      geom_density(aes(fill = clust15), alpha = 0.5) +
      facet_wrap(~marker, ncol = 4, scales = 'free', 
                 labeller = labeller(marker = marker_scale_lab_vec)) +
      labs(title = paste0("Cluster ", clust)) +
      scale_fill_brewer(palette = "Set1") +
      theme(legend.title = element_blank()) +
      labs(x = "Marker expression", y = "Density")
  })%>%
    setNames(paste0("p markers c", unique(flowsom_out$clust)))
  
  analysispipeline::save_objects(obj_list = p_list, 
                               dir_proj = dir_save, 
                               width = 20, 
                               height = ceiling(length(chnl_sel)/4) * 5, 
                               empty = FALSE)
  
  rm('p_list_scaled')
  rm('p_list_unscaled')
  
  invisible(TRUE)
}

#' @title Plot cluster frequencies per cyt+ and all cells, by stim and prog
#' 
#' @export
plot_fs_cluster_freq <- function(flowsom_freq, dir_save){
  
  # ------------------------
  # Plot freq_cyt - clusters together
  # ------------------------
  
  mult_stim <- length(setdiff(flowsom_freq$stim, "uns")) > 1
  
  # vector to label stims
  stim_lab_vec <- c("uns" = "Unstim", 
                    "ebv" = "EBV/CMV", 
                    "mtb" = "Mtb Aux", 
                    "p1" = "P1", 
                    "p4" = "P4")
  
  # across all stims
  plot_tbl <- flowsom_freq 
  plot_tbl_no_uns <- plot_tbl %>% dplyr::filter(stim != 'uns')

  stability_tbl <- flowsom_freq %>%
    dplyr::group_by(clust) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(clust, stability) %>%
    dplyr::mutate(stability = round(stability, 2))
  
  stability_lab_vec <- setNames(paste0(stability_tbl$clust, " (", stability_tbl$stability,
                                       ")"), stability_tbl$clust)
 # factor(stability_lab_vec, levels = as.character(sort(as.numeric(names(stability_lab_vec)))))
  
  p_freq_cyt <- ggplot(plot_tbl_no_uns, 
                       aes(x = Progressor, y = freq_cyt_stim)) + 
    geom_hline(yintercept = 0) +
    cowplot::background_grid(major = 'y') +
    geom_boxplot(aes(fill = Progressor), outlier.size = -1) +
    geom_point(size = 0.5) +
    
    scale_fill_manual(values = c("Progressor" = "orange",
                                 "Control" = "dodgerblue")) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          legend.title = element_blank()) +
    labs(x = "Progressor status", y = "Frequency of cytokine-positive cells")
  
  p_freq_cyt_free <- p_freq_cyt + 
    facet_wrap(~clust, scales = 'free')
  p_freq_cyt_free_asn <- p_freq_cyt_free +
    scale_y_continuous(trans = 'asn')
  p_freq_cyt_fixed <- p_freq_cyt + 
    facet_wrap(~clust, scales = 'fixed')
  
  if(mult_stim){
    
    stim_lab_vec <- c("uns" = "Unstim", 
                      "ebv" = "EBV/CMV", 
                      "mtb" = "Mtb Aux", 
                      "p1" = "P1", 
                      "p4" = "P4")
    p_freq_cyt_stim <- ggplot(plot_tbl_no_uns, 
                              aes(x = Progressor, y = freq_cyt_stim)) + 
      geom_hline(yintercept = 0) +
      cowplot::background_grid(major = 'y') +
      geom_boxplot(aes(fill = stim), outlier.size = -1, position =) +
      #geom_point(size = 0.5) +
      scale_fill_manual(values = c("uns" = "#a6cee3", 
                                   "ebv" = "#1f78b4", 
                                   "mtbaux" = "#b2df8a", 
                                   "mtb" = "#b2df8a", 
                                   "p1" = "#33a02c"), 
                        labels = stim_lab_vec) +
      theme(legend.title = element_blank(), 
            legend.position = 'bottom') +
      labs(x = "Progressor status", y = "Frequency of cytokine-positive cells")
    
    p_freq_cyt_stim_free <- p_freq_cyt_stim + 
      facet_wrap(~clust, scales = 'free_y', ncol = 3)
    p_freq_cyt_stim_free_asn <- p_freq_cyt_stim_free + 
      scale_y_continuous(trans = 'asn')
    p_freq_cyt_stim_fixed <- p_freq_cyt_stim + 
      facet_wrap(~clust, scales = 'fixed')
  }
  
  # ------------------------
  # Plot freq_tot - clusters separate
  # ------------------------
  
  # across all stims
  p_freq_tot <- ggplot(plot_tbl_no_uns, 
                       aes(x = Progressor, y = freq_tot_stim)) + 
    geom_hline(yintercept = 0) +
    cowplot::background_grid(major = 'y') +
    geom_boxplot(aes(fill = Progressor), outlier.size = -1) +
    geom_point(size = 0.5) +
    
    scale_fill_manual(values = c("Progressor" = "orange",
                                 "Control" = "dodgerblue")) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          legend.title = element_blank(), 
          legend.position = 'bottom') +
    labs(x = "Progressor status", y = "Frequency of all cells")
  
  p_freq_tot_free <- p_freq_tot + 
    facet_wrap(~clust, scales = 'free')
  p_freq_tot_free_asn <- p_freq_tot_free +
    scale_y_continuous(trans = 'asn')
  p_freq_tot_fixed <- p_freq_tot + 
    cowplot::background_grid(minor = 'y') +
    facet_wrap(~clust, scales = 'fixed') +
    scale_y_continuous(trans = 'asn')
  
  if(mult_stim){
    # across all stims
    p_freq_tot_stim <- ggplot(plot_tbl, 
                              aes(x = Progressor, y = freq_tot_stim)) + 
      geom_hline(yintercept = 0) +
      cowplot::background_grid(major = 'y') +
      geom_boxplot(aes(fill = stim), outlier.size = 0.25, outlier.colour = 'gray65') +
      scale_fill_manual(values = c("Progressor" = "orange",
                                   "Control" = "dodgerblue")) +
      theme(legend.title = element_blank(), 
            legend.position = 'bottom') +
      scale_fill_manual(values = c("uns" = "#a6cee3", 
                                   "ebv" = "#1f78b4", 
                                   "mtbaux" = "#b2df8a", 
                                   "mtb" = "#b2df8a", 
                                   "p1" = "#33a02c"), 
                        labels = stim_lab_vec) +
      labs(x = "Progressor status", y = "Frequency of all cells")
    
    p_freq_tot_stim_free <- p_freq_tot_stim + 
      facet_wrap(~clust, scales = 'free_y', ncol = 3) 
    p_freq_tot_stim_free_asn <- p_freq_tot_stim_free + 
      scale_y_continuous(trans = 'asn')
    p_freq_tot_stim_fixed <- p_freq_tot_stim + 
      cowplot::background_grid(minor = 'y') +
      facet_wrap(~clust, scales = 'fixed') +
      scale_y_continuous(trans = 'asn')
    
    
    analysispipeline::save_objects(p_freq_tot_stim_fixed = p_freq_tot_stim_fixed, 
                                 p_freq_cyt_stim_fixed = p_freq_cyt_stim_fixed,
                                 dir_proj = dir_save, 
                                 width = 20, 
                                 height = length(stability_lab_vec), 
                                 empty = FALSE)
    
    analysispipeline::save_objects(p_freq_cyt_stim_free = p_freq_cyt_stim_free, 
                                 p_freq_cyt_stim_free_asn = p_freq_cyt_stim_free_asn,
                                 p_freq_tot_stim_free = p_freq_tot_stim_free, 
                                 p_freq_tot_stim_free_asn = p_freq_tot_stim_free_asn,
                                 dir_proj = dir_save, 
                                 width = 20, 
                                 height = length(stability_lab_vec) * 1.25, 
                                 empty = FALSE)
  }
  
  analysispipeline::save_objects(p_freq_tot_fixed = p_freq_tot_fixed, 
                               p_freq_cyt_fixed = p_freq_cyt_fixed,
                               dir_proj = dir_save, 
                               width = 20, 
                               height = length(stability_lab_vec), 
                               empty = FALSE)
  
  analysispipeline::save_objects(p_freq_cyt_free = p_freq_cyt_free,
                               p_freq_cyt_free_asn = p_freq_cyt_free_asn,
                               p_freq_tot_free = p_freq_tot_free, 
                               p_freq_tot_free_asn = p_freq_tot_free_asn,
                               dir_proj = dir_save, 
                               width = 20, 
                               height = length(stability_lab_vec) * 1.25, 
                               empty = FALSE)
  
  # ------------------------
  # Plot freq_tot_bs - clusters together
  # ------------------------
  
  if('uns' %in% plot_tbl$stim){
    p_freq_bs_tot <- ggplot(plot_tbl_no_uns, 
                            aes(x = Progressor, y = freq_tot_bs)) + 
      geom_hline(yintercept = 0) +
      cowplot::background_grid(major = 'y') +
      geom_boxplot(aes(fill = Progressor), outlier.size = -1) +
      geom_point(size = 0.5) +
      
      scale_fill_manual(values = c("Progressor" = "orange",
                                   "Control" = "dodgerblue")) +
      theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.title = element_blank(), 
            legend.position = 'bottom') +
      labs(x = "Progressor status", y = "Background-subtracted frequency of all cells") +
      expand_limits(y = 0)
    
    p_freq_bs_tot_free <- p_freq_bs_tot + 
      facet_wrap(~clust, scales = 'free')
    p_freq_bs_tot_fixed <- p_freq_bs_tot + 
      cowplot::background_grid(minor = 'y') +
      facet_wrap(~clust, scales = 'fixed') +
      scale_y_continuous(trans = 'asn')
    
    if(mult_stim){
      # across all stims
      p_freq_bs_tot_stim <- ggplot(plot_tbl_no_uns, 
                                   aes(x = Progressor, y = freq_tot_bs)) + 
        geom_hline(yintercept = 0, size = 1) +
        cowplot::background_grid(major = 'y') +
        geom_hline(yintercept = 0, linetype = 'solid') +
        geom_boxplot(aes(fill = stim), outlier.size = 0.25, outlier.colour = 'gray65') +
        scale_fill_manual(values = c("Progressor" = "orange",
                                     "Control" = "dodgerblue")) +
        theme(legend.title = element_blank(), 
              legend.position = 'bottom') +
        scale_fill_manual(values = c("uns" = "#a6cee3", 
                                     "ebv" = "#1f78b4", 
                                     "mtbaux" = "#b2df8a", 
                                     "mtb" = "#b2df8a", 
                                     "p1" = "#33a02c"), 
                          labels = stim_lab_vec) +
        labs(x = "Progressor status", y = "Background-subtracted frequency of all cells") +
        expand_limits(y = 0)
      
      p_freq_bs_tot_stim_free <- p_freq_bs_tot_stim + 
        facet_wrap(~clust, scales = 'free_y', ncol = 3) 
      p_freq_bs_tot_stim_free_asn <- p_freq_bs_tot_stim_free + 
        scale_y_continuous(trans = 'asn')
      p_freq_bs_tot_stim_fixed <- p_freq_bs_tot_stim + 
        cowplot::background_grid(minor = 'y') +
        facet_wrap(~clust, scales = 'fixed') +
        scale_y_continuous(trans = 'asn')
      
      analysispipeline::save_objects(p_freq_bs_tot_stim_fixed = p_freq_bs_tot_stim_fixed,
                                   dir_proj = dir_save, 
                                   width = 20, 
                                   height = length(stability_lab_vec), 
                                   empty = FALSE)
      
      analysispipeline::save_objects(p_freq_bs_tot_stim_free = p_freq_bs_tot_stim_free,
                                   p_freq_bs_tot_stim_free_asn = p_freq_bs_tot_stim_free_asn, 
                                   dir_proj = dir_save, 
                                   width = 20, 
                                   height = length(stability_lab_vec) * 1.25, 
                                   empty = FALSE)
    }
    
    analysispipeline::save_objects(p_freq_bs_tot_free = p_freq_bs_tot_free, 
                                 p_freq_bs_tot_fixed = p_freq_bs_tot_fixed,
                                 dir_proj = dir_save, 
                                 width = 20, 
                                 height = length(stability_lab_vec), 
                                 empty = FALSE)
  }
  
  
  # ------------------------
  # Plot freq_cyt - clusters separate
  # ------------------------
  
  p_freq_cyt_list_free <- purrr::map(unique(plot_tbl_no_uns$clust), function(clust){
    ggplot(plot_tbl_no_uns %>%
             dplyr::filter(clust == .env$clust), 
           aes(x = Progressor, y = freq_cyt_stim)) + 
      geom_hline(yintercept = 0) +
      cowplot::background_grid(major = 'y') +
      geom_boxplot(aes(fill = Progressor), outlier.size = -1) +
      geom_point(size = 0.5) +
      
      scale_fill_manual(values = c("Progressor" = "orange",
                                   "Control" = "dodgerblue")) +
      theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.title = element_blank()) +
      labs(x = "Progressor status", y = "Frequency of\ncytokine-positive cells")
  }) %>%
    setNames(paste0("p cyt prog c ", unique(plot_tbl_no_uns$clust), " free"))
  
  p_freq_cyt_list_fixed <- purrr::map(p_freq_cyt_list_free, function(p){
    range <- diff(range(plot_tbl_no_uns$freq_cyt_stim))
    p + lims(y = c(0, max(plot_tbl_no_uns$freq_cyt_stim) + 0.025 * range))
  }) %>%
    setNames(paste0("p cyt prog c ", unique(plot_tbl_no_uns$clust), " fixed"))
  
  analysispipeline::save_objects(obj_list = p_freq_cyt_list_free %>%
                                 append(p_freq_cyt_list_fixed), 
                               dir_proj = dir_save, 
                               width = 10, 
                               height = 6, 
                               empty = FALSE)
  
  
  if(mult_stim){
    p_freq_cyt_stim_list_free <- purrr::map(unique(plot_tbl$clust), function(clust){
      ggplot(plot_tbl_no_uns %>% dplyr::filter(clust == .env$clust), 
             aes(x = Progressor, y = freq_cyt_stim)) + 
        geom_hline(yintercept = 0) +
        cowplot::background_grid(major = 'y') +
        geom_boxplot(aes(fill = stim), outlier.size = 0.25, outlier.colour = 'gray66') +
        scale_fill_manual(values = c("uns" = "#a6cee3", 
                                     "ebv" = "#1f78b4", 
                                     "mtbaux" = "#b2df8a", 
                                     "mtb" = "#b2df8a", 
                                     "p1" = "#33a02c"), 
                          labels = stim_lab_vec) +
        theme(legend.title = element_blank()) +
        labs(x = "Progressor status", y = "Frequency of\ncytokine-positive cells")
    }) %>%
      setNames(paste0("p cyt prog stim c ", unique(plot_tbl$clust), " free"))
    p_freq_cyt_stim_list_fixed <- purrr::map(p_freq_cyt_stim_list_free, function(p){
      range <- diff(range(plot_tbl$freq_cyt_stim))
      p + lims(y = c(0, max(plot_tbl$freq_cyt_stim) + 0.025 * range))
    }) %>%
      setNames(paste0("p cyt prog stim c ", unique(plot_tbl$clust), " fixed"))
    analysispipeline::save_objects(obj_list = p_freq_cyt_stim_list_free %>%
                                   append(p_freq_cyt_stim_list_fixed), 
                                 dir_proj = dir_save, 
                                 width = 10, 
                                 height = 6, 
                                 empty = FALSE)
  }
  
  # ------------------------
  # Plot freq_tot - clusters separate
  # ------------------------
  
  p_freq_tot_list_free <- purrr::map(unique(plot_tbl$clust), function(clust){
    ggplot(plot_tbl %>%
             dplyr::filter(clust == .env$clust), 
           aes(x = Progressor, y = freq_tot_stim)) + 
      geom_hline(yintercept = 0) +
      cowplot::background_grid(major = 'y') +
      geom_boxplot(aes(fill = Progressor), outlier.size = -1) +
      geom_point(size = 0.5) +
      
      scale_fill_manual(values = c("Progressor" = "orange",
                                   "Control" = "dodgerblue")) +
      theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.title = element_blank()) +
      labs(x = "Progressor status", y = "Frequency of\ntotal cells")
  }) %>%
    setNames(paste0("p tot prog c ", unique(plot_tbl$clust), " free"))
  
  p_freq_tot_list_fixed <- purrr::map(p_freq_tot_list_free, function(p){
    range <- diff(range(plot_tbl$freq_tot_stim))
    p + lims(y = c(0, max(plot_tbl$freq_tot_stim) + 0.025 * range))
  }) %>%
    setNames(paste0("p tot prog c ", unique(plot_tbl$clust), " fixed"))
  
  analysispipeline::save_objects(obj_list = p_freq_tot_list_free %>%
                                 append(p_freq_tot_list_fixed), 
                               dir_proj = dir_save, 
                               width = 10, 
                               height = 6, 
                               empty = FALSE)
  
  
  if(mult_stim){
    p_freq_tot_stim_list_free <- purrr::map(unique(plot_tbl$clust), function(clust){
      ggplot(plot_tbl %>% dplyr::filter(clust == .env$clust), 
             aes(x = Progressor, y = freq_tot_stim)) + 
        geom_hline(yintercept = 0) +
        cowplot::background_grid(major = 'y') +
        geom_boxplot(aes(fill = stim), outlier.size = 0.25, outlier.colour = 'gray66') +
        scale_fill_manual(values = c("uns" = "#a6cee3", 
                                     "ebv" = "#1f78b4", 
                                     "mtbaux" = "#b2df8a", 
                                     "mtb" = "#b2df8a", 
                                     "p1" = "#33a02c"), 
                          labels = stim_lab_vec) +
        theme(legend.title = element_blank()) +
        labs(x = "Progressor status", y = "Frequency of\ntotal cells")
    }) %>%
      setNames(paste0("p tot prog stim c ", unique(plot_tbl$clust), " free"))
    
    p_freq_tot_stim_list_fixed <- purrr::map(p_freq_tot_stim_list_free, function(p){
      range <- diff(range(plot_tbl$freq_tot_stim))
      p + lims(y = c(0, max(plot_tbl$freq_tot_stim) + 0.025 * range))
    }) %>%
      setNames(paste0("p tot prog stim c ", unique(plot_tbl$clust), " fixed"))
    analysispipeline::save_objects(obj_list = p_freq_tot_stim_list_free %>%
                                   append(p_freq_tot_stim_list_fixed), 
                                 dir_proj = dir_save, 
                                 width = 10, 
                                 height = 6, 
                                 empty = FALSE)
  }
  
  # ------------------------
  # Plot freq_tot_bs - clusters separate
  # ------------------------
  
  if('uns' %in% unique(plot_tbl$stim)){
    p_freq_tot_bs_list_free <- purrr::map(unique(plot_tbl_no_uns$clust), function(clust){
      ggplot(plot_tbl_no_uns %>%
               dplyr::filter(clust == .env$clust), 
             aes(x = Progressor, y = freq_tot_bs)) + 
        cowplot::background_grid(major = 'y') +
        geom_boxplot(aes(fill = Progressor), outlier.size = -1) +
        geom_point(size = 0.5) +
        
        scale_fill_manual(values = c("Progressor" = "orange",
                                     "Control" = "dodgerblue")) +
        theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              legend.title = element_blank()) +
        labs(x = "Progressor status", y = "Background-subtracted\nfrequency of total cells")
    }) %>%
      setNames(paste0("p bs prog c ", unique(plot_tbl_no_uns$clust), " free"))
    
    p_freq_tot_bs_list_fixed <- purrr::map(p_freq_tot_list_free, function(p){
      range <- diff(range(plot_tbl_no_uns$freq_tot_bs))
      p + lims(y = c(0, max(plot_tbl_no_uns$freq_tot_bs) + 0.025 * range))
    }) %>%
      setNames(paste0("p bs prog c ", unique(plot_tbl_no_uns$clust), " fixed"))
    
    analysispipeline::save_objects(obj_list = p_freq_tot_bs_list_free %>%
                                   append(p_freq_tot_bs_list_fixed), 
                                 dir_proj = dir_save, 
                                 width = 10, 
                                 height = 6, 
                                 empty = FALSE)
    
    
    if(mult_stim){
      p_freq_tot_bs_stim_list_free <- purrr::map(unique(plot_tbl_no_uns$clust), function(clust){
        ggplot(plot_tbl_no_uns %>% dplyr::filter(clust == .env$clust), 
               aes(x = Progressor, y = freq_tot_bs)) + 
          cowplot::background_grid(major = 'y') +
          geom_boxplot(aes(fill = stim), outlier.size = 0.25, outlier.colour = 'gray66') +
          scale_fill_manual(values = c("uns" = "#a6cee3", 
                                       "ebv" = "#1f78b4", 
                                       "mtbaux" = "#b2df8a", 
                                       "mtb" = "#b2df8a", 
                                       "p1" = "#33a02c"), 
                            labels = stim_lab_vec) +
          theme(legend.title = element_blank()) +
          labs(x = "Progressor status", y = "Background-subtracted\nfrequency of total cells")
      }) %>%
        setNames(paste0("p bs prog stim c ", unique(plot_tbl_no_uns$clust), " free"))
      
      p_freq_tot_bs_stim_list_fixed <- purrr::map(p_freq_tot_bs_stim_list_free, function(p){
        range <- diff(range(plot_tbl_no_uns$freq_tot_bs))
        p + lims(y = c(0, max(plot_tbl_no_uns$freq_tot_bs) + 0.025 * range))
      }) %>%
        setNames(paste0("p bs prog stim c ", unique(plot_tbl_no_uns$clust), " fixed"))
      analysispipeline::save_objects(obj_list = p_freq_tot_bs_stim_list_free %>%
                                     append(p_freq_tot_bs_stim_list_fixed), 
                                   dir_proj = dir_save, 
                                   width = 12, 
                                   height = 7.25, 
                                   empty = FALSE)
    }
  }
  
  invisible(TRUE)
  
}

#' @title Plot t-SNE and UMAP dim reds with clusters highlighted
#' @export
plot_dim_red <- function(flowsom_out, dir_save, stim_plot){
  
  flowsom_out_stim <- flowsom_out %>% dplyr::filter(stim %in% stim_plot)
  
  clust_vec <- sort(as.numeric(unique(flowsom_out$clust)))
  
  # t-SNE
  tsne_clust_list <- purrr::map(clust_vec, function(clust){
    
    plot_tbl <- flowsom_out %>%
      dplyr::mutate(clust = ifelse(.data$clust == .env$clust, 
                            paste0("Cluster ", clust), "Other")) %>%
      dplyr::select(tsne1, tsne2, clust, stim)
    
    # remove cells from cluster of interest if they do not belong to a specific stim 
    # of interest (if not plotting all stims used to cluster)
    if(!stim_plot %in% c('all', 'all_u')){
      plot_tbl %<>%
        dplyr::filter(clust == 'Other' | stim %in% stim_plot)
    }
    
    ggplot(plot_tbl, 
           aes(x = tsne1, y = tsne2, col = clust)) + 
      cowplot::theme_cowplot(font_size = 10) +
      geom_point(alpha = 0.5, size = 0.25, show.legend = FALSE)+
      scale_colour_manual(values = setNames(c("red", "gray75"), 
                                            c(paste0("Cluster ", clust), "Other"))) +
      theme(legend.title = element_blank()) +
      labs(x = "t-SNE 1", y = "t-SNE 2", title = paste0("Cluster ", clust)) +
      coord_equal() 
  }) %>%
    setNames(paste0("p tsne c", clust_vec))
  
  tsne_clust_grid <- cowplot::plot_grid(plotlist = purrr::map(tsne_clust_list, function(x){
    x +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), 
            axis.title = element_blank())
  }), 
  ncol = 5)
  
  analysispipeline::save_objects(obj_list = tsne_clust_list, 
                               dir_proj = dir_save, 
                               width = 6 * 3, 
                               height = 6 * 3, 
                               empty = FALSE)
  
  
  analysispipeline::save_objects(tsne_clust_grid = tsne_clust_grid, 
                               dir_proj = dir_save, 
                               width = 19.85 * 2, 
                               height = ceiling(length(clust_vec) / 5) * 4.5 * 2, 
                               empty = FALSE)
  
  rm('tsne_clust_list')
  rm('tsne_clust_grid')
  
  # umap
  umap_clust_list <- purrr::map(clust_vec, function(clust){
    
    plot_tbl <- flowsom_out %>%
      dplyr::mutate(clust = ifelse(.data$clust == .env$clust, 
                            paste0("Cluster ", clust), "Other")) %>%
      dplyr::select(umap1, umap2, clust, stim)
    
    if(!stim_plot %in% c('all', 'all_u')){
      plot_tbl %<>%
        dplyr::filter(clust == 'Other' | stim %in% stim_plot)
    }
    
    ggplot(plot_tbl, 
           aes(x = umap1, y = umap2, col = clust)) + 
      cowplot::theme_cowplot(font_size = 10) +
      geom_point(alpha = 0.5, size = 0.25, show.legend = FALSE)+
      scale_colour_manual(values = setNames(c("red", "gray75"), 
                                            c(paste0("Cluster ", clust), "Other"))) +
      theme(legend.title = element_blank()) +
      labs(x = "UMAP 1", y = "UMAP 2", title = paste0("Cluster ", clust)) +
      coord_equal() 
  }) %>%
    setNames(paste0("p umap c",clust_vec))
  
  umap_clust_grid <- cowplot::plot_grid(plotlist = purrr::map(umap_clust_list, function(x){
    x +
      theme(axis.text = element_blank(), axis.ticks = element_blank(), 
            axis.title = element_blank())
  }), 
  ncol = 5)
  
  analysispipeline::save_objects(obj_list = umap_clust_list, 
                               dir_proj = dir_save, 
                               width = 6 * 3, 
                               height = 6 * 3, 
                               empty = FALSE)


  
  analysispipeline::save_objects(umap_clust_grid = umap_clust_grid, 
                               dir_proj = dir_save, 
                               width = 19.85 * 2, 
                               height = ceiling(length(clust_vec) / 5) * 4.5 * 2, 
                               empty = FALSE)
  
  rm('umap_clust_list')
  rm('umap_clust_grid')
  
  invisible(TRUE)
}

# TODO: this function
plot_fs_post_probs <- function(fs_post_probs){
  
}