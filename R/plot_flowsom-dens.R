#' @title Plot marker density for each cluster
#'
#' @param stim_plot character: "all", "all_u", or individual stim.
#' If "all" or "all_u", then
#' all cells in flowsom_out belonging to a particular cluster
#' are plotted for that cluster, including unstim.
#' If an individual stim, then only cells belonging
#' to that stim are plotted for that cluster.
#' When plotting non-cluster cells, all
#' non-cluster cells are plotted regardless of the value for
#' \code{stim_clust}.
#' @export
plot_marker_density_by_cluster <- function(flowsom_out,
                                           path_flowsom_out,
                                           chnl_lab,
                                           dir_save,
                                           stim_plot,
                                           font_size) {
  unlink(dir_save, recursive = TRUE)
  if (!dir.exists(dir_save)) {
    dir.create(dir_save, recursive = TRUE)
  }
  chnl_sel <- attr(flowsom_out, "par_fs")$chnl_sel
  scale_fs <- attr(flowsom_out, "par_fs")$scale
  clust_vec <- unique(flowsom_out$clust)
  inc_exc_stim_vec <- c("inc_exc_stim", "exc_exc_stim")
  if (identical(setdiff(c("p1", "ebv", "mtb"), stim_plot), ""[-1])) {
    inc_exc_stim_vec <- inc_exc_stim_vec[1]
  }

  for (inc_exc_stim in inc_exc_stim_vec) {
    for (clust in clust_vec) {
      fn <- paste0("p-", inc_exc_stim, "-", clust)
      fn <- file.path(dir_save, fn)
      path_p <- saver::save_rds_eval(
        fn_or_call = .plot_mk_dens_clust_rds,
        .data = path_flowsom_out,
        .data_nm = "flowsom_out",
        filename = fn,
        return_obj = FALSE,
        eval = FALSE,
        chnl_sel = chnl_sel,
        chnl_lab = chnl_lab,
        stim_plot = stim_plot,
        font_size = font_size,
        inc_exc_stim = inc_exc_stim == "inc_exc_stim",
        clust = clust
      )
      p <- saver::load_rds_eval(filename = path_p)
      for (gd in c("png", "pdf")) {
        cowplot::ggsave2(
          filename = paste0(fn, ".", gd),
          plot = p,
          height = ceiling(length(chnl_sel) / 4) * 5,
          width = 20,
          units = "cm"
        )
      }
    }
  }

  invisible(TRUE)
}

.plot_mk_dens_clust_rds <- function(flowsom_out, # provide path here
                                    clust,
                                    chnl_sel,
                                    chnl_lab,
                                    stim_plot,
                                    font_size,
                                    inc_exc_stim) {
  plot_tbl <- flowsom_out %>%
    dplyr::filter(
      (stim %in% stim_plot &
        .data$clust == .env$clust) |
        .data$clust != .env$clust
    )

  if (!inc_exc_stim) {
    plot_tbl <- plot_tbl %>%
      dplyr::filter(
        (stim %in% stim_plot &
          .data$clust != .env$clust) |
          .data$clust == .env$clust
      )
  }

  plot_tbl <- plot_tbl %>%
    dplyr::mutate(
      clust = ifelse(.data$clust == .env$clust,
        paste0("Cluster ", clust), "Other"
      )
    )

  marker_sel_vec <- sort(chnl_lab[chnl_sel])
  plot_tbl <- plot_tbl[, colnames(plot_tbl) %in%
    c("clust", marker_sel_vec)]

  plot_tbl <- plot_tbl %>%
    tidyr::pivot_longer(
      -clust,
      names_to = "marker",
      values_to = "expr"
    )

  ggplot(plot_tbl, aes(x = expr)) +
    cowplot::theme_cowplot(font_size = font_size) +
    geom_density(aes(fill = clust), alpha = 0.5) +
    facet_wrap(~marker,
      ncol = 4, scales = "free"
    ) +
    labs(title = paste0("Cluster ", clust)) +
    scale_fill_brewer(palette = "Set1") +
    theme(
      legend.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    ) +
    labs(x = "Marker expression", y = "Density")
}