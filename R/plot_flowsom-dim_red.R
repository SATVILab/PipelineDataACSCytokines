#' @title Plot t-SNE and UMAP dim reds with clusters highlighted
#'
#' @param col_clust
#' character.
#' Colour of in-cluster points.
#'
#' @export
plot_dim_red <- function(flowsom_out,
                         path_flowsom_out,
                         dir_save,
                         stim_plot,
                         col_clust,
                         font_size = 14,
                         width = 6,
                         height = 6) {

  unlink(dir_save, recursive = TRUE)
  if (!dir.exists(dir_save)) {
    dir.create(dir_save, recursive = TRUE)
  }

  method_vec <- c("tsne", "umap")
  inc_exc_stim_vec <- c("inc_exc_stim", "exc_exc_stim")
  if (identical(setdiff(c("p1", "ebv", "mtb"), stim_plot), ""[-1])) {
    inc_exc_stim_vec <- inc_exc_stim_vec[1]
  }

  purrr::walk(method_vec, function(method) {
    purrr::walk(inc_exc_stim_vec, function(inc_exc_stim) {
      clust_vec <- flowsom_out %>%
        dplyr::filter(
          stim %in% stim_plot
        ) %>%
        dplyr::pull("clust") %>%
        unique() %>%
        as.numeric() %>%
        sort()

      purrr::walk(clust_vec, function(clust) {
        fn <- paste0(
          "p-", method, "-", inc_exc_stim, "-", clust
        )
        fn <- file.path(
          dir_save, fn
        )
        path_p <- saver::save_rds_eval(
          fn_or_call = .plot_dim_red_clust_rds,
          .data = path_flowsom_out,
          .data_nm = "flowsom_out",
          filename = fn,
          return_obj = FALSE,
          eval = FALSE,
          stim_plot = stim_plot,
          clust = clust,
          method = method,
          col_clust = col_clust,
          inc_exc_stim = inc_exc_stim == "inc_exc_stim"
        )
        p <- saver::load_rds_eval(filename = path_p)
        purrr::walk(c("png", "pdf"), function(gd) {
          cowplot::ggsave2(
            filename = paste0(fn, ".", gd),
            plot = p,
            height = height,
            width = width,
            units = "cm"
          )
        })
      })
    })
  })
}


#' @title Plot dim red view of data
#'
#' @description Plot dimension-reduced view of data,
#' for all stims and UMAP and t-SNE and all clusters
#' individually.
#'
#' @param inc_exc_stim logical.
#' If \code{TRUE}, then non-cluster
#' cells not from the stim condition of concern are plotted.
#' If \code{FALSE}, then they are excluded (and only cells
#' that belong to the stim are plotted).
#' @param col_clust character.
#' Colour of in-cluster points.
#' Default is \code{"red"}.
#' Out-of-cluster points are always \code{"gray75"}.
.plot_dim_red_clust_rds <- function(flowsom_out,
                                    stim_plot,
                                    clust,
                                    method,
                                    col_clust = NULL,
                                    font_size = 10,
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

  dim_vec <- paste0(method, 1:2)

  plot_tbl <- plot_tbl %>%
    dplyr::mutate(
      clust = ifelse(.data$clust == .env$clust,
        paste0("Cluster ", clust), "Other"
      )
    ) %>%
    dplyr::select(
      .data[[dim_vec[1]]],
      .data[[dim_vec[2]]],
      clust,
      stim
    )

  axis_lab <- switch(method,
    "tsne" = paste0("t-SNE ", 1:2),
    "umap" = paste0("UMAP ", 1:2),
  )

  ggplot(
    plot_tbl,
    aes(
      x = .data[[dim_vec[1]]],
      y = .data[[dim_vec[2]]],
      col = clust
    )
    ) +
    cowplot::theme_cowplot(font_size = font_size) +
    geom_point(
      aes(alpha = clust),
      size = 0.25, show.legend = FALSE
      ) +
    scale_colour_manual(values = setNames(
      c(col_clust, "gray75"),
      c(paste0("Cluster ", clust), "Other")
    )) +
    scale_alpha_manual(values = setNames(
      c(0.4, 0.75),
      c(paste0("Cluster ", clust), "Other")
    )) +
    theme(
      legend.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
     ) +
    labs(title = paste0("Cluster ", clust)) +
    coord_equal()
}
