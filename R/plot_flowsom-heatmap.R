#' @title Get list of vectors ordering markers and clusters along heat map
get_order_list <- function(flowsom_out,
                           ecdf_list = NULL,
                           chnl_lab,
                           chnl_sel = NULL) {
  flowsom_out_long <- flowsom_out %>%
    tidyr::pivot_longer(
      cols = -c(clust:tsne2),
      names_to = "marker",
      values_to = "expr"
    )

  if (is.null(ecdf_list)) {
    ecdf_list <- get_ecdf_list(
      flowsom_out = flowsom_out,
      chnl_lab = chnl_lab,
      chnl_sel = chnl_sel
    )
  }

  if (is.null(chnl_sel)) {
    marker_sel_vec <- purrr::map_chr(
      ecdf_list,
      function(x) unique(x$marker)
    ) %>%
      unlist() %>%
      unique()
  } else {
    marker_sel_vec <- chnl_lab[chnl_sel]
  }

  flowsom_out_long_grid <- purrr::map_df(
    unique(flowsom_out_long$clust), function(clust) {
      flowsom_out_long_clust <- flowsom_out_long %>%
        dplyr::filter(.data$clust == .env$clust)
      flowsom_out_long_clust_non <- flowsom_out_long %>%
        dplyr::filter(.data$clust != .env$clust)
      purrr::map_df(marker_sel_vec, function(marker) {
        flowsom_out_long_clust_marker <- flowsom_out_long_clust %>%
          dplyr::filter(.data$marker == .env$marker)
        med <- median(flowsom_out_long_clust_marker$expr)
        ecdf_other <- ecdf_list[[clust]][[marker]]
        flowsom_out_long_clust_marker[1, ] %>%
          dplyr::mutate(perc = ecdf_other(med) * 1e2) %>%
          dplyr::select(clust, marker, perc)
      })
    }
  )


  expr_tbl <- flowsom_out_long_grid %>%
    tidyr::pivot_wider(names_from = marker, values_from = perc)

  clust_vec <- expr_tbl$clust

  expr_mat_clust <- expr_tbl %>%
    dplyr::select(-clust) %>%
    as.matrix()

  rownames(expr_mat_clust) <- clust_vec

  expr_mat_marker <- t(expr_mat_clust)


  hclust_obj_clust <- hclust(dist(expr_mat_clust))
  hclust_obj_marker <- hclust(dist(expr_mat_marker))

  order_vec_list <- list(
    "marker" = rownames(expr_mat_marker)[hclust_obj_marker$order],
    "cluster" = rownames(expr_mat_clust)[hclust_obj_clust$order]
  )

  order_vec_list
}

#' @title Get ecdf for marker expression on non-cluster cells for each cluster
get_ecdf_list <- function(flowsom_out, chnl_lab, chnl_sel) {
  flowsom_out_long <- flowsom_out %>%
    tidyr::pivot_longer(
      cols = -c(clust:tsne2),
      names_to = "marker",
      values_to = "expr"
    )

  ecdf_list <- purrr::map(unique(flowsom_out_long$clust), function(clust) {
    flowsom_out_long_clust_non <- flowsom_out_long %>%
      dplyr::filter(.data$clust != .env$clust)
    purrr::map(chnl_lab[chnl_sel], function(marker) {
      ecdf(flowsom_out_long_clust_non %>%
        dplyr::filter(.data$marker == .env$marker) %>%
        dplyr::pull(expr))
    }) %>%
      setNames(chnl_lab[chnl_sel])
  }) %>%
    setNames(unique(flowsom_out_long$clust))
}

#' @title Plot heat map of marker expression per cluster
#'
#' @param stability_show logical.
#' If TRUE, then stability of each cluster is shown
#' beneath its text on the x-axis.
#' @param stability_min [0,1].
#' Only clusters with stability greater than this value are shown.
plot_clust_map <- function(flowsom_out,
                           path_flowsom_out,
                           dir_save,
                           stim_plot,
                           even_col = TRUE,
                           expand_mid_col = FALSE,
                           stability_show = FALSE,
                           stability_min = 0,
                           col_scheme = "RdYlBu",
                           chnl_lab,
                           chnl_sel = NULL) {
  unlink(dir_save, recursive = TRUE)
  if (!dir.exists(dir_save)) {
    dir.create(dir_save, recursive = TRUE)
  }
  chnl_sel <- attr(flowsom_out, "par_fs")$chnl_sel

  inc_exc_stim_vec <- c("inc_exc_stim", "exc_exc_stim")
  if (identical(
    setdiff(c("p1", "ebv", "mtb"), stim_plot),
    ""[-1]
  )) {
    inc_exc_stim_vec <- inc_exc_stim_vec[1]
  }

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

  for (inc_exc_stim in inc_exc_stim_vec) {
    for (stab_min in stability_min) {
      fn <- paste0(
        "p-", inc_exc_stim, "-",
        ifelse(stab_min == 0, "stb_n", "stb")
      )
      fn <- file.path(
        dir_save, fn
      )
      path_p <- saver::save_rds_eval(
        fn_or_call = .plot_fs_clust_map_rds,
        .data = path_flowsom_out,
        .data_nm = "flowsom_out",
        filename = fn,
        return_obj = FALSE,
        eval = FALSE,
        stim_plot = stim_plot,
        order_list = order_list,
        ecdf_list = ecdf_list,
        even_col = even_col,
        expand_mid_col = expand_mid_col,
        stability_show = stability_show,
        stab_min = stab_min,
        col_scheme = col_scheme,
        chnl_lab = chnl_lab,
        chnl_sel = chnl_sel
      )
      p <- saver::load_rds_eval(filename = path_p)

      n_clust_plot <- length(unique(p$data$clust))
      n_marker_plot <- length(unique(p$data$marker))
      for (gd in c("png", "pdf")) {
        cowplot::ggsave2(
          filename = paste0(fn, ".", gd),
          plot = p,
          height = n_marker_plot * 0.4 + 3,
          width = min(15, n_clust_plot) * 1.5 + 4,
          units = "cm"
        )
      }
    }
  }
}

.plot_fs_clust_map_rds <- function(flowsom_out,
                                   order_list = NULL,
                                   ecdf_list = NULL,
                                   even_col = TRUE,
                                   expand_mid_col = FALSE,
                                   stability_show = TRUE,
                                   stab_min = 0,
                                   col_scheme = "RdYlBu",
                                   chnl_lab,
                                   chnl_sel = NULL,
                                   stim_plot,
                                   dir_save) {
  if (stab_min > 0) {
    if (stab_min > 1) {
      warning("stability_mean greater than 1 in plot_flowsom_grid_quantile,
      and so not applied.")
      stab_min <- 0
    } else {
      if ("stability" %in% colnames(flowsom_out)) {
        flowsom_out <- flowsom_out %>%
          dplyr::filter(.data$stability > stab_min)
      }
      if (nrow(flowsom_out) == 0) {
        warning(paste0(
          "No clusters had stability greater than ",
          stab_min, ", and so a blank plot is returned."
        ))
        return(ggplot2::ggplot())
      }
    }
  }

  flowsom_out_long <- flowsom_out %>%
    tidyr::pivot_longer(
      cols = -c(clust:tsne2),
      names_to = "marker",
      values_to = "expr"
    )

  # get empirical cdfs to calculate scaling if not supplied
  if (is.null(ecdf_list)) {
    ecdf_list <- get_ecdf_list(
      flowsom_out = flowsom_out,
      chnl_sel = chnl_sel,
      chnl_lab = chnl_lab
    )
  }

  # get order of markers along rows and clusters along columns
  if (is.null(order_list)) {
    order_list <- get_order_list(
      flowsom_out = flowsom_out,
      ecdf_list = ecdf_list,
      chnl_lab = chnl_lab,
      chnl_sel = chnl_sel
    )
  }

  marker_sel_vec <- purrr::map(ecdf_list, function(x) names(x)) %>%
    unlist() %>%
    unique()

  # get relative expression for each marker for each clusters
  flowsom_out_long_grid <- purrr::map_df(
    unique(flowsom_out_long$clust),
    function(clust) {
      flowsom_out_long_clust <- flowsom_out_long %>%
        dplyr::filter(.data$clust == .env$clust)

      if (!stim_plot %in% c("all", "all_u")) {
        flowsom_out_long_clust %<>%
          dplyr::filter(stim %in% stim_plot)
      }

      purrr::map_df(marker_sel_vec, function(marker) {
        flowsom_out_long_clust_marker <- flowsom_out_long_clust %>%
          dplyr::filter(.data$marker == .env$marker)
        med <- median(flowsom_out_long_clust_marker %>%
          dplyr::pull(expr))
        ecdf_other <- ecdf_list[[clust]][[marker]]

        flowsom_out_long_clust_marker[1, ] %>%
          dplyr::mutate(perc = ecdf_other(med) * 1e2) %>%
          dplyr::select(clust, marker, perc)
      })
    }
  )

  #
  order_vec_marker <- order_list$marker
  order_vec_cluster <- order_list$cluster
  if (stab_min > 0) {
    clust_vec_fs <- unique(flowsom_out$clust)
    order_vec_cluster <- order_vec_cluster[order_vec_cluster %in% clust_vec_fs]
  }

  # get colours for heat map
  fill_vec <- RColorBrewer::brewer.pal(11, name = col_scheme)
  if (col_scheme == "OrRd") fill_vec <- rev(fill_vec)
  if (expand_mid_col) {
    n_col <- length(fill_vec)
    med_ind <- median(seq_along(n_col))
    if (as.integer(med_ind) == med_ind) {
      fill_vec[c(med_ind - 1, med_ind + 1)] <- fill_vec[med_ind]
    }
  }

  # create base plots
  p <- ggplot(
    flowsom_out_long_grid %>%
      dplyr::filter(marker %in% marker_sel_vec) %>%
      dplyr::mutate(
        clust = factor(clust, levels = order_vec_cluster),
        marker = factor(marker, levels = order_vec_marker)
      ) %>%
      dplyr::rename(`Percentile on non-cluster cells` = perc),
    aes(y = marker, x = clust)
  ) +
    geom_raster(aes(fill = `Percentile on non-cluster cells`))


  stability_tbl <- flowsom_out_long %>%
    dplyr::group_by(clust) %>%
    dplyr::slice(1) %>%
    dplyr::select(clust, stability) %>%
    dplyr::ungroup()

  stability_lab_vec <- setNames(
    round(stability_tbl$stability, 2), stability_tbl$clust
  )
  if (stability_show) {
    x_text_lab_vec <- setNames(
      paste0(
        order_vec_cluster, "\n(", stability_lab_vec[order_vec_cluster], ")"
      ),
      order_vec_cluster
    )
    p <- p +
      scale_x_discrete(labels = x_text_lab_vec)
  }

  p <- p + labs(x = "Cluster", y = "Marker")

  # if colours are evenly spaced
  if (even_col) {
    p <- p +
      scale_fill_gradientn(
        colours = rev(fill_vec),
        name = "Relative\nexpression",
        limits = c(0, 100)
      )
  } else {
    # if colours around 50 are expanded to be set to 30 and 70
    p <- p +
      scale_fill_gradientn(
        colours = rev(fill_vec),
        name = "Relative\nexpression",
        limits = c(0, 100),
        values = c(
          0, 7.5, 15, 22.5, 30,
          50,
          70, 77.5, 85, 92.5, 100
        ) / 1e2
      )
  }
  p
}
