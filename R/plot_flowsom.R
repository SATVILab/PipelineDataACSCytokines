#' @title Plot FlowSOM results
#'
#'
#' @param responders_list list.
#' Output from\code{get_post_probs} function for the current combination
#' of excluded cytokine combinations, gate method and
#' stimulation (all, all_u or a specific stim).
#' @param chnl_lab named character vector.
#' If not \code{NULL}, channels are renamed to markers using it.
#' Default is \code{NULL}
#' @param stats_combn_tbl dataframe.
#' Contains total classified cells for each sample
#' (SubjectID and VisitType) for
#' each stimulation in a column called \code{n_cell_stim}
#' and also has the total number of classified cells
#' in the corresponding unstim sample in a column
#' called \code{n_cell_uns}.
#' @param uns_chr character.
#' Name for unstim in \code{stim_vec_u}.
#' Default is "uns".
#' @param reuse,dir_save_cluster_results logical and character, respectively.
#' If \code{reuse == TRUE} and \code{!is.null(dir_save_cluster_results)},
#' then if flowsom_freq is the same as flowsom_freq found
#' @export
plot_flowsom <- function(flowsom_out,
                         flowsom_freq,
                         flowsom_cluster_post_probs,
                         stim,
                         dir_save_base,
                         chnl_lab,
                         stats_combn_tbl,
                         uns_chr = "uns",
                         reuse = FALSE,
                         dir_save_cluster_results = NULL,
                         chnl_sel) {
  if (reuse && !is.null(dir_save_cluster_results)) {
    path_fs_freq_old <- file.path(
      dir_save_cluster_results, "fs_freq_old.rds"
    )
    if (file.exists(path_fs_freq_old)) {
      flowsom_freq_old <- try(readRDS(path_fs_freq_old))
      if (identical(flowsom_freq_old, flowsom_freq)) {
        print(
          "Skipping fs plots as reuse == TRUE
          and flowsom_freq has not changed."
        )
        return(invisible(TRUE))
      }
    }
  }

  flowsom_out <- flowsom_out |>
    dplyr::mutate(stim = .data$stim |>
      str_remove("all_u-") |>
      str_remove("all-"))

  par_list <- attr(flowsom_out, "par")
  dir_save_base <- file.path(
    dir_save_base,
    par_list$scale,
    par_list$n_clust,
    paste0(stim, collapse = "_")
  )
  if (!dir.exists(dir_save_base)) {
    dir.create(dir_save_base, recursive = TRUE)
  }

  path_flowsom_out <- file.path(
    dir_save_base,
    "data_rds",
    "flowsom_out.rds"
  )
  if (!dir.exists(dirname(path_flowsom_out))) {
    dir.create(dirname(path_flowsom_out))
  }
  saveRDS(
    flowsom_out,
    file = path_flowsom_out
  )

  # ------------------
  # Preparation
  # ------------------

  # get stims to plot
  stim_vec_plot <- unique(c(as.character(flowsom_out$stim), stim)) |>
    setdiff(uns_chr)

  # ------------------
  # create plots for each stimulation
  # ------------------

  purrr::walk(stim_vec_plot, function(stim_curr) {
    # get directory to save to
    # -----------

    # get how to recognise stim in path
    stim_add <- ifelse(stim_curr %in% c("all", "all_u"), paste0(stim, "-all"),
      paste0(stim, "-", stim_curr)
    )

    # get overall directory
    dir_save_stim <- file.path(
      dir_save_base,
      stim_add
    )

    if (!dir.exists(dir_save_stim)) {
      dir.create(dir_save_stim, recursive = TRUE)
    }
    # choose which stimulations to plot
    stim_plot <- switch(stim_curr,
      "all_u" = setdiff(stim_vec_plot, "all_u"),
      "all" = setdiff(stim_vec_plot, "all"),
      stim_curr
    )

    # do dimension reduction plots only
    # if the current stim is all/all_u or if no stim is all or all_u
    plot_dim_red <- stim_curr %in% c("all", "all_u") ||
      !any(stim %in% c("all", "all_u"))
    plot_dim_red <- TRUE

    # add unstim if available

    # summarise and plot flowsom results for selected stim
    .plot_flowsom(
      flowsom_out = flowsom_out,
      path_flowsom_out = path_flowsom_out,
      flowsom_freq = flowsom_freq,
      dir_save = dir_save_stim,
      chnl_lab = chnl_lab,
      chnl_sel = chnl_sel,
      stim_plot = stim_plot,
      dim_red = plot_dim_red
    )
  })

  invisible(TRUE)
}



#' @title Plot FlowSOM results for individual stim
#' @param cell_count_lab named numeric vector.
#' Each name indicates sample id and stimulation
#' (paste0(SubjectID, "_", VisitType,
#' "_", stim)), and each element is the number of classified
#' cells for that combination of sample id and stimulation.
#' @param stim_plot character vector.
#' Stimulations to plot.
#' If "uns" and other stimulations are also in flowsom_out$stim,
#' then it will be plotted
#' where appropriate as well.
#' @param stim_clust character: "all", "all_u" or individual stim.
#' Specifies stimulations that were used to cluster.
#' @export
.plot_flowsom <- function(flowsom_out,
                          path_flowsom_out,
                          flowsom_freq,
                          dir_save,
                          chnl_lab,
                          chnl_sel,
                          stim_plot,
                          dim_red = TRUE) {
  # =======================
  # prep
  # =======================

  # identify specific stimulations needed
  stim_filter_vec <- switch(as.character(any(stim_plot %in% c("all", "all_u"))),
    "TRUE" = setdiff(unique(flowsom_out$stim), "uns"),
    "FALSE" = stim_plot
  )

  dir_save_rds <- file.path(
    dir_save, "data_rds"
  )
  if (!dir.exists(dir_save_rds)) {
    dir.create(dir_save_rds, recursive = TRUE)
  }

  # save plotting data to refer to later
  flowsom_freq_all <- flowsom_freq
  path_fs_freq_all <- file.path(
    dir_save_rds,
    "flowsom_freq_all.rds"
  )
  saveRDS(
    flowsom_freq_all,
    path_fs_freq_all
  )

  flowsom_freq_sel <- flowsom_freq |>
    dplyr::filter(stim %in% c(stim_filter_vec, "uns"))
  path_fs_freq_sel <- file.path(
    dir_save_rds, "flowsom_freq_sel.rds"
  )
  saveRDS(
    flowsom_freq_sel,
    path_fs_freq_sel
  )

  flowsom_freq_no_uns <- flowsom_freq_sel |>
    dplyr::filter(.data$stim != "uns")
  path_fs_freq_no_uns <- file.path(
    dir_save_rds, "flowsom_freq_no_uns.rds"
  )
  saveRDS(
    flowsom_freq_no_uns,
    path_fs_freq_no_uns
  )

  # colours
  # ---------------------
  stim_col_vec <- c(
    "uns" = "#a6cee3",
    "ebv" = "#1f78b4",
    "mtbaux" = "#b2df8a",
    "mtb" = "#b2df8a",
    "p1" = "#33a02c",
    "all" = "red"
  )
  prog_col_vec <- c(
    "yes" = "orange",
    "no" = "dodgerblue"
  )

  col_list <- list(
    "stim" = stim_col_vec,
    "prog" = prog_col_vec
  )

  # labels
  # ---------------------

  #' @title Labelling vector for stimulations
  stim_lab_vec_full <- c(
    "mtb" = "Live Mtb",
    "p1" = "Secreted Mtb proteins",
    "p4" = "Non-secreted Mtb proteins",
    "ebv" = "EBV and CMV",
    "uns" = "Unstim"
  )

  #' @title Labelling vector for stimulations
  stim_lab_vec_full_break <- c(
    "mtb" = "Live Mtb",
    "p1" = "Secreted\nMtb proteins",
    "p4" = "Non-secreted\nMtb proteins",
    "ebv" = "EBV and CMV",
    "uns" = "Unstim"
  )

  #' @title Labelling vector for stimulations
  stim_lab_vec_full_break_ebv <- c(
    "mtb" = "Live Mtb",
    "p1" = "Secreted\nMtb proteins",
    "p4" = "Non-secreted\nMtb proteins",
    "ebv" = "EBV and\nCMV",
    "uns" = "Unstim"
  )

  #' @title Labelling vector for stimulations when saving
  stim_lab_vec_full_save <- c(
    "mtb" = "Live Mtb",
    "p1" = "Secreted Mtb proteins",
    "p4" = "Non-secreted Mtb proteins",
    "ebv" = "EBV and CMV",
    "uns" = "Unstim"
  )

  prog_lab_vec <- c(
    "yes" = "Progressor",
    "no" = "Non-progressor"
  )

  lab_list <- list(
    "prog" = prog_lab_vec,
    "stim" = stim_lab_vec_full
  )

  # =======================
  # prep
  # =======================

  # frequencies
  # -----------------------

  # plot for unstim as well
  plot_freq(
    data = path_fs_freq_sel,
    data_no_uns = path_fs_freq_no_uns,
    dir_save = file.path(dir_save, "freq"),
    fill_col = col_list,
    fill_lab = lab_list
  )

  # =============================
  # Plot marker density by cluster
  # =============================

  # plot only for stim to plot
  plot_marker_density_by_cluster(
    flowsom_out = flowsom_out,
    path_flowsom_out = path_flowsom_out,
    chnl_lab = chnl_lab,
    dir_save = file.path(dir_save, "dens"),
    stim_plot = stim_plot,
    font_size = 10
  )

  # plot t-sne and umap
  # -----------------------

  if (dim_red) {
    plot_dim_red(
      flowsom_out = flowsom_out,
      path_flowsom_out = path_flowsom_out,
      dir_save = file.path(dir_save, "dim_red"),
      stim_plot = stim_plot,
      col_clust = switch(as.character(length(stim_plot) > 1),
        "TRUE" = col_list[["stim"]]["all"],
        "FALSE" = col_list[["stim"]][stim_plot]
      ),
      font_size = 10,
      width = 6,
      height = 6
    )
  }

  # =============================
  # Plot heat maps
  # =============================

  plot_clust_map(
    flowsom_out = flowsom_out,
    path_flowsom_out = path_flowsom_out,
    dir_save = file.path(dir_save, "h_map"),
    stim_plot = stim_plot,
    even_col = TRUE,
    stability_min = c(0, 0.5),
    expand_mid_col = TRUE,
    stability_show = TRUE,
    chnl_lab = chnl_lab
  )

  invisible(TRUE)
}

# TODO: this function
plot_fs_post_probs <- function(fs_post_probs) {
  1
}
