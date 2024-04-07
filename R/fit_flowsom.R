#' @title Fit FlowSOM to flowSet
#'
#' @description Also calculates cluster stability and extracts sample-specific information.
#'
#' @param fs flowSet.
#' @param n_clust "n_m" or "n_f". If \code{"n_m"}, then 15 clusters are used.
#' If \code{"n_f"}, then 12 clusters are used for CD4 and CD8 datasets but only 7 for
#' the TCRgd dataset.
#' @param scale 'scale' or 'n_s'. If \code{scale}, then columns are centred and standardised.
#' @param chnl_sel character vector. Channels to apply FlowSOM to.
#' @param ds_name character. Name of dataset (as used to name datasets in output of data package).
#' @param dir_save_base character. Base directory to save to.
#' @param jacc_n Number of bootstrapped samples to apply cluster algorithm to to calculate cluster stability. Default is 0.
#' @param jacc_reuse logical. If \code{TRUE}, then old data is reused if available. Default is \code{FALSE}.
#' @param fs_seed,jacc_seed integer. Seeds to use for the FlowSOM and cluster stability algorithms.
#' @param chnl_lab_vec named character vector. If not \code{NULL}, channels are renamed to markers using it.
#' Default is \code{NULL}.
#' @param tsne_tbl dataframe. Has columns SubjectID, VisitType, tsne1, tsne2, umap1 and umap2.
#' @param stim 'all_u', 'all', 'p1', 'p4', 'ebv' or 'mtbaux'. Specifies stimulations.
#'
#' @return A dataframe with cells along rows and cell-level information ranging from sample information to
#' cluster and tsne information to expression information.
#' @export
fit_flowsom <- function(fs, n_clust, scale, chnl_sel, ds_name, dir_save_base,
                        jacc_reuse = FALSE, jacc_n = 0, fs_seed = NULL, jacc_seed = NULL,
                        chnl_lab_vec, tsne_tbl, stim) {
  # ===========================================
  # cluster all data together using FlowSOM
  # ===========================================

  nClus <- switch(n_clust,
    "n_m" = 15,
    "n_f" = purrr::map_dbl(ds_name, function(x) {
      if (str_detect(ds_name, "cd4|cd8")) {
        return(12)
      }
      if (str_detect(ds_name, "tcrgd")) {
        return(7)
      }
      if (str_detect(ds_name, "bcell")) {
        return(6)
      }
      if (str_detect(ds_name, "nk")) {
        return(10)
      }
      12
    }),
    "n_f_f" = purrr::map_dbl(ds_name, function(x) {
      if (str_detect(ds_name, "cd4|cd8")) {
        return(8)
      }
      if (str_detect(ds_name, "tcrgd")) {
        return(5)
      }
      if (str_detect(ds_name, "nk")) {
        return(6)
      }
      if (str_detect(ds_name, "bcell")) {
        return(3)
      }
      8
    }),
    "n_meta" = NULL,
    stop("n_clust value not either n_m or n_f or n_meta")
  )

  # browser()
  # debugonce(FlowS)
  flowsom <- FlowSOM::FlowSOM(
    input = fs,
    transform = FALSE,
    colsToUse = chnl_sel,
    nClus = nClus,
    maxMeta = 20,
    scale = scale == "scale",
    seed = fs_seed
  )


  # ===========================================
  # Calculate cluster stability
  # ===========================================

  # directory to save stability to/read from
  stim_add <- switch(as.character(str_detect(stim, "all")),
    "TRUE" = paste0(stim, "-", "all"),
    "FALSE" = stim
  )
  dir_save_jacc <- file.path(dir_save_base, stim_add, scale, n_clust)

  # get old cluster stability results if available
  # --------------------------

  # if allowed to reuse
  if (jacc_reuse) {
    # try to load
    jacc_tbl <- tryCatch(readRDS(file.path(dir_save_jacc, "jacc_tbl.rds")),
      error = function(cond) NULL
    )
    jacc_lab_vec <- tryCatch(readRDS(file.path(dir_save_jacc, "jacc_lab_vec.rds")),
      error = function(cond) NULL
    )
    # set to NULL if number of previous samples was less than current
    if (!is.null(jacc_lab_vec)) {
      # check that we don't have too few samples and that we have
      # the right number of clusters
      if (nrow(jacc_tbl) < jacc_n || ncol(jacc_tbl) != (nClus + 1)) jacc_lab_vec <- NULL
    }
  } else {
    jacc_lab_vec <- NULL
  } # set to NULL if not available

  # calculate again if no old data available and number of new boot sample is greater than 0
  # ---------------------------
  if (is.null(jacc_lab_vec) && jacc_n > 0) {
    # get jacc coef for each sample
    jacc_tbl <- calc_jacc(
      fs = fs,
      fsom_orig = flowsom,
      seed = jacc_seed,
      boot = jacc_n,
      colsToUse = chnl_sel,
      nClus = nClus,
      stim_and_prog_equal = TRUE,
      sample_by = "cell and sample",
      scale = scale == "scale"
    )

    # get mean over all samples
    jacc_tbl_mean <- jacc_tbl |>
      tidyr::pivot_longer(-boot,
        names_to = "cluster",
        values_to = "jacc"
      ) |>
      dplyr::group_by(cluster) |>
      dplyr::summarise(jacc_mean = mean(jacc[jacc >= 0]))

    # get labelling vector for matching clusters to stability
    jacc_lab_vec <- stats::setNames(
      jacc_tbl_mean$jacc_mean,
      jacc_tbl_mean$cluster
    )

    # save scores
    analysispipeline::save_objects(
      jacc_tbl_mean = jacc_tbl_mean,
      jacc_tbl = jacc_tbl,
      jacc_lab_vec = jacc_lab_vec,
      dir_proj = dir_save_jacc,
      empty = FALSE
    )
  }


  # =================================================
  # Merge flowsom output with other information
  # =================================================

  # preparation of object itself
  # --------------------

  # relevant data
  mapping <- flowsom$FlowSOM$map$mapping
  meta <- flowsom$FlowSOM$metaData

  # joined table
  flowsom_out <- flowsom$FlowSOM$data |>
    tibble::as_tibble() |>
    dplyr::bind_cols(tibble::tibble(clust = flowsom$metaclustering[flowsom$FlowSOM$map$mapping[, 1]] |>
      as.character()))

  # rescale
  if (scale == "scale") {
    for (chnl in chnl_sel) {
      flowsom_out[, chnl] <- flowsom_out[, chnl] * flowsom$FlowSOM$scaled.scale[chnl] +
        flowsom$FlowSOM$scaled.center[chnl]
    }
  }

  # get sample-specific information from file names
  # -----------

  # get file name for each cell
  flowsom_out <- flowsom_out |> dplyr::mutate(fcs = NA_character_)
  for (i in seq_along(meta)) {
    fn_long <- names(meta)[i]
    slash_loc_tbl <- str_locate_all(fn_long, "/")[[1]]
    last_slash_loc <- slash_loc_tbl[, "end"][[nrow(slash_loc_tbl)]]
    fn_short <- str_sub(fn_long, start = last_slash_loc + 1)
    ind_vec <- seq(meta[[i]][1], meta[[i]][2])
    flowsom_out[ind_vec, "fcs"] <- fn_short
  }

  # get file names
  fcs_vec <- flowsom_out$fcs

  # get sample info from each file name
  sample_info_tbl <- purrr::map_df(unique(fcs_vec), function(fcs) {
    split_vec <- str_split(str_sub(fcs, end = -5), "_")[[1]]
    tibble::tibble(
      SubjectID = split_vec[1],
      VisitType = split_vec[2],
      Progressor = split_vec[4],
      timeToTB = split_vec[3],
      stim = split_vec[5],
      fcs = fcs
    )
  })

  # bind onto flowsom data
  flowsom_out <- flowsom_out |> dplyr::left_join(sample_info_tbl, by = c("fcs"))

  # rename columns in flowsom_out
  # ---------------------------

  # only do this if renaming vector is supplied
  if (!is.null(chnl_lab_vec)) {
    for (chnl in names(chnl_lab_vec)) {
      if (chnl %in% colnames(flowsom_out)) {
        colnames(flowsom_out)[which(colnames(flowsom_out) == chnl)] <- chnl_lab_vec[chnl]
      }
    }
    # remove CD45
    flowsom_out <- flowsom_out[, which(colnames(flowsom_out) != "CD45")]
  }

  # add cluster stability, if available
  # -------------------------

  if (!is.null(jacc_lab_vec)) {
    flowsom_out <- flowsom_out |>
      dplyr::mutate(stability = jacc_lab_vec[clust])
  }

  # add tsne coords, if available
  # -----------------------------

  # do this separately for each stimulation
  flowsom_out <- purrr::map_df(unique(flowsom_out$stim), function(stim_tsne) {
    # subset tsne tbl
    tsne_tbl_curr <- tsne_tbl |> dplyr::filter(.data$stim == .env$stim_tsne)

    # subset flowsom tbl
    flowsom_out <- flowsom_out[, !colnames(flowsom_out) == "CD45"] |>
      dplyr::filter(stim == stim_tsne)

    # create cell_ind column in flowsom_out
    flowsom_out <- flowsom_out |>
      dplyr::group_by(SubjectID, VisitType, stim) |>
      dplyr::mutate(cell_ind = 1:n()) |>
      dplyr::ungroup()

    # create cell_ind column in tsne_tbl
    tsne_tbl_curr <- tsne_tbl_curr |>
      dplyr::group_by(SubjectID, VisitType) |>
      dplyr::mutate(cell_ind = 1:n()) |>
      dplyr::ungroup()

    # bind flowsom and tsne tbls together (umap assumed here)
    flowsom_out <- flowsom_out |>
      dplyr::left_join(
        tsne_tbl_curr |>
          dplyr::select(SubjectID, VisitType, cell_ind, tsne1, tsne2, umap1, umap2),
        by = c("SubjectID", "VisitType", "cell_ind")
      )

    # return
    flowsom_out
  })

  # Add relevant attributes
  # ------------------------

  # add scaling factors if available
  attr(flowsom_out, "scale") <- switch(as.character(scale == "scale"),
    "TRUE" = flowsom$FlowSOM$scaled.scale,
    "FALSE" = NULL
  )

  # add scale and n_clust parameters
  attr(flowsom_out, "par_fs") <- list(
    "n_clust" = n_clust, "scale" = scale,
    "chnl_sel" = chnl_sel
  )

  # make cluster and stim variables factors
  # -----------------------

  # convert mtb to mtbaux if mtbaux is what's in flowsom_out
  if ("mtbaux" %in% flowsom_out$stim) {
    flowsom_out <- flowsom_out |> dplyr::mutate(stim = ifelse(stim == "mtbaux", "mtb", stim))
  }

  # get order vector (assuming unstim is available)
  stim_order_vec <- switch(stim,
    "all_u" = c("uns", "ebv", "mtb", "p1"),
    "all" = c("ebv", "mtb", "p1"),
    c("uns", stim)
  )

  # remove unstim if it is available
  if (!"uns" %in% flowsom_out$stim) {
    stim_order_vec <- stim_order_vec |> setdiff("uns")
  }

  # make stim and cluster factors
  stim_factor_vec <- factor(flowsom_out$stim, stim_order_vec)

  flowsom_out <- flowsom_out |>
    dplyr::mutate(
      stim = stim_factor_vec,
      clust = factor(as.character(clust), levels = as.character(sort(as.numeric(unique(clust)))))
    )

  # arrange
  flowsom_out <- flowsom_out |> dplyr::arrange(SubjectID, VisitType, stim, clust, cell_ind)

  flowsom_out
}
