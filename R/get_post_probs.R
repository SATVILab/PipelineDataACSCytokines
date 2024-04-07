#' @title Get posterior probabilities of response for each sample for each individual stim and per-stim list of responders
#'
#' @description Uses MIMOSA.
#'
#' @param stats_combn_tbl df.
#' @param stim 'all_u', 'all', 'p1', 'p4', 'ebv' or 'mtbaux'. Stimulations for which posterior probabilities should
#' be calculated
#' @param post_prob_min 0 to 1. Minimum posterior probability to count as a responder.
#' @param chnl_lab_cyt_custom named character vector. Converts cytokine chnl names to custom marker names.
#' If \code{NULL}, then preset labelling vector is used. Default is \code{NULL}.
#' @param chnl_lab_cyt_fcs named character vector. Converts cytokine chnl names to FCS-derived marker names. No default, and
#' must be supplied.
#' @param stim_vec_u character vector. Consists of names of individual stimulations. If \code{NULL}, then
#' set to default value. Default is \code{NULL}.
#' @param uns_chr character. Name for unstim in \code{stim_vec_u}. Default is "uns".
#' @param conv_mtbaux_to_mtb logical. If \code{TRUE}, then when looping over stimulations to
#' calculate posterior probabilities, "mtbaux" is converted to "mtb" in the loop index.
#' @param gn character. Name of gate in \code{stats_combn_tbl}. If \code{NULL}, then set to "locb0.15_min_clust".
#' Default is \code{NULL}.
#' @param exc character vector. Contains channel names of cytokine combination(s) separated by underscores ("_")
#' that are not included as part of the summed response of cytokine-positive cells.
#'
#' @return List with one named element for each stimulation, each containing
#' a list with named elements 'post_prob' and 'responders', which are the per-sample posterior probabilities
#' and the responding samples, respectively.
#' @export
get_post_probs <- function(stats_combn_tbl, stim, exc,
                           uns_chr = "uns", conv_mtbaux_to_mtb = TRUE,
                           gn = NULL, chnl_lab_cyt_custom = NULL, chnl_lab_cyt_fcs,
                           post_prob_min) {
  # ==============================
  # Preparation
  # ==============================

  # named character vector for converting chnl names to custom (i.e. non-FCS-derived) marker names
  if (is.null(chnl_lab_cyt_custom)) {
    chnl_lab_cyt_custom <- c(
      "Ho165Di" = "IFNg", "Gd158Di" = "IL2",
      "Nd146Di" = "TNF", "Dy164Di" = "IL17",
      "Nd150Di" = "IL22", "Gd156Di" = "IL6"
    )
  }
  # gn
  if (is.null(gn)) gn <- "locb0.15_min_clust"

  if (conv_mtbaux_to_mtb && any(stim == "mtbaux")) stim[which(stim == "mtbaux")] <- "mtb"

  # ==============================
  # Get posterior probabilities and responders for each stim
  # ==============================
  out_list <- purrr::map(stim, function(stim_ind) {
    # ---------------------------
    # prepare stats_combn_tbl
    # ---------------------------

    # filter out for selected stim and gate
    ds <- stats_combn_tbl |>
      dplyr::filter(
        .data$stim == stim_ind,
        gate_name == gn
      )

    # remove pre-specified cytokine combinations
    # ---------------------------

    # identify what pre-specified cytokine combinations are
    exc_cyt_single_vec <- purrr::map(stringr::str_split(exc, "_")[[1]], function(x) {
      chnl_vec <- names(chnl_lab_cyt_fcs)
      for (chnl in chnl_vec) {
        if (str_detect(x, chnl)) {
          return(chnl)
        }
      }
      NULL
    }) |>
      purrr::compact() |>
      unlist()

    cyt_vec <- stats_combn_tbl$cyt_combn[[1]] |>
      stringr::str_remove_all("[[!]]") |>
      stringr::str_split("&") |>
      magrittr::extract2(1)

    cyt_lab_cyt_vec <- stats::setNames(
      names(chnl_lab_cyt_custom),
      chnl_lab_cyt_custom
    )

    chnl_vec <- cyt_lab_cyt_vec[cyt_vec]

    # remove pre-specified cytokine combinations
    exc_cyt_combn_vec <- purrr::map_chr(exc_cyt_single_vec, function(x) {
      if (x == chnl_vec[1]) {
        return(paste0(
          chnl_vec[1], "&!",
          paste0(cyt_vec[-1], collapse = "&!")
        ))
      }
      cyt_combn <- paste0("!", cyt_vec[[1]])
      for (i in seq_along(chnl_vec)[-1]) {
        if (chnl_vec[i] == x) {
          cyt_combn <- paste0(cyt_combn, "&", cyt_vec[i])
        } else {
          cyt_combn <- paste0(cyt_combn, "&!", cyt_vec[i])
        }
      }
      cyt_combn
    })


    for (cyt_combn in exc_cyt_combn_vec) {
      ds <- ds |> dplyr::filter(
        !.data$cyt_combn == .env$cyt_combn
      )
    }

    # remove all-negative cytokine combination
    # ---------------------------

    # get all-negative cytokine combination in correct format
    all_neg_cyt_combn <- paste0("!", paste0(cyt_vec, collapse = "&!"))

    # remove all-negative cytokine combination
    ds <- ds |> dplyr::filter(!.data$cyt_combn == all_neg_cyt_combn)

    # create SampleID column
    ds <- ds |>
      dplyr::mutate(SampleID = paste0(SubjectID, "_", VisitType))

    # prepare count data for MIMOSA
    # ---------------------------

    # calculate summed counts
    ds <- ds |>
      group_by(
        gate_name, SampleID, SubjectID, VisitType, stim, n_cell_stim,
        n_cell_uns
      ) |>
      dplyr::summarise(
        count_stim = sum(count_stim),
        count_uns = sum(count_uns),
        .groups = "drop"
      )

    # calculate count_pos and count_neg columns
    ds <- ds |>
      dplyr::rename(stim_old = stim) |>
      tidyr::pivot_longer(n_cell_stim:n_cell_uns,
        names_to = "stim",
        values_to = "total_cells_count"
      ) |>
      dplyr::mutate(stim = ifelse(str_detect(stim, "stim"), stim_old, "uns")) |>
      dplyr::select(-stim_old) |>
      dplyr::mutate(
        count_neg = ifelse(stim == "uns",
          total_cells_count - count_uns,
          total_cells_count - count_stim
        ),
        count_pos = ifelse(stim == "uns",
          count_uns,
          count_stim
        )
      ) |>
      dplyr::select(-c(count_stim, count_uns, total_cells_count))

    # create MIMOSA-specific columns
    # ---------------------------

    # add feature variable (neeed for MIMOSA)
    ds <- ds |>
      dplyr::mutate(feature = "any_pos_but_exc") |>
      dplyr::select(feature, everything())

    # create ag column
    ds <- ds |>
      dplyr::mutate(ag = setdiff(stim, uns_chr))

    # create MIMOSA expression set
    # ---------------------------

    # browser()
    # library(MIMOSA)
    mimosa_expr_set <- MIMOSA::ConstructMIMOSAExpressionSet(
      thisdata = ds,
      reference = stim %in% "uns",
      measure.columns = c("count_pos", "count_neg"),
      other.annotations = c("gate_name", "SubjectID", "VisitType", "feature"),
      .variables = plyr::.(SampleID)
    )

    # ---------------------------
    # Run MIMOSA
    # ---------------------------


    mim_res <- MIMOSA::MIMOSA(
      formula = count_neg + count_pos ~ SampleID | feature,
      data = mimosa_expr_set,
      method = "EM"
    )

    # detach('package:MIMOSA', unload = TRUE)
    # detach('package:MASS', unload = TRUE)
    # detach('package:plyr', unload = TRUE)
    # detach('package:reshape', unload = TRUE)
    # detach('package:Biobase', unload = TRUE)

    # ---------------------------
    # Extract results from MIMOSA
    # ---------------------------

    # get posterior probabilities
    mim_res_ft <- mim_res$`any_pos_but_exc`
    post_prob_vec <- mim_res_ft@z[, 1]

    # get sample id and subject id
    sampleid_vec <- str_sub(names(post_prob_vec), 11)
    subjectid_vec <- str_sub(sampleid_vec, end = 6)
    sampleid_vec <- str_sub(sampleid_vec, start = 8)
    slash_loc_tbl <- str_locate(sampleid_vec, "_")[, 1]
    visittype_vec <- str_sub(sampleid_vec, end = slash_loc_tbl - 1)

    # get named vector of posterior probabilities per sample
    post_prob_lab_vec <- stats::setNames(post_prob_vec, paste0(subjectid_vec, "_", visittype_vec))

    # get sampleids with post probs that exceed minimum
    r_o_sampleid_vec <- names(post_prob_lab_vec)[post_prob_lab_vec > post_prob_min]

    list(
      "post_prob" = post_prob_lab_vec,
      "responders" = r_o_sampleid_vec
    )
  }) |>
    stats::setNames(stim)

  # add minimum posterior probability to count as a responder here
  attr(out_list, "pars") <- list("post_prob_min" = post_prob_min)

  out_list
}
