#' @title Get posterior probability of a response for each cluster
#'
#' @param reuse logical. If \code{TRUE} (and dir_save != NULL), then previous output will be reused if found in
#' \code{dir_save}.
#' @param dir_save_cluster_results character. Directory to save output to, if created, and to get past output from (if reuse == TRUE).
#'
#' @details
#' If the MIMOSA fit fails for a given stim and cluster combination
#' with the EM algorithm, then MCMC will
#' be attempted. If this too fails, then a single row dataframe
#' with NAs except for stim and cluster will be returned for that stim
#' and cluster combination.
#'
#' @return
#' If there is an 'uns' stim, thena dataframe with columns clust, stim
#' SubjectID, VisitType, Progressor and prob is returned.
#' If there is no 'uns' stim, then NULL is returned.
get_fs_cluster_post_probs <- function(flowsom_freq = flowsom_freq,
                                      dir_save_cluster_results = NULL, reuse = FALSE) {
  if (!"uns" %in% flowsom_freq$stim) {
    return(NULL)
  }

  if (reuse && !is.null(dir_save_cluster_results)) {
    path_fs_freq_old <- file.path(dir_save_cluster_results, "fs_freq_old.rds")
    if (file.exists(path_fs_freq_old)) {
      flowsom_freq_old <- try(readRDS(path_fs_freq_old))
      if (identical(flowsom_freq_old, flowsom_freq)) {
        path_fs_probs_old <- file.path(dir_save_cluster_results, "fs_post_probs.rds")
        if (file.exists(path_fs_probs_old)) {
          print("Skipping fs plots as reuse == TRUE, flowsom_freq has not changed and previous fs_post_probs available.")
          return(readRDS(path_fs_probs_old))
        }
      }
    }
  }

  if (!is.null(dir_save_cluster_results)) {
    if (!dir.exists(dir_save_cluster_results)) dir.create(dir_save_cluster_results, recursive = TRUE)
    path_fn_probs_save <- file.path(dir_save_cluster_results, "fs_post_probs.rds")
  }


  fs_freq_uns <- flowsom_freq |> dplyr::filter(stim == "uns")
  out_tbl <- purrr::map_df(setdiff(unique(flowsom_freq$stim), "uns"), function(stim) {
    fs_freq_stim <- flowsom_freq |> dplyr::filter(stim == .env$stim)
    # print(stim)
    clust_vec <- unique(fs_freq_stim$clust)
    # clust_vec <- setdiff(clust_vec, c("3"))
    clust_tbl <- purrr::map_df(clust_vec, function(clust) {
      # print(clust)
      fs_freq_stim_clust <- fs_freq_stim |> dplyr::filter(clust == .env$clust)
      fs_freq_uns_clust <- fs_freq_uns |> dplyr::filter(clust == .env$clust)
      fs_freq_clust <- try(fs_freq_stim_clust |>
        dplyr::bind_rows(fs_freq_uns_clust) |>
        dplyr::mutate(
          count_pos = count_stim,
          count_neg = n_cell_tot_stim - count_stim
        ) |>
        dplyr::mutate(feature = paste0("c", clust)) |>
        dplyr::select(SubjectID, VisitType, Progressor, feature, stim, count_pos, count_neg) |>
        dplyr::mutate(SampleID = paste0(SubjectID, "_", VisitType)) |>
        dplyr::select(-c(SubjectID, VisitType)))

      mimosa_expr_set <- try(MIMOSA::ConstructMIMOSAExpressionSet(
        thisdata = fs_freq_clust,
        reference = stim %in% "uns",
        measure.columns = c("count_pos", "count_neg"),
        other.annotations = c("feature"),
        .variables = plyr::.(SampleID)
      ))

      mim_res <- try(MIMOSA::MIMOSA(
        formula = count_neg + count_pos ~ SampleID | feature,
        data = mimosa_expr_set,
        method = "EM"
      ))

      # get posterior probabilities
      post_prob_vec <- try(mim_res[[1]]@z[, "z2"])
      if (class(post_prob_vec)[1] == "try-error") {
        return(tibble::tibble(
          clust = clust,
          stim = stim,
          SubjectID = NA,
          VisitType = NA,
          Progressor = NA,
          prob = NA
        ))
      }

      # get sample id and subject id
      sampleid_vec <- str_sub(names(post_prob_vec), 11) |>
        str_remove("_count") |>
        str_remove(paste0("_c", clust))
      subjectid_vec <- str_sub(sampleid_vec, end = 6)
      visittype_vec <- str_sub(sampleid_vec, start = 8)
      prog_lab_tbl <- flowsom_freq |>
        # mutate(SampleID = paste0(SubjectID, '_', VisitType)) |>
        dplyr::group_by(SubjectID) |>
        dplyr::slice(1)
      prog_lab_vec <- setNames(prog_lab_tbl$Progressor, prog_lab_tbl$SubjectID)

      # get named vector of posterior probabilities per sample
      post_prob_lab_vec <- setNames(post_prob_vec, paste0(subjectid_vec, "_", visittype_vec))

      # get sampleids with post probs that exceed minimum
      # r_o_sampleid_vec <- names(post_prob_lab_vec)[post_prob_lab_vec > post_prob_min]

      out_tbl_ind <- tibble::tibble(
        clust = clust, stim = stim, SubjectID = subjectid_vec, VisitType = visittype_vec, Progressor = prog_lab_vec[subjectid_vec],
        prob = post_prob_lab_vec
      )

      if (all(is.na(out_tbl_ind$prob))) {
        return(tibble::tibble(
          clust = clust,
          stim = stim,
          SubjectID = NA,
          VisitType = NA,
          Progressor = NA,
          prob = NA
        ))
      }
      out_tbl_ind
    })

    clust_tbl
  })

  if (!is.null(dir_save_cluster_results)) saveRDS(out_tbl, path_fn_probs_save)

  out_tbl
}
