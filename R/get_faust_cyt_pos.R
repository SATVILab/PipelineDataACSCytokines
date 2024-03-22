#' @export
get_faust_cyt_pos <- function(dir_faust,
                              dir_save,
                              pop_to_defn, # pop_list
                              pop_to_gate, # pop_to_gate_tbl
                              chnl,
                              force = FALSE) {
  if (!dir.exists(dir_save)) dir.create(dir_save)

  # faust-specific information
  # ======================
  fcs_vec_faust <- list.dirs(
    file.path(dir_faust, "sampleData"),
    recursive = FALSE,
    full.names = FALSE
  )

  # helps go from batch id to faust directory containing
  # annotations for each sample
  tidy_to_faust_fcs <- setNames(
    fcs_vec_faust,
    DataTidyACSCyTOFFAUST::clean_fcs_for_matching(fcs_vec_faust)
  )

  for (pop in names(pop_to_defn)) {
    dir_save_pop <- file.path(
      dir_save,
      pop
    )
    if (force) {
      unlink(dir_save_pop, recursive = TRUE)
    }
    if (!dir.exists(dir_save_pop)) dir.create(dir_save_pop)

    # Non-FAUST prep
    # =====================

    # get pop definition
    # ---------------------
    pop_defn <- pop_to_defn[[pop]]

    # get gate table
    # ----------------------
    gate_tbl <- pop_to_gate[[pop]]

    # get expression data
    # ----------------------

    dir_gs <- normalizePath(
      file.path(
        dir_faust,
        "gsData",
        paste0("gs_cytof_acs_", gsub("_dim|_bright", "", pop))
      )
    )
    gs <- flowWorkspace::load_gs(dir_gs)

    # get ind_batch_list
    # ----------------------
    pop_gate <- "root"
    ind_in_batch_lab_vec <- c(
      "1" = "ebv", "2" = "mtb", "3" = "p1",
      "4" = "p4", "5" = "uns"
    )
    ind_in_batch_uns <- 5

    fcs_stim_tbl <- purrr::map_df(seq_along(gs), function(i) {
      fn <- flowCore::description(
        flowWorkspace::gh_pop_get_data(gs[[i]])
      )$FILENAME
      last_slash_loc_tbl <- stringr::str_locate_all(fn, "/")[[1]]
      last_slash_loc <- last_slash_loc_tbl[nrow(last_slash_loc_tbl), "end"][[1]]
      fn <- stringr::str_sub(fn, last_slash_loc + 1) |>
        stringr::str_remove(".fcs")
      second_underscore_loc <- stringr::str_locate_all(fn, "_")[[1]][2, "end"][[1]]
      fcs <- stringr::str_sub(fn, end = second_underscore_loc - 1)
      stim <- stringr::str_split(fn, "_")[[1]]
      stim <- stim[[length(stim)]]
      tibble::tibble(fcs = fcs, stim = stim)
    })

    ind_batch_list <- outliergatev1:::.get_ind_batch_list(
      data = gs, batch_size = 5,
      cytof_fcs_to_clin_map = NULL,
      fcs = fcs_stim_tbl$fcs
    )


    # get cytokine-combination information
    # ----------------------------

    n_chnl <- length(chnl)

    # this is used lower down
    combn_ind_pos_list <- purrr::map(seq_len(length(chnl)), function(n_pos) {
      gtools::combinations(n = n_chnl, r = n_pos)
    }) |>
      setNames(seq_len(length(chnl)))
    n_combn <- purrr::map_dbl(combn_ind_pos_list, nrow) |>
      sum()
    combn_nm_list <- purrr::map(names(combn_ind_pos_list), function(n_pos_nm) {
      combn_mat <- combn_ind_pos_list[[n_pos_nm]]
      purrr::map_chr(1:nrow(combn_mat), function(i) {
        chnl_pos <- cyt_vec[combn_mat[i, , drop = TRUE]]
        chnl_neg <- cyt_vec[setdiff(seq_len(n_chnl), combn_mat[i, , drop = TRUE])]
        purrr::map_chr(cyt_vec, function(chnl_curr) {
          ifelse(chnl_curr %in% chnl_pos, paste0(chnl_curr, "~+~"), paste0(chnl_curr, "~-~"))
        }) |>
          paste0(collapse = "")
      })
    }) |>
      setNames(names(combn_ind_pos_list))


    # loop over batches (just to make sure you load unstim correctly)
    batch <- names(ind_batch_list)[5]
    for (batch in names(ind_batch_list)) {
      # sample id
      subject_id <- gsub("_D\\d+", "", batch)
      visit_type <- gsub("\\w+_", "", batch)

      # indices in batch
      ind_batch <- ind_batch_list[[batch]]

      # get expression data
      ex_list <- purrr::map(ind_batch, function(ind) {
        ind_in_batch <- ifelse(ind %% 5 == 0, 5, ind - (ind %/% 5) * 5)
        flowWorkspace::gh_pop_get_data(
          gs[[ind]]
        ) |>
          flowCore::exprs() |>
          tibble::as_tibble() |>
          dplyr::mutate(
            batch = batch,
            ind = ind,
            stim = ind_in_batch_lab_vec[ind_in_batch]
          )
      })

      # get general info, inc. fcs name, for this sample
      # so that you can match it to the FAUST fcs
      # to get the FAUST annotation
      match_tbl <- DataTidyACSCyTOFFAUST::fcs_to_sampleid_map |>
        dplyr::filter(
          SubjectID == subject_id,
          VisitType == visit_type
        )

      ex_ind <- 1
      for (ex_ind in seq_along(ex_list[-ind_in_batch_uns])) {
        # if (debug) print(paste0("ex_ind: ", ex_ind))

        # get faust labels - stim
        # ============================

        # get faust fcs name from ind
        # to get faust directory containing
        # sample's faust annotation
        ind_in_batch <- ifelse(ex_ind %% 5 == 0, 5, ex_ind - (ex_ind %/% 5) * 5)
        stim <- ind_in_batch_lab_vec[ind_in_batch][[1]]
        match_tbl_stim <- match_tbl |>
          dplyr::filter(Stim == stim)
        ind_to_fcs_match <- match_tbl_stim$MatchFCSName
        faust_fcs <- tidy_to_faust_fcs[ind_to_fcs_match]

        # create directory to save to
        dir_save_pop_ind <- file.path(
          dir_save_pop,
          paste0(match_tbl_stim$SampleID, "_", match_tbl_stim$Stim)
        )
        if (!dir.exists(dir_save_pop_ind)) {
          dir.create(dir_save_pop_ind)
        }
        if (!force & file.exists(file.path(dir_save_pop_ind, "completed.rds"))) {
          next
        }
        if (file.exists(file.path(dir_save_pop_ind, "completed.rds"))) {
          file.remove(file.path(dir_save_pop_ind, "completed.rds"))
        }

        # get sample faust directory
        dir_sample <- file.path(
          dir_faust,
          "sampleData",
          faust_fcs
        )

        # get faust annotations
        ann_tbl_faust <- readr::read_csv(
          file.path(
            dir_sample,
            "faustAnnotation.csv"
          ),
          col_names = FALSE,
          col_types = "c"
        )

        # subset to only pop
        ind_list <- purrr::map(seq_along(pop_defn), function(i) {
          pttrn <- paste0(names(pop_defn)[i], "~", pop_defn[[i]])
          if (length(pttrn) == 1) {
            return(grepl(pttrn, ann_tbl_faust[[1]]))
          }
          purrr::map(pttrn, function(x) {
            grepl(x, ann_tbl_faust[[1]])
          }) |>
            purrr::reduce(
              .f = function(x, y) x | y
            )
        })

        ind_vec <- purrr::reduce(
          .x = ind_list,
          .f = function(x, y) x & y
        )

        ann_tbl_faust_pop <- ann_tbl_faust[ind_vec, ]

        n_cell <- nrow(ann_tbl_faust_pop)
        for (i in seq_along(combn_nm_list)) {
          combn_vec <- combn_nm_list[[i]]
          for (j in seq_along(combn_vec)) {
            tbl_add <- tibble::tibble(X = rep(NA_real_, n_cell))
            colnames(tbl_add) <- combn_vec[j]
            ann_tbl_faust_pop <- ann_tbl_faust_pop |>
              dplyr::bind_cols(tbl_add)
          }
        }

        # get faust labels - uns
        # ============================

        # get faust fcs name from ind
        match_tbl_uns <- match_tbl |>
          dplyr::filter(Stim == "uns")
        ind_to_fcs_match_uns <- match_tbl_uns$MatchFCSName
        faust_fcs_uns <- tidy_to_faust_fcs[ind_to_fcs_match_uns]


        # get sample faust directory
        dir_sample_uns <- file.path(
          dir_faust,
          "sampleData",
          faust_fcs_uns
        )

        # get faust annotations
        ann_tbl_faust_uns <- readr::read_csv(
          file.path(
            dir_sample_uns,
            "faustAnnotation.csv"
          ),
          col_names = FALSE,
          col_types = "c"
        )

        # subset to only pop
        ind_list_uns <- purrr::map(seq_along(pop_defn), function(i) {
          pttrn <- paste0(names(pop_defn)[i], "~", pop_defn[[i]])
          if (length(pttrn) == 1) {
            return(grepl(pttrn, ann_tbl_faust_uns[[1]]))
          }
          purrr::map(pttrn, function(x) {
            grepl(x, ann_tbl_faust_uns[[1]])
          }) |>
            purrr::reduce(
              .f = function(x, y) x | y
            )
        })

        ind_vec_uns <- purrr::reduce(
          .x = ind_list_uns,
          .f = function(x, y) x & y
        )

        ann_tbl_faust_pop_uns <- ann_tbl_faust_uns[ind_vec_uns, ]

        n_cell_uns <- nrow(ann_tbl_faust_pop_uns)
        for (i in seq_along(combn_nm_list)) {
          combn_vec <- combn_nm_list[[i]]
          for (j in seq_along(combn_vec)) {
            tbl_add <- tibble::tibble(X = rep(NA_real_, n_cell_uns))
            colnames(tbl_add) <- combn_vec[j]
            ann_tbl_faust_pop_uns <- ann_tbl_faust_pop_uns |>
              dplyr::bind_cols(tbl_add)
          }
        }

        # get expression data
        # ----------------------
        # get it from GatingSet, as
        # the FAUST expression data doesn't
        # include markers FAUST didn't use
        ex <- ex_list[-ind_in_batch_uns][[ex_ind]]
        if (nrow(ex) != nrow(ann_tbl_faust_pop)) {
          stop("rows in ex and ann_tbl_faust_pop do not match")
        }
        ex_uns <- ex_list[[ind_in_batch_uns]]
        if (nrow(ex_uns) != nrow(ann_tbl_faust_pop_uns)) {
          stop("rows in ex and ann_tbl_faust_pop do not match")
        }

        # get gate data for this individual
        # -----------------------
        gate_tbl_ind <- gate_tbl |> dplyr::filter(ind == ex$ind[1])

        n_pos_nm <- "2"
        for (n_pos_nm in names(combn_ind_pos_list)) {
          # if (debug) print(paste0("n_pos_nm: ", n_pos_nm))

          # indices of cytokines that should
          # be positive
          combn_ind_pos_mat <- combn_ind_pos_list[[n_pos_nm]]
          combn_nm_vec <- combn_nm_list[[n_pos_nm]]

          i <- 1
          for (i in seq_len(nrow(combn_ind_pos_mat))) {
            # if (debug) print(paste0("i: ", i))
            chnl_pos <- chnl[combn_ind_pos_mat[i, , drop = TRUE]]
            chnl_neg <- chnl[setdiff(seq_len(n_chnl), combn_ind_pos_mat[i, , drop = TRUE])]
            combn_nm <- combn_nm_vec[i]
            # now just add to the correct column,
            # then bind to unstim and to get the unstim counts
            combn_ind_pos <- outliergatev1:::.get_pos_ind_cyt_combn(
              ex = ex, gate_tbl = gate_tbl_ind,
              chnl_pos = chnl_pos, chnl_neg = chnl_neg,
              chnl_alt = NULL,
              gate_type_cyt_pos = "cyt",
              gate_type_single_pos = "single"
            )
            ann_tbl_faust_pop[, combn_nm] <- combn_ind_pos

            combn_ind_pos_uns <- outliergatev1:::.get_pos_ind_cyt_combn(
              ex = ex_uns, gate_tbl = gate_tbl_ind,
              chnl_pos = chnl_pos, chnl_neg = chnl_neg,
              chnl_alt = NULL,
              gate_type_cyt_pos = "cyt",
              gate_type_single_pos = "single"
            )
            ann_tbl_faust_pop_uns[, combn_nm] <- combn_ind_pos_uns
          }
        }

        ann_tbl_faust_pop <- ann_tbl_faust_pop |>
          dplyr::mutate(n_cell_pop = n_cell) |>
          dplyr::group_by(X1) |>
          dplyr::mutate(n_cell_pheno = dplyr::n()) |>
          dplyr::ungroup() |>
          tidyr::pivot_longer(
            -c(X1, n_cell_pop, n_cell_pheno),
            names_to = "combn",
            values_to = "pos_ind"
          ) |>
          dplyr::select(n_cell_pop, X1, n_cell_pheno, combn, pos_ind) |>
          dplyr::group_by(n_cell_pop, X1, n_cell_pheno, combn) |>
          dplyr::summarise(
            count = sum(pos_ind),
            .groups = "drop"
          )

        ann_tbl_faust_pop_uns <- ann_tbl_faust_pop_uns |>
          dplyr::mutate(n_cell_pop_uns = n_cell_uns) |>
          dplyr::group_by(X1) |>
          dplyr::mutate(n_cell_pheno_uns = dplyr::n()) |>
          dplyr::ungroup() |>
          tidyr::pivot_longer(
            -c(X1, n_cell_pop_uns, n_cell_pheno_uns),
            names_to = "combn",
            values_to = "pos_ind"
          ) |>
          dplyr::select(n_cell_pop_uns, X1, n_cell_pheno_uns, combn, pos_ind) |>
          dplyr::group_by(n_cell_pop_uns, X1, n_cell_pheno_uns, combn) |>
          dplyr::summarise(
            count_uns = sum(pos_ind),
            .groups = "drop"
          )

        ann_tbl <- ann_tbl_faust_pop |>
          dplyr::left_join(ann_tbl_faust_pop_uns,
            by = c("X1", "combn")
          )
        ann_tbl <- ann_tbl |>
          dplyr::rename(pheno = X1)

        pheno_vec_uni <- unique(ann_tbl$pheno)
        long_to_short_pheno <- purrr::map_chr(pheno_vec_uni, function(pheno) {
          for (i in seq_along(pop_defn)) {
            pttrn <- paste0(names(pop_defn)[i], "~", pop_defn[[i]], "~\\d~")
            for (pt in pttrn) {
              pheno <- gsub(pt, "", pheno)
            }
          }

          level_loc <- stringr::str_locate(pheno, "~\\d~\\d~")
          while (length(level_loc) > 0 & !is.na(level_loc[[1]])) {
            lvl <- ifelse(stringr::str_sub(pheno,
              start = level_loc[[1]] + 1,
              end = level_loc[[1]] + 1
            ) == "1",
            "-", "+"
            )
            str_start <- stringr::str_sub(pheno, end = level_loc[1] - 1)
            str_end <- stringr::str_sub(pheno, start = level_loc[2] + 1)
            pheno <- paste0(str_start, lvl, str_end)
            level_loc <- stringr::str_locate(pheno, "~\\d~\\d~")
          }

          pheno
        }) |>
          setNames(pheno_vec_uni)
        ann_tbl <- ann_tbl |>
          dplyr::mutate(
            pheno = long_to_short_pheno[pheno]
          )

        ann_tbl <- ann_tbl |>
          dplyr::group_by(pheno, combn) |>
          dplyr::summarise(
            n_cell_pop = n_cell_pop[1],
            n_cell_pop_uns = n_cell_pop_uns[1],
            n_cell_pheno = sum(n_cell_pheno),
            count = sum(count),
            n_cell_pheno_uns = sum(n_cell_pheno_uns),
            count_uns = sum(count_uns),
            .groups = "drop"
          )


        ann_tbl <- ann_tbl |>
          dplyr::mutate(
            combn = gsub("~", "", combn)
          )

        sample_info_tbl <- tibble::tibble(
          SubjectID = match_tbl_stim$SubjectID[1],
          VisitType = match_tbl_stim$VisitType[1],
          stim = match_tbl_stim$Stim[1],
          OrigFCSName = match_tbl_stim$OrigFCSName[1],
          MatchFCSName = match_tbl_stim$OrigFCSName[1],
          dir_faust = dir_sample,
          ind = ex$ind[1],
          pop = pop
        )
        saveRDS(
          object = ann_tbl,
          file = file.path(
            dir_save_pop_ind,
            "ann_tbl.rds"
          )
        )
        saveRDS(
          object = sample_info_tbl,
          file = file.path(
            dir_save_pop_ind,
            "sample_info_tbl.rds"
          )
        )
        saveRDS(
          object = TRUE,
          file = file.path(
            dir_save_pop_ind,
            "completed.rds"
          )
        )
      }
    }
  }

  out_tbl <- purrr::map_df(names(pop_to_defn), function(pop) {
    dir_save_pop <- file.path(dir_save, pop)
    dir_vec <- list.dirs(dir_save_pop, recursive = FALSE, full.names = FALSE)
    purrr::map(dir_vec, function(sample) {
      ann_tbl <- try(
        readRDS(
          file.path(dir_save_pop, sample, "ann_tbl.rds")
        ),
        silent = TRUE
      )
      sample_info_tbl <- try(
        readRDS(
          file.path(dir_save_pop, sample, "sample_info_tbl.rds")
        ),
        silent = TRUE
      )
      if (class(ann_tbl) == "try-error" ||
        class(sample_info_tbl)[1] == "try-error") {
        return(NULL)
      }
      sample_info_tbl |>
        dplyr::bind_cols(ann_tbl)
    }) |>
      purrr::compact() |>
      dplyr::bind_rows()
  })

  out_tbl
}
