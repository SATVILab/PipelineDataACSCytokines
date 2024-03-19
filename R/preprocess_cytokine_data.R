#' @title Preprocess cytokine data
#'
#' @param path_faust_dir character. Path to FAUST project directory.
#' @param label_cytokines function. A function that takes
#' the cytokine names as input and returns the
#' value desired to be plotted along the y-axis
#' of the grid plot.
#' If not \code{NULL}, then
#' it is supplied when creating the grid plot to
#' the \code{scale_y_continuous} function
#' via the \code{labels} parameter.
#' For example, if we have two cytokines,
#' \code{"IFNg"} and
#' \code{"TNF"}, but we want to
#' display \code{"IFNg"} with
#' the Greek gamma symbol, then we can
#' set \code{cyt_lab} equal to the following:
#' \code{cyt_lab = function(cyt)
#' purrr::map(cyt, function(cyt_ind) {switch(cyt_ind,
#' "IFNg" = bquote(paste(plain(paste("IFN")), gamma)),
#' cyt_ind)})}.
#' This will change the label for
#' \code{"IFNg"} but leave all the others as is.
#'
#' @export
#'
#' @importFrom flowCore exprs exprs<-
#' @importFrom rlang !!
#' @importFrom magrittr %>% %<>%
#' @import stringr
#' @import ggplot2
#' @import dplyr
preprocess_cytokine_data <-
  function(path_faust_dir,
           path_output,
           gs_nm = "gs_cytof_acs_nk",
           chnls_gate = c(
             "Ho165Di", "Gd158Di", "Nd146Di",
             "Dy164Di", "Gd156Di", "Nd150Di"
           ),
           chnls_analyse = "",
           gate = FALSE,
           fast = FALSE,
           stats_only = FALSE,
           plots_only = FALSE,
           save_cyt_pos_fcs = TRUE,
           combn_exc_list = list(
             NULL,
             c("Nd146Di", "Gd156Di", "Nd150Di")
           ),
           hladr = FALSE,
           hladr_bulk = FALSE,
           label_cytokines = PipelineDataACSCytokines::label_cytokines,
           compass = TRUE,
           compass_stim = c("p1", "mtb", "ebv", "p4"),
           compass_fast = FALSE,
           downstream_exc_list = NULL,
           tsne = TRUE,
           tsne_plot = TRUE,
           tsne_plot_force = TRUE,
           tsne_save = TRUE,
           tsne_reuse = TRUE,
           post_probs_bulk = TRUE,
           post_probs_bulk_stim = c("p1", "mtb", "ebv"),
           flowsom = TRUE,
           flowsom_stim = c("all_u"),
           flowsom_p4_exc = TRUE,
           flowsom_scale = FALSE,
           flowsom_n_clust = "n_f",
           flowsom_r_o = c("r_o", "all"),
           jacc_n = 1e2,
           jacc_reuse = TRUE,
           flowsom_post_probs = FALSE,
           flowsom_plot = TRUE,
           flowsom_plot_force = TRUE) {
    # ======================
    # Preparation
    # ======================

    # libraries
    on.exit(ggplot2::theme_set(theme_gray()))
    ggplot2::theme_set(cowplot::theme_cowplot())

    # data
    dir_gs <- normalizePath(file.path(
      path_faust_dir, "faustData",
      "gsData", gs_nm
    ))
    gs <- flowWorkspace::load_gs(dir_gs)
    if (fast) gs <- gs[1:10]

    # lists of markers for each cytokine
    # ------------------------

    # prep
    fdr_vec <- c(0.3)
    tol <- 1e-6

    gate_combn_list <- list("min" = c("loc"))
    gate_quant_list <- list(
      "min" = c(0.25, 0.5),
      "no" = c(0.2, 0.35),
      "prejoin" = c(0.15, 0.25)
    )
    marker_2d_vec <- c(
      "Eu151Di", "Ho165Di", "Gd158Di", "Nd146Di",
      "Lu175Di", "Gd156Di", "Dy164Di", "Nd150Di"
    )
    marker_2d_vec <- c(
      marker_2d_vec, "Sm154Di", "Nd145Di",
      "Er168Di", "Dy161Di", "Gd160Di", "Nd148Di",
      "Er166Di"
    )
    marker_2d_vec <- setNames(
      rep(1, length(marker_2d_vec)),
      marker_2d_vec
    )

    # actual marker list
    marker_list <- purrr::map(chnls_gate, function(chnl) {
      switch(chnl,
        "Ho165Di" = list(
          cut = "Ho165Di",
          pop_man_sub = NULL,
          high = marker_2d_vec,
          cps_scp = c(90, 95),
          pop_plot_sub = c(""),
          fdr = fdr_vec,
          tol = tol,
          gate_combn = gate_combn_list
        ),
        "Gd158Di" = list(
          cut = "Gd158Di",
          pop_man_sub = NULL,
          high = marker_2d_vec,
          cps_scp = c(90, 95),
          pop_plot_sub = c(""),
          fdr = fdr_vec,
          tol = tol,
          gate_combn = gate_combn_list
        ),
        "Nd146Di" = list(
          cut = "Nd146Di",
          pop_man_sub = NULL,
          high = marker_2d_vec,
          cps_scp = c(90, 95),
          pop_plot_sub = c(""),
          fdr = fdr_vec,
          tol = tol,
          gate_combn = gate_combn_list
        ),
        "Dy164Di" = list(
          cut = "Dy164Di",
          pop_man_sub = NULL,
          high = marker_2d_vec,
          cps_scp = c(90, 95),
          pop_plot_sub = c(""),
          fdr = fdr_vec,
          tol = tol,
          gate_combn = gate_combn_list
        ),
        "Gd156Di" = list(
          cut = "Gd156Di",
          pop_man_sub = NULL,
          high = marker_2d_vec,
          cps_scp = c(90, 95),
          pop_plot_sub = c(""),
          fdr = fdr_vec,
          tol = tol,
          gate_combn = gate_combn_list
        ),
        "Nd150Di" = list(
          cut = "Nd150Di",
          pop_man_sub = NULL,
          high = marker_2d_vec,
          cps_scp = c(90, 95),
          pop_plot_sub = c(""),
          fdr = fdr_vec,
          tol = tol,
          gate_combn = gate_combn_list
        ),
        stop(paste0("chnl ", chnl, " not recognised"))
      )
    })

    # get stim and fcs labels for each sample

    pop_gate <- "root"
    ind_in_batch_lab_vec <- c(
      "1" = "ebv", "2" = "mtb", "3" = "p1",
      "4" = "p4", "5" = "uns"
    )

    fcs_stim_tbl <- purrr::map_df(seq_along(gs), function(i) {
      fn <- flowCore::description(
        flowWorkspace::gh_pop_get_data(gs[[i]])
      )$FILENAME
      last_slash_loc_tbl <- str_locate_all(fn, "/")[[1]]
      last_slash_loc <- last_slash_loc_tbl[nrow(last_slash_loc_tbl), "end"][[1]]
      fn <- str_sub(fn, last_slash_loc + 1) %>%
        str_remove(".fcs")
      second_underscore_loc <- str_locate_all(fn, "_")[[1]][2, "end"][[1]]
      fcs <- str_sub(fn, end = second_underscore_loc - 1)
      stim <- str_split(fn, "_")[[1]]
      stim <- stim[[length(stim)]]
      tibble::tibble(fcs = fcs, stim = stim)
    })

    fcs_vec_gs <- fcs_stim_tbl$fcs
    fcs_lab_vec <- setNames(fcs_stim_tbl$fcs, seq_len(nrow(fcs_stim_tbl)))
    sampleid_lab_vec <- setNames(fcs_stim_tbl$fcs, seq_len(nrow(fcs_stim_tbl)))
    stim_lab_vec <- setNames(fcs_stim_tbl$stim, seq_len(nrow(fcs_stim_tbl)))

    # ========================
    # Identify gates
    # ========================

    if (gate || stats_only || plots_only) {
      message("gating")
      outliergatev1:::gate(
        data = gs,
        data_name = gs_nm,
        pop_gate = pop_gate,
        marker = marker_list,
        ind_in_batch_lab_vec = ind_in_batch_lab_vec,
        ind_in_batch_uns = 5,
        ind_in_batch_gate = 1:5,
        batch_size = 5,
        plot = FALSE,
        noise_sd = 2,
        bias_uns = 0.15, # original is 5/200
        min_bw = 0.3, # original is 7/200
        cp_min = 0,
        remove_beads = FALSE,
        boot_n = NULL,
        ind_skip = NULL,
        cytof_fcs_to_clin_map = NULL,
        fcs = fcs_vec_gs,
        max_pos_prob_x = 5,
        gate_quant = gate_quant_list,
        tol_gate = 1e-7,
        tol_gate_single = 1e-8,
        calc_cyt_pos_gates = TRUE,
        calc_single_pos_gates = TRUE,
        gate_name_plot = "locb0.15_min_clust",
        gate_name_stats = "locb0.15_min_clust",
        debug_stats = FALSE,
        sampleid_lab = sampleid_lab_vec,
        stim_lab = stim_lab_vec,
        stats_only = stats_only,
        plots_only = plots_only
      )
    }

    # ================================
    # Extract stats_combn_tbl
    # ================================

    gate_name_vec <- c("locb0.15_min_clust")
    chnl_lab <- UtilsCytoRSV::chnl_lab(
      flowWorkspace::gh_pop_get_data(gs[[1]])
    )
    ind_batch_list <- outliergatev1:::.get_ind_batch_list(
      data = gs, batch_size = 5, fcs = fcs_vec_gs
    )
    params_env <- list(
      pop_gate = "root",
      chnl_lab = chnl_lab,
      ind_in_batch_lab_vec = ind_in_batch_lab_vec,
      ind_in_batch_gate = 1:5, data_name = gs_nm,
      fcs = fcs_vec_gs,
      ind_in_batch_uns = 5,
      ind_batch_list = ind_batch_list,
      data = gs
    )
    path_base <- outliergatev1:::.create_dir_base(
      params = params_env %>%
        append(list(cut = marker_list[[1]]$cut)),
      dir_base_init = file.path(
        path_output,
        "datalarge"
      )
    ) %>%
      dirname()
    path_cyt_combn_base_stats <- file.path(
      path_base,
      paste0(purrr::map_chr(
        marker_list, function(x) x$cut
      ), collapse = "_")
    )

    path_stats <- file.path(
      path_cyt_combn_base_stats,
      "tables", "gate_stats.rds"
    )
    stats_combn_tbl <- readRDS(path_stats)
    saveRDS(
      object = paste0(chnls_gate),
      file = file.path(
        path_cyt_combn_base_stats,
        "chnls_gate.rds"
      )
    )

    path_cyt_combn_base <- file.path(
      path_base,
      paste0(chnls_analyse, collapse = "_")
    )
    if (!dir.exists(path_cyt_combn_base)) {
      dir.create(path_cyt_combn_base, recursive = TRUE)
    }
    readr::write_csv(
      x = data.frame(
        chnl = chnls_gate, marker = chnl_lab[chnls_gate]
      ),
      file = file.path(
        path_cyt_combn_base,
        "chnls_gate.csv"
      )
    )

    # ================================
    # Save cyt+ FCS files
    # ================================

    # directory for analysise

    gate_tbl <- readRDS(file.path(
      path_cyt_combn_base_stats, "tables", "gate_tbl.rds"
    ))

    if (save_cyt_pos_fcs) {
      message("saving cyt+ FCS files")
      purrr::walk(combn_exc_list, function(combn_exc) {
        purrr::walk(gate_name_vec, function(gn) {
          outliergatev1:::get_stim_pos_fcs(
            data = gs,
            data_name = gs_nm,
            dir_base = path_base,
            pop_gate = pop_gate,
            chnl = chnls_analyse,
            ind_in_batch_lab_vec = ind_in_batch_lab_vec,
            gate_name = gn,
            ind_in_batch_gate = 1:5,
            cytof_fcs_to_clin_map = NULL,
            ind_in_batch_uns = 5,
            gate_tbl = gate_tbl,
            gate_uns = TRUE,
            trans_fn = function(x) 5 * sinh(x),
            trans_chnl = c(
              "Dy161Di", "Dy162Di", "Dy163Di",
              "Dy164Di", "Er166Di", "Er167Di",
              "Er168Di", "Er170Di", "Eu151Di",
              "Eu153Di", "Gd155Di", "Gd156Di",
              "Gd158Di", "Gd160Di", "Ho165Di",
              "Lu175Di", "Lu176Di", "Nd142Di",
              "Nd143Di", "Nd144Di", "Nd145Di",
              "Nd146Di", "Nd148Di", "Nd150Di",
              "Pr141Di", "Sm147Di", "Sm149Di",
              "Sm152Di", "Sm154Di", "Tb159Di",
              "Tm169Di", "Yb171Di", "Yb172Di",
              "Yb173Di", "Yb174Di"
            ),
            combn_exc = combn_exc
          )
        })
      })
    }

    # ================================
    # Calculate HLADR expression on cyt+ cells
    # ================================

    if (hladr) {
      message("calcating HLADR expression on cyt+ cells")
      chnl_list <- list(c("Ho165Di", "Nd146Di"))
      mult_vec <- c(TRUE)
      hladr_tbl <- purrr::map_df(gate_name_vec, function(gate_name) {
        print(gate_name)
        purrr::map_df(chnl_list, function(chnl_vec) {
          print(chnl_vec)
          purrr::map_df(mult_vec, function(mult) {
            print(mult)
            outliergatev1:::get_hladr_med_diff(
              data = gs,
              data_name = gs_nm,
              pop_gate = pop_gate,
              chnl = chnl_vec,
              mult = mult,
              ind_in_batch_lab_vec = ind_in_batch_lab_vec,
              gate_name = gate_name,
              ind_in_batch_gate = 1:5,
              cytof_fcs_to_clin_map = NULL,
              ind_in_batch_uns = 5
            ) %>%
              dplyr::mutate(mult = mult)
          }) %>%
            dplyr::mutate(chnl = paste0(chnl_vec, collapse = "_"))
        }) %>%
          dplyr::mutate(gate_name = .env$gate_name)
      })
      analysispipeline::save_objects(
        "hladr" = hladr_tbl,
        dir_proj = path_cyt_combn_base,
        dir_sub = c("data_out", "hladr"),
        empty = TRUE
      )
    }

    # ================================
    # Calculate HLADR expression on bulk cells
    # ================================

    if (hladr_bulk) {
      message("calcating HLADR expression on bulk cells")
      hladr_tbl_bulk <- purrr::map_df(seq_along(gs), function(i) {
        fr <- flowWorkspace::gh_pop_get_data(gs[[i]])
        ex <- exprs(fr) %>% tibble::as_tibble()
        tibble::tibble(SampleID = fcs_lab_vec[i], stim = stim_lab_vec[i]) %>%
          dplyr::mutate(
            SubjectID = str_sub(SampleID, end = 6),
            VisitType = str_sub(SampleID, start = 8)
          ) %>%
          dplyr::select(SubjectID, VisitType, stim) %>%
          dplyr::mutate(
            med = median(ex[["Eu151Di"]]),
            mean = mean(ex[["Eu151Di"]]),
            q75 = quantile(ex[["Eu151Di"]], 0.75),
            q90 = quantile(ex[["Eu151Di"]], 0.9),
            q25 = quantile(ex[["Eu151Di"]], 0.25)
          )
      })
      analysispipeline::save_objects(
        "hladr_bulk" = hladr_tbl_bulk,
        dir_proj = path_cyt_combn_base,
        dir_sub = c("data_out", "hladr_bulk"),
        empty = TRUE
      )
    }

    # ================================
    # Format stats_combn_tbl
    # ================================

    chnl_vec_all <- stringr::str_split(
      stats_combn_tbl$cyt_combn[1],
      "~[+-]~"
    )[[1]]
    chnl_vec_all <- chnl_vec_all[chnl_vec_all != ""]
    chnl_vec_exc <- setdiff(chnl_vec_all, chnls_analyse)
    for (i in seq_along(chnl_vec_exc)) {
      stats_combn_tbl <- stats_combn_tbl %>%
        dplyr::mutate(
          cyt_combn = stringr::str_remove(
            cyt_combn,
            paste0(chnl_vec_exc[i], "~[+-]~")
          )
        )
    }
    if (length(chnl_vec_exc) > 0) {
      stats_combn_tbl <- stats_combn_tbl %>%
        dplyr::group_by(
          gate_name, ind, SampleID, stim, cyt_combn
        ) %>%
        dplyr::summarise(
          count_stim = sum(count_stim),
          n_cell_stim = n_cell_stim[1],
          count_uns = sum(count_uns),
          n_cell_uns = n_cell_uns[1],
          prop_stim = sum(prop_stim),
          prop_uns = sum(prop_uns),
          prop_bs = sum(prop_bs),
          freq_stim = sum(freq_stim),
          freq_uns = sum(freq_uns),
          freq_bs = sum(freq_bs),
          .groups = "drop"
        )
    }

    chnl_lab_vec <- c(
      "Ho165Di" = "IFNg", "Gd158Di" = "IL2", "Nd146Di" = "TNF",
      "Dy164Di" = "IL17", "Nd150Di" = "IL22", "Gd156Di" = "IL6"
    )

    cyt_combn_rep_lab_vec <- purrr::map_chr(
      unique(stats_combn_tbl$cyt_combn),
      function(cco) {
        chnl_vec <- str_split(cco, "~[-+]~")[[1]]
        chnl_vec <- chnl_vec[!chnl_vec == ""]
        pos_ind_vec <- rep(NA, length(chnl_vec))
        cco_rem <- cco
        for (i in seq_along(chnl_vec)) {
          cco_rem %<>% str_remove(chnl_vec[i])
          pos_ind_vec[i] <- str_sub(cco_rem, 2 + (i - 1) * 3, 2 + (i - 1) * 3)
        }
        pos_ind_lab_vec <- setNames(c("&", "&!"), c("+", "-"))
        cyt_combn_rep_end <- paste0(
          pos_ind_lab_vec[pos_ind_vec[-1]],
          chnl_lab_vec[chnl_vec[-1]],
          sep = "", collapse = ""
        )
        cyt_combn_rep_start <- paste0(
          setNames(c("", "!"), c("+", "-"))[pos_ind_vec[1]],
          chnl_lab_vec[chnl_vec[1]]
        )
        paste0(cyt_combn_rep_start, cyt_combn_rep_end)
      }
    ) %>%
      setNames(unique(stats_combn_tbl$cyt_combn))

    stats_combn_tbl <- stats_combn_tbl %>%
      dplyr::mutate(
        stim = stim_lab_vec[as.character(ind)],
        batch_sh = fcs_lab_vec[as.character(ind)]
      ) %>%
      dplyr::mutate(
        SubjectID = purrr::map_chr(
          batch_sh,
          function(x) str_split(x, "_")[[1]][[1]]
        ),
        VisitType = purrr::map_chr(
          batch_sh,
          function(x) str_split(x, "_")[[1]][[2]]
        )
      ) %>%
      dplyr::select(
        gate_name, ind, batch_sh, SubjectID, VisitType,
        stim, cyt_combn, count_stim, n_cell_stim,
        count_uns, n_cell_uns
      )

    stats_combn_tbl <- stats_combn_tbl %>% dplyr::mutate(
      cyt_combn = cyt_combn_rep_lab_vec[cyt_combn]
    )

    all_neg_cyt_combn <- str_remove_all(cyt_combn_rep_lab_vec[[1]], "[!]") %>%
      str_replace_all("[&]", "&!")
    all_neg_cyt_combn <- paste0("!", all_neg_cyt_combn)

    stats_combn_tbl <- stats_combn_tbl %>%
      dplyr::filter(!cyt_combn == all_neg_cyt_combn)

    stats_combn_tbl_neg <- stats_combn_tbl %>%
      dplyr::group_by(ind, gate_name, batch_sh, SubjectID, VisitType, stim) %>%
      dplyr::summarise(
        n_cell_stim = n_cell_stim[1],
        n_cell_uns = n_cell_uns[1],
        count_stim = n_cell_stim[1] - sum(count_stim),
        count_uns = n_cell_uns - sum(count_uns)
      ) %>%
      dplyr::mutate(cyt_combn = all_neg_cyt_combn) %>%
      dplyr::ungroup()

    stats_combn_tbl <- stats_combn_tbl %>% dplyr::bind_rows(stats_combn_tbl_neg)

    stats_combn_tbl <- stats_combn_tbl %>%
      dplyr::mutate(stim = ifelse(stim == "mtbaux", "mtb", stim))

    # ========================
    # Save results thus far and extract old required results
    # ========================

    # create name of datset to save final results in package under
    # --------------------

    ds_name <- stringr::str_remove(gs_nm, "gs_cytof_acs_")
    if (all(c("Ho165Di", "Gd158Di", "Nd146Di") %in%
      chnls_analyse)) {
      ds_name <- paste0(ds_name, "_th1")
    } else {
      if ("Ho165Di" %in% chnls_analyse) ds_name <- ds_name %>% paste0("_ifng")
      if ("Gd158Di" %in% chnls_analyse) ds_name <- ds_name %>% paste0("_il2")
      if ("Nd146Di" %in% chnls_analyse) ds_name <- ds_name %>% paste0("_tnf")
    }
    if ("Dy164Di" %in% chnls_analyse) ds_name <- ds_name %>% paste0("_il17")
    if ("Nd150Di" %in% chnls_analyse) ds_name <- ds_name %>% paste0("_il22")
    if ("Gd156Di" %in% chnls_analyse) ds_name <- ds_name %>% paste0("_il6")

    # save results
    # ---------------------

    out_list <- list("stats_combn_tbl" = stats_combn_tbl)

    if (!compass) {
      path_compass_old <- file.path(
        path_cyt_combn_base, "data_out", "compass", "compass.rds"
      )
      if (file.exists(path_compass_old)) {
        out_list$compass <- readRDS(path_compass_old)
      }
    }

    if (!hladr) {
      path_hladr_old <- file.path(
        path_cyt_combn_base, "data_out", "hladr", "hladr.rds"
      )
      if (file.exists(path_hladr_old)) {
        out_list$hladr <- readRDS(path_hladr_old)
      }
    } else {
      out_list$hladr <- hladr_tbl
    }


    if (!hladr_bulk) {
      path_hladr_old <- file.path(
        path_cyt_combn_base, "data_out", "hladr_bulk", "hladr_bulk.rds"
      )
      if (file.exists(path_hladr_old)) {
        out_list$hladr_bulk <- readRDS(path_hladr_old)
      }
    } else {
      out_list$hladr_bulk <- hladr_tbl_bulk
    }

    if (!flowsom) {
      path_flowsom_old <- file.path(
        path_cyt_combn_base, "data_out", "flowsom", "flowsom.rds"
      )
      if (file.exists(path_flowsom_old)) {
        out_list$flowsom <- readRDS(path_flowsom_old)
      }
    }

    # ==========================
    # COMPASS
    # ==========================

    if (compass) {
      message("running COMPASS")
      compass_list <- purrr::map(gate_name_vec, function(gate_name_curr) {
        print(gate_name_curr)
        purrr::map(compass_stim, function(stim_curr) {
          print(stim_curr)
          stats_combn_tbl_curr <- stats_combn_tbl %>%
            dplyr::filter(
              gate_name == gate_name_curr,
              stim == stim_curr
            )

          # =============================
          # Get the SampleID and other metatadata
          # =============================

          full_tbl <- stats_combn_tbl_curr %>%
            dplyr::select(
              SubjectID, VisitType, stim,
              gate_name, cyt_combn, n_cell_stim,
              count_stim, n_cell_uns, count_uns
            ) %>%
            dplyr::left_join(
              DataTidyACSClinical::clinical_data %>%
                dplyr::select(
                  SubjectID, Progressor, timeToTB,
                  timeToTBFromVisit, SampleID, VisitType
                ),
              by = c("SubjectID", "VisitType")
            )

          # =============================
          # Get data_uns
          # =============================

          data_uns <- full_tbl %>%
            dplyr::select(
              SampleID, cyt_combn,
              count_uns, n_cell_uns
            ) %>%
            dplyr::rename(
              count = count_uns,
              n_cell = n_cell_uns
            )

          data_uns <- data_uns %>%
            dplyr::select(-n_cell) %>%
            tidyr::pivot_wider(
              id_cols = "SampleID",
              names_from = cyt_combn,
              values_from = count
            )

          id_vec_uns <- data_uns$SampleID

          data_uns <- data_uns %>%
            dplyr::select(-SampleID) %>%
            dplyr::mutate_all(as.numeric)

          data_uns <- data_uns %>% as.matrix()
          rownames(data_uns) <- id_vec_uns

          # =============================
          # Get data_stim
          # =============================

          data_stim <- full_tbl %>%
            dplyr::select(SampleID, cyt_combn, count_stim, n_cell_stim) %>%
            dplyr::rename(count = count_stim, n_cell = n_cell_stim)

          data_stim <- data_stim %>%
            dplyr::select(-n_cell) %>%
            tidyr::pivot_wider(
              id_cols = "SampleID",
              names_from = cyt_combn,
              values_from = count
            )

          id_vec_stim <- data_stim$SampleID

          data_stim <- data_stim %>%
            dplyr::select(-SampleID) %>%
            dplyr::mutate_all(as.numeric)

          data_stim <- data_stim %>% as.matrix()
          rownames(data_stim) <- id_vec_stim

          # =============================
          # Get data_meta
          # =============================


          data_stim <- data_stim %>% as.matrix()
          rownames(data_stim) <- id_vec_stim

          if (!identical(rownames(data_stim), rownames(data_uns))) {
            stop("stim and uns not identically ordered")
          }

          data_meta <- full_tbl %>%
            dplyr::select(
              SubjectID, VisitType, SampleID,
              gate_name, stim, timeToTB, Progressor
            ) %>%
            dplyr::group_by(
              SampleID, SubjectID, VisitType,
              gate_name, stim, timeToTB, Progressor
            ) %>%
            dplyr::slice(1) %>%
            dplyr::ungroup() %>%
            as.data.frame()

          data_meta <- purrr::map_df(rownames(data_stim), function(id) {
            data_meta %>%
              dplyr::filter(.data$SampleID == id)
          })

          if (!identical(data_meta$SampleID, rownames(data_uns))) {
            stop("meta and uns not identically ordered")
          }

          # =============================
          # Run COMPASS
          # =============================

          fit <- COMPASS::SimpleCOMPASS(
            n_s = data_stim,
            n_u = data_uns,
            meta = data_meta,
            individual_id = "SampleID",
            iterations = ifelse(
              compass_fast, 1e2, 4e4
            ),
            replications = ifelse(
              compass_fast, 2, 8
            )
          )

          # data to NULL
          fit$data$n_s <- NULL
          fit$data$n_u <- NULL
          fit$data$categories <- NULL

          # fit elements to NULL
          fit$fit$alpha_s <- NULL
          fit$fit$alpha_u <- NULL
          fit$fit$gamma <- NULL
          fit$fit$A_gamma <- NULL
          fit$fit$A_alphas <- NULL
          fit$fit$A_alphau <- NULL
          # fit$fit$categories <- NULL
          fit$fit$model <- NULL
          fit$data$counts_s <- NULL
          fit$data$counts_u <- NULL

          fit
        }) %>%
          setNames(compass_stim)
      }) %>%
        setNames(gate_name_vec)

      analysispipeline::save_objects(
        "compass" = compass_list,
        dir_proj = path_cyt_combn_base,
        dir_sub = c("data_out", "compass"),
        empty = TRUE
      )

      out_list <- out_list %>% append(list("compass" = compass_list))
    }

    # ==========================
    # Expression extraction and t-SNE
    # ==========================


    # preparation
    # ----------------------

    # channel label vectors
    df_lab <- flowCore::parameters(flowWorkspace::gh_pop_get_data(gs[[1]]))@data
    chnl_lab_vec_fcs <- setNames(df_lab$desc, df_lab$name)
    chnl_lab_vec_fcs <- chnl_lab_vec_fcs[!is.na(chnl_lab_vec_fcs)]

    # get exclusion options for downstream t-SNE,
    # post probs and flowSOM (as requested)
    if (is.null(downstream_exc_list)) {
      exc_vec <- list.dirs(
        path_cyt_combn_base,
        full.names = FALSE,
        recursive = FALSE
      )
      exc_vec <- exc_vec[stringr::str_detect(exc_vec, "exc-")]
    } else { # if specified, use these channels
      exc_vec <- purrr::map_chr(downstream_exc_list, function(downstream_exc) {
        if (is.null(downstream_exc)) {
          return("exc-none")
        }
        paste0("exc-", paste0(downstream_exc, collapse = "_&_"))
      })
    }

    # stop if there is nothing further to work with
    if (length(exc_vec) == 0) {
      # save_final
      # --------------
      assign(ds_name, out_list)
      eval(parse(text = paste0(
        "usethis::use_data(", ds_name,
        ", overwite = TRUE)"
      )))
      return(invisible(TRUE))
    }

    # run expression extraction and possibly t-SNE
    # ----------------------------

    if (tsne || flowsom) {
      message("calculating tsne")
      ex_ag_list <- purrr::map(exc_vec, function(exc) {
        print(exc)

        dir_fcs_exc <- file.path(path_cyt_combn_base, exc)
        if (!dir.exists(dir_fcs_exc)) {
          warning(paste0(
            dir_fcs_exc,
            " does not exist for getting cyt+ expression data"
          ))
          return(NULL)
        }
        gn_vec <- list.dirs(dir_fcs_exc, full.names = FALSE, recursive = FALSE)
        if (length(gn_vec) == 0) {
          if (!dir.exists(dir_fcs_exc)) {
            warning(paste0(
              "No sub-directories in ",
              dir_fcs_exc,
              " for getting cyt+ expression data"
            ))
            return(NULL)
          }
        }
        purrr::map(gn_vec, function(gn) {
          print(gn)
          dir_fcs_exc_gn_fcs <- file.path(dir_fcs_exc, gn, "fcs")
          if (!dir.exists(dir_fcs_exc_gn_fcs)) {
            return(NULL)
          }

          # prepare data
          # -------------------------

          # load
          fcs_vec <- list.files(dir_fcs_exc_gn_fcs)
          if (length(fcs_vec) == 0) {
            return(NULL)
          }
          ex_tbl <- purrr::map(fcs_vec, function(x) {
            ex <- try(flowCore::exprs(
              flowCore::read.FCS(file.path(dir_fcs_exc_gn_fcs, x))
            ) %>%
              tibble::as_tibble() %>%
              dplyr::mutate(fcs = x))
            if (identical(class(ex), "try-error")) {
              return(NULL)
            }
            ex
          }) %>%
            purrr::compact() %>%
            dplyr::bind_rows()

          # format
          chnl_lab_vec_fcs <- c(chnl_lab_vec_fcs, c("fcs" = "fcs"))
          marker_lab_vec_fcs <- c(names(chnl_lab_vec_fcs), chnl_lab_vec_fcs)

          # remove columnns without descriptors and CD45
          ex_tbl <- ex_tbl[, colnames(ex_tbl) %in% names(chnl_lab_vec_fcs)]
          ex_tbl <- ex_tbl[
            ,
            !colnames(ex_tbl) %in%
              (names(chnl_lab_vec_fcs)[chnl_lab_vec_fcs == "CD45"])
          ]


          # rename columns
          colnames(ex_tbl) <- chnl_lab_vec_fcs[colnames(ex_tbl)]

          # transform columns
          ex_tbl <- ex_tbl %>%
            dplyr::select(-fcs) %>%
            dplyr::mutate_all(function(x) asinh(x / 5)) %>%
            dplyr::mutate(fcs = ex_tbl$fcs)

          # get sample info from FCS name

          fcs_vec <- ex_tbl$fcs

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

          ex_tbl <- ex_tbl %>%
            dplyr::mutate(fcs = fcs_vec) %>%
            dplyr::left_join(sample_info_tbl)

          ex_tbl
        }) %>%
          setNames(gn_vec) %>%
          purrr::compact()
      }) %>%
        setNames(exc_vec) %>%
        purrr::compact()

      # chnls to include in tsne
      # ------------------
      df_lab <- flowCore::parameters(flowWorkspace::gh_pop_get_data(gs[[1]]))@data
      chnl_vec_sel <- c(
        "Er166Di", "Ho165Di", "Gd155Di", "Dy163Di", "Dy164Di",
        "Lu175Di", "Tm169Di", "Yb174Di", "Sm149Di", "Eu153Di",
        "Nd143Di", "Nd144Di", "Lu176Di", "Yb172Di", "Pr141Di",
        "Yb171Di", "Er167Di", "Nd148Di", "Dy162Di", "Nd150Di",
        "Gd158Di", "Nd146Di", "Eu151Di", "Gd156Di", "Tb159Di",
        "Yb173Di", "Er170Di", "Nd142Di"
      )
      chnl_lab_vec_fcs <- setNames(df_lab$desc %>%
        str_remove("-beads") %>%
        str_remove("-Beads"), df_lab$name)
      chnl_lab_vec_fcs <- chnl_lab_vec_fcs[!is.na(chnl_lab_vec_fcs)]
      chnl_lab_vec_fcs["Nd146Di"] <- "TNF"
      chnl_lab_vec_fcs["Lu175Di"] <- "Perforin"
      chnl_lab_vec_fcs_no_adj <- setNames(df_lab$desc, df_lab$name)
      chnl_lab_vec_fcs <- chnl_lab_vec_fcs[!is.na(chnl_lab_vec_fcs)]
      marker_vec_sel <- chnl_lab_vec_fcs[chnl_vec_sel]


      tsne_list <- purrr::map(names(ex_ag_list), function(exc) {
        purrr::map(names(ex_ag_list[[exc]]), function(gn) {
          if (!tsne_reuse) {
            return(NULL)
          }
          path_tsne_old <- file.path(
            path_cyt_combn_base, exc, gn,
            "data_int", "tsne", "tsne.rds"
          )
          if (!file.exists(path_tsne_old)) {
            return(NULL)
          }
          readRDS(path_tsne_old)
        }) %>%
          setNames(names(ex_ag_list[[exc]]))
      }) %>%
        setNames(names(ex_ag_list))

      # run t-SNE plots
      # ------------------

      tsne_list <- purrr::map(exc_vec, function(exc) {
        gn_vec <- names(ex_ag_list[[exc]])
        gn <- gn_vec[1]
        purrr::map(gn_vec, function(gn) {
          # whether to rerun t-SNE or not
          # --------------------

          if (tsne_reuse) {
            tsne_old <- tsne_list[[exc]][[gn]]
            tsne_rerun <- is.null(tsne_old)
          } else {
            tsne_rerun <- TRUE
          }

          # get tsne data
          # ---------------------

          if (tsne_rerun) {
            ex_tbl <- ex_ag_list[[exc]][[gn]]
            for (cn in colnames(ex_tbl)) {
              if (!cn %in% marker_vec_sel) {
                ex_tbl <- ex_tbl[, -which(colnames(ex_tbl) == cn)]
              }
            }
            ex_mat <- ex_tbl %>% as.matrix()
            rownames(ex_mat) <- ex_ag_list[[exc]][[gn]]$fcs

            # apply umap
            time_umap_start <- proc.time()[3]
            set.seed(21)
            umap_mat <- uwot::tumap(ex_mat,
              n_threads = 12
            )
            time_umap_total <- proc.time()[3] - time_umap_start
            time_umap_total

            # apply t-SNE
            time_tsne_start <- proc.time()[3]
            tsne_mat <- fftRtsne(ex_mat,
              rand_seed = 210692,
              perplexity = 35,
              max_iter = 3e3
            )
            time_tsne_total <- proc.time()[3] - time_tsne_start
            time_tsne_total

            tsne_tbl <- ex_ag_list[[exc]][[gn]] %>%
              dplyr::mutate(
                tsne1 = tsne_mat[, 1],
                tsne2 = tsne_mat[, 2],
                umap1 = umap_mat[, 1],
                umap2 = umap_mat[, 2]
              )

            analysispipeline::save_objects(
              "tsne" = tsne_tbl,
              dir_proj = path_cyt_combn_base,
              dir_sub = c(exc, gn, "data_int", "tsne"),
              empty = TRUE
            )
          } else {
            tsne_tbl <- tsne_list[[exc]][[gn]]
          }

          # plot tsne data
          # -------------------

          # Whether to plot or not. can force it,
          # or can happen if t-SNE recalculated
          tsne_plot_force_ctb <- tsne_plot_force
          tsne_plot_ctb <- tsne_plot && tsne_rerun
          tsne_plot_final <- tsne_plot_force_ctb || tsne_plot_ctb

          .plot_tsne_bulk <- function(dir_base, exc, gn) {

          }

          if (tsne_plot_final) {
            plot_tbl <- tsne_tbl %>%
              dplyr::select(
                CD33:`IFNg-beads`,
                `Perf-beads`:CCR6,
                CD20:CD38,
                fcs:umap2
              ) %>%
              tidyr::pivot_longer(CD33:CD38,
                names_to = "marker",
                values_to = "expr"
              )


            if (!tsne_plot) {
              return(tsne_tbl)
            }

            p_marker_tsne <- ggplot(
              plot_tbl %>%
                dplyr::mutate(
                  marker = str_remove(marker, "-beads") %>%
                    str_remove("-Beads")
                ) %>%
                dplyr::mutate(expr = log(pmax(expr, 0.5))),
              aes(x = tsne1, y = tsne2, col = expr)
            ) +
              coord_equal() +
              geom_point(size = 0.25) +
              scale_colour_viridis_c() +
              facet_wrap(~marker, ncol = 5) +
              theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.title = element_blank()
              )

            analysispipeline::save_objects(
              p_marker_tsne = p_marker_tsne,
              dir_proj = path_cyt_combn_base,
              dir_sub = c(exc, gn, "tsne"),
              height = 45,
              width = 50,
              empty = FALSE
            )

            p_marker_umap <- ggplot(
              plot_tbl %>%
                dplyr::mutate(marker = str_remove(marker, "-beads") %>%
                  str_remove("-Beads")) %>%
                dplyr::mutate(expr = log(pmax(expr, 0.5))),
              aes(x = umap1, y = umap2, col = expr)
            ) +
              coord_equal() +
              geom_point(size = 0.2) +
              scale_colour_viridis_c() +
              facet_wrap(~marker, ncol = 5) +
              theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.title = element_blank()
              )

            analysispipeline::save_objects(
              p_marker_umap = p_marker_umap,
              dir_proj = path_cyt_combn_base,
              dir_sub = c(exc, gn, "tsne"),
              height = 45,
              width = 50,
              empty = FALSE
            )

            p_prog_stim_equal_cells <- ggplot(
              tsne_tbl %>%
                dplyr::select(fcs:tsne2),
              aes(x = tsne1, y = tsne2)
            ) +
              coord_equal() +
              geom_hex(
                bins = 32,
                data = tsne_tbl %>%
                  dplyr::select(fcs:tsne2) %>%
                  dplyr::group_by(Progressor, stim) %>%
                  dplyr::mutate(n_cell = n()) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(min_n_cell = min(n_cell)) %>%
                  dplyr::group_by(Progressor, stim) %>%
                  dplyr::sample_n(min_n_cell[1]) %>%
                  dplyr::ungroup()
              ) +
              geom_density_2d(
                col = "orange",
                size = 0.5
              ) +
              scale_fill_viridis_c(trans = "log") +
              facet_wrap(Progressor ~ stim, ncol = 4) +
              theme(
                axis.text = element_blank(),
                axis.ticks = element_blank()
              ) +
              labs(x = "t-SNE 1", y = "t-SNE 2")

            p_prog_stim_raw <- ggplot(
              tsne_tbl %>%
                dplyr::select(fcs:tsne2),
              aes(x = tsne1, y = tsne2)
            ) +
              coord_equal() +
              geom_hex(bins = 32) +
              geom_density_2d(
                col = "orange",
                size = 0.5
              ) +
              scale_fill_viridis_c(trans = "log") +
              facet_wrap(Progressor ~ stim, ncol = 4) +
              theme(
                axis.text = element_blank(),
                axis.ticks = element_blank()
              ) +
              labs(x = "t-SNE 1", y = "t-SNE 2")

            p_prog_stim_equal_cells_hex_only <- ggplot(
              tsne_tbl %>%
                dplyr::select(fcs:tsne2),
              aes(x = tsne1, y = tsne2)
            ) +
              coord_equal() +
              geom_hex(
                bins = 32,
                data = tsne_tbl %>%
                  dplyr::select(fcs:tsne2) %>%
                  dplyr::group_by(Progressor, stim) %>%
                  dplyr::mutate(n_cell = n()) %>%
                  dplyr::ungroup() %>%
                  dplyr::mutate(min_n_cell = min(n_cell)) %>%
                  dplyr::group_by(Progressor, stim) %>%
                  dplyr::sample_n(min_n_cell[1]) %>%
                  dplyr::ungroup()
              ) +
              scale_fill_viridis_c(trans = "log") +
              facet_wrap(Progressor ~ stim, ncol = 4) +
              theme(
                axis.text = element_blank(),
                axis.ticks = element_blank()
              ) +
              labs(x = "t-SNE 1", y = "t-SNE 2")

            p_prog_stim_raw_hex_only <- ggplot(
              tsne_tbl %>%
                dplyr::select(fcs:tsne2),
              aes(x = tsne1, y = tsne2)
            ) +
              coord_equal() +
              geom_hex(bins = 32) +
              scale_fill_viridis_c(trans = "log") +
              facet_wrap(Progressor ~ stim, ncol = 4) +
              theme(
                axis.text = element_blank(),
                axis.ticks = element_blank()
              ) +
              labs(x = "t-SNE 1", y = "t-SNE 2")

            p_prog_stim_dens_only <- ggplot(
              tsne_tbl %>%
                dplyr::select(fcs:tsne2),
              aes(
                x = tsne1, y = tsne2,
                col = Progressor
              )
            ) +
              coord_equal() +
              geom_density_2d(size = 0.5) +
              scale_colour_manual(values = c(
                "Progressor" = "orange2",
                "Control" = "dodgerblue"
              )) +
              facet_wrap(Progressor ~ stim, ncol = 4) +
              theme(
                axis.text = element_blank(),
                axis.ticks = element_blank()
              ) +
              labs(x = "t-SNE 1", y = "t-SNE 2")

            analysispipeline::save_objects(
              p_prog_stim_dens_only = p_prog_stim_dens_only,
              p_prog_stim_raw_hex_only = p_prog_stim_raw_hex_only,
              p_prog_stim_raw = p_prog_stim_raw,
              p_prog_stim_equal_cells_hex_only =
                p_prog_stim_equal_cells_hex_only,
              p_prog_stim_equal_cells = p_prog_stim_equal_cells,
              dir_proj = path_cyt_combn_base,
              dir_sub = c(exc, gn, "tsne"),
              height = 30,
              width = 35,
              empty = FALSE
            )

            analysispipeline::save_objects(
              p_prog_stim_dens_only = p_prog_stim_dens_only,
              p_prog_stim_raw_hex_only =
                p_prog_stim_raw_hex_only,
              p_prog_stim_raw = p_prog_stim_raw,
              p_prog_stim_equal_cells_hex_only =
                p_prog_stim_equal_cells_hex_only,
              p_prog_stim_equal_cells = p_prog_stim_equal_cells,
              dir_proj = path_cyt_combn_base,
              dir_sub = c(exc, gn, "tsne"),
              height = 30,
              width = 35,
              empty = FALSE
            )

            p_list_save <- list(
              p_prog_stim_dens_only =
                p_prog_stim_dens_only,
              p_prog_stim_raw_hex_only =
                p_prog_stim_raw_hex_only,
              p_prog_stim_raw =
                p_prog_stim_raw,
              p_prog_stim_equal_cells_hex_only =
                p_prog_stim_equal_cells_hex_only,
              p_prog_stim_equal_cells =
                p_prog_stim_equal_cells
            )
            p_list_save <- purrr::map(p_list_save, function(p) {

            })
          }

          # return output
          tsne_tbl
        }) %>%
          setNames(names(ex_ag_list[[exc]]))
      }) %>%
        setNames(names(ex_ag_list))

      if (tsne_save) {
        analysispipeline::save_objects(
          "tsne" = tsne_list,
          dir_proj = path_cyt_combn_base,
          dir_sub = c("data_int", "tsne"),
          empty = TRUE
        )
      }
    }

    # ==========================
    # Bulk probability of responding per individual
    # ==========================

    if (flowsom || post_probs_bulk) {
      message("calculating responder probabilities")

      exc <- exc_vec[1]
      responders_list <- purrr::map(exc_vec, function(exc) {
        # create cyt_combn from exc
        if (exc != "exc-none") {
          exc_wk <- exc %>% stringr::str_remove("exc-")
          chnl_combn_exc_vec_raw <- stringr::str_split(exc_wk, "&")[[1]]
          cyt_vec <- stats_combn_tbl$cyt_combn[1] %>%
            stringr::str_split("&") %>%
            magrittr::extract2(1)
          cyt_vec <- purrr::map_chr(cyt_vec, function(x) {
            if (!stringr::str_sub(x, end = 1) == "!") {
              return(x)
            }
            stringr::str_sub(x, start = 2)
          })
          marker_lab_vec_cyt <- setNames(
            names(chnl_lab_vec), chnl_lab_vec
          )
          chnl_vec <- marker_lab_vec_cyt[cyt_vec]

          cyt_combn_exc_vec <- purrr::map_chr(
            chnl_combn_exc_vec_raw,
            function(x) {
              chnl_vec_exc <- stringr::str_split(x, "_")[[1]]
              chnl_vec_exc <- chnl_vec_exc[!chnl_vec_exc == ""]
              if (!all(chnl_vec_exc %in% chnl_vec)) {
                stop("not all channels to exclude in modelled channels")
              }
              purrr::map_chr(chnl_vec, function(chnl) {
                if (chnl %in% chnl_vec_exc) {
                  return(
                    paste0("!", chnl_lab_vec[chnl])
                  )
                }
                chnl_lab_vec[chnl]
              }) %>%
                paste0(collapse = "&")
            }
          )
        } else {
          cyt_combn_exc_vec <- NULL
        }


        dir_fcs_exc <- file.path(path_cyt_combn_base, exc)
        gn_vec <- list.dirs(dir_fcs_exc, full.names = FALSE, recursive = FALSE)
        if (length(gn_vec) == 0) {
          return(NULL)
        }
        gn <- gn_vec[1]
        purrr::map(gn_vec, function(gn) {
          compass_list <- out_list$compass[[gn]]
          purrr::map(names(compass_list), function(stim) {
            UtilsCompassSV::response_prob(
              c_obj = compass_list[[stim]],
              exc = cyt_combn_exc_vec
            ) %>%
              dplyr::mutate(stim = .env$stim)
          }) %>%
            setNames(names(compass_list))
        }) %>%
          setNames(gn_vec)
      }) %>%
        setNames(exc_vec)

      out_list <- out_list %>%
        append(list("post_probs_bulk" = responders_list))
    }

    # ===========================
    # FlowSOM
    # ===========================

    if (flowsom) {
      message("running FlowSOM")

      stim_vec_u <- c("p1", "mtb", "ebv", "p4", "uns")
      if (flowsom_p4_exc) stim_vec_u <- stim_vec_u[!stim_vec_u == "p4"]
      marker_lab_vec_fcs <- setNames(names(chnl_lab_vec_fcs), chnl_lab_vec_fcs)

      responders_only <- flowsom_r_o[1]
      for (responders_only in flowsom_r_o) {
        print(responders_only)

        exc <- exc_vec[1]
        for (exc in exc_vec) {
          print(exc)

          dir_fcs_exc <- file.path(path_cyt_combn_base, exc)
          gn_vec <- list.dirs(dir_fcs_exc, full.names = FALSE, recursive = FALSE)
          if (length(gn_vec) == 0) next
          gn <- gn_vec[1]

          for (gn in gn_vec) {
            print(gn)
            dir_fcs_exc_gn_fcs <- file.path(dir_fcs_exc, gn, "fcs")
            if (!dir.exists(dir_fcs_exc_gn_fcs)) next

            # load
            fcs_vec <- list.files(dir_fcs_exc_gn_fcs)
            if (length(fcs_vec) == 0) next
            dir_fcs_exc_gn_fcs <- file.path(path_cyt_combn_base, exc, gn, "fcs")

            stim <- flowsom_stim[1]
            for (stim in flowsom_stim) {
              print(stim)

              # -----------------
              # Prepare data
              # -----------------

              # get list of responders
              # -----------------

              responders_list_resp_stim <- switch(as.character(responders_only == "r_o"),
                "FALSE" = NULL,
                "TRUE" = purrr::map(setdiff(stim_vec_u, "uns"), function(stim) {
                  responders_list[[exc]][[gn]][[stim]] %>%
                    dplyr::filter(prob > 0.75) %>%
                    magrittr::extract2("sampleid")
                }) %>%
                  setNames(setdiff(stim_vec_u, "uns")),
              )

              # load flowSet
              # -----------------
              fs <- load_fs_for_flowsom(
                dir_fcs = dir_fcs_exc_gn_fcs,
                stim = stim,
                stim_vec_u = stim_vec_u,
                responders_list = responders_list_resp_stim,
                conv_mtbaux_to_mtb = TRUE,
                chnl_sel = chnl_vec_sel
              )

              scale <- flowsom_scale[1]
              for (scale in flowsom_scale) {
                print(scale)

                n_clust <- flowsom_n_clust[1]
                for (n_clust in flowsom_n_clust) {
                  print(n_clust)

                  # ------------------------
                  # Fit FlowSOM and merge with other data
                  # ------------------------

                  dir_save_base <- file.path(
                    str_sub(dir_fcs_exc_gn_fcs, end = -5),
                    "fs", responders_only
                  )

                  dir_save_cluster_results <- file.path(
                    dir_save_base,
                    ifelse(stim %in% c("all", "all_u"),
                      paste0(stim, "-all"),
                      stim
                    ),
                    ifelse(scale, "scale", "no_s"),
                    n_clust
                  )
                  if (!dir.exists(dir_save_cluster_results)) {
                    dir.create(
                      dir_save_cluster_results,
                      recursive = TRUE
                    )
                  }

                  print("clustering")

                  flowsom_out <- fit_flowsom(
                    fs = fs, n_clust = n_clust,
                    scale = ifelse(scale, "scale", "no_s"),
                    dir_save_base = dir_save_base,
                    ds_name = ds_name,
                    jacc_reuse = jacc_reuse, jacc_n = jacc_n,
                    fs_seed = 2106, jacc_seed = 2107,
                    chnl_lab_vec = chnl_lab_vec_fcs,
                    tsne_tbl = tsne_list[[exc]][[gn]],
                    chnl_sel = chnl_vec_sel,
                    stim = stim
                  )

                  # ------------------------
                  # Get per-cluster frequencies
                  # ------------------------

                  print("getting per-cluster frequencies")
                  flowsom_freq <- get_fs_cluster_freq(
                    flowsom_out = flowsom_out,
                    dir_save_cluster_results = dir_save_cluster_results,
                    stats_combn_tbl = stats_combn_tbl,
                    stim = stim
                  )

                  # ------------------------
                  # Get per-cluster MIMOSA response probabilities
                  # ------------------------

                  print("getting per-cluster post probs")

                  if (flowsom_post_probs) {
                    flowsom_cluster_post_probs <- get_fs_cluster_post_probs(
                      flowsom_freq = flowsom_freq,
                      dir_save_cluster_results = dir_save_cluster_results,
                      reuse = TRUE
                    )
                  } else {
                    flowsom_cluster_post_probs <- NULL
                  }

                  # ------------------
                  # Summarise and plot flowsom results
                  # ------------------

                  print("making plots")

                  if (flowsom_plot) {
                    plot_flowsom(
                      flowsom_out = flowsom_out, stim = stim,
                      flowsom_freq = flowsom_freq,
                      flowsom_cluster_post_probs = flowsom_cluster_post_probs,
                      dir_save_base = dir_save_base,
                      dir_save_cluster_results = dir_save_cluster_results,
                      reuse = !flowsom_plot_force,
                      chnl_lab = chnl_lab_vec_fcs,
                      chnl_sel = chnl_vec_sel,
                      stats_combn_tbl = stats_combn_tbl %>%
                        dplyr::filter(gate_name == gn)
                    )
                  }

                  flowsom_tbl_curr <- flowsom_freq %>%
                    dplyr::select(-c(Progressor, timeToTB)) %>%
                    dplyr::mutate(
                      responders_only = responders_only,
                      exc = exc, gn = gn,
                      stim = stim, scale = scale, n_clust = n_clust
                    )

                  if (!is.null(flowsom_cluster_post_probs)) {
                    flowsom_tbl_curr %<>%
                      dplyr::left_join(
                        flowsom_cluster_post_probs,
                        by = c("clust", "stim", "SubjectID", "VisitType")
                      )
                  } else {
                    flowsom_tbl_curr %<>% dplyr::mutate(prob = NA)
                  }

                  flowsom_tbl_curr %<>%
                    dplyr::select(
                      responders_only:n_clust, stim,
                      clust, stability,
                      SubjectID, VisitType,
                      n_cell_tot_stim, n_cell_cyt_stim,
                      prob, count_stim:freq_tot_bs, everything()
                    ) %>%
                    dplyr::select(-fcs)

                  if (stim %in% c("all", "all_u")) {
                    flowsom_tbl_curr <- flowsom_tbl_curr %>%
                      dplyr::mutate(stim = paste0(.env$stim, "-", stim))
                  }

                  # output
                  if (!exists("flowsom_tbl")) {
                    flowsom_tbl <- flowsom_tbl_curr
                  } else {
                    flowsom_tbl <- flowsom_tbl %>%
                      dplyr::bind_rows(flowsom_tbl_curr)
                  }
                }
              }
            }
          }
        }
      }


      # save and add flowsom results
      # ----------------
      analysispipeline::save_objects(
        "flowsom" = flowsom_tbl,
        dir_proj = path_cyt_combn_base,
        dir_sub = c("data_out", "flowsom"),
        empty = TRUE
      )
    } else {
      flowsom_tbl <- readRDS(
        file.path(
          path_cyt_combn_base, "data_out",
          "flowsom", "flowsom.rds"
        )
      )
    }

    out_list <- out_list %>%
      append(list("flowsom" = flowsom_tbl))

    # save_final
    # --------------

    # prep
    dir_data_rda <- file.path(
      here::here(), "data"
    )
    if (!dir.exists(dir_data_rda)) {
      dir.create(dir_data_rda, recursive = TRUE)
    }
    path_data_rda <- file.path(
      dir_data_rda, paste0(ds_name, ".rda")
    )
    envir <- environment()
    assign(x = ds_name, value = out_list, envir = envir)
    save(
      list = ds_name,
      file = path_data_rda,
      version = 2,
      envir = envir,
      compress = "bzip2"
    )
    return(invisible(TRUE))
  }
