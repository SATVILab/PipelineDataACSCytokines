cd8_full_vec_bright <- c('CD3' = 2, "CD8-IgD" = 3, "CD4" = 1,
                  "CD33" = 1, "CD20" = 1, "CD14" = 1,
                  "TCRgd-CD19" = 1) # yes
cd8_full_vec_dim <- c('CD3' = 2, "CD8-IgD" = 2, "CD4" = 1,
                  "CD33" = 1, "CD20" = 1, "CD14" = 1,
                  "TCRgd-CD19" = 1) # yes
cd4_full_vec <- c('CD3' = 2, "CD4" = 2, "CD8-IgD" = 1,
                  "CD33" = 1, "CD20" = 1, "CD14" = 1,
                  "TCRgd-CD19" = 1) # yes

tcrgd_vec <- c('CD3' = 2, "TCRgd-CD19" = 2, "CD33" = 1, "CD20" = 1, "CD14" = 1) # yes

nkt_vec <- c('CD3' = 2, "TCRgd-CD19" = 1, "CD33" = 1, "CD20" = 1, "CD14" = 1,
             'CD16' = 2)

mait_vec <- c('CD3' = 2, "TCRgd-CD19" = 1, "CD33" = 1, "CD20" = 1, "CD14" = 1,
              'CD16' = 1, 'CD161' = 2)

# cd4_basic_vec <- c('CD3' = 2, "CD4" = 2,  "CD8-IgD" = 1)
# cd8_basic_vec <- c('CD3' = 2, "CD8-IgD" = 2, "CD4" = 1)
pre_nk_vec <- c('CD3' = 1, 'CD33' = 1, 'CD14' = 1) # yes

nk_vec <- c('CD3' = 1, 'CD33' = 1, 'CD7' = 2, 'CD14' = 1, 'CD20' = 1, 'TCRgd-CD19' = 1) # yes

monocyte_vec <- c('CD3' = 1, 'CD33' = 2, 'HLA-DR-beads' = 2, 'TCRgd-CD19' = 1,
                  "CD8-IgD" = 1, 'CD14' = 2, "CD20" = 1)
bcell_vec_cd19 <- c('CD3' = 1, 'CD33' = 1, 'CD7' = 1,
                    'HLA-DR-beads' = 2, 'CD14' = 1, 'TCRgd-CD19' = 2) # yes


pop_list <- list('cd8_bright' = cd8_full_vec_bright,
                 'cd8_dim' = cd8_full_vec_dim,
                 'cd4' = cd4_full_vec,
                 'tcrgd' = tcrgd_vec, #'nkt' = nkt_vec,
                 'nk' = nk_vec, 'monocyte' = monocyte_vec,
                 'bcell' =  bcell_vec_cd19)

dir_base_tcell <- file.path(
  dirname(here::here()), 
  "DataTidyACSCyTOFCytokinesTCells", 
  "datalarge"
)
dir_base_nkbcell <- file.path(
  dirname(here::here()), 
  "DataTidyACSCyTOFCytokinesNKBCells", 
  "datalarge"
)

pop_to_path_pop <- c(
  "cd8_bright" = file.path(
    dir_base_tcell, 
    "gs_cytof_acs_cd8", 
    "root", 
    "ebvmtbp1p4uns"
  ),
  "cd8_dim" = file.path(
    dir_base_tcell, 
    "gs_cytof_acs_cd8", 
    "root", 
    "ebvmtbp1p4uns"
  ),
  "cd4" = file.path(
    dir_base_tcell, 
    "gs_cytof_acs_cd4", 
    "root", 
    "ebvmtbp1p4uns"
  ),
  "tcrgd" = file.path(
    dir_base_tcell, 
    "gs_cytof_acs_tcrgd", 
    "root", 
    "ebvmtbp1p4uns"
  )
)

cyt_vec <- c("IFNg-beads", "IL2", "IL6", "IL17", "IL22", "TNFa")
pop_to_gate_tbl <- purrr::map(pop_to_path_pop, function(path_pop) {
  purrr::map_df(cyt_vec, function(cyt){
    # get base directory
    path_cyt <- file.path(path_pop, cyt)
    # get stats tbl
    gate_tbl <- readRDS(file.path(path_cyt, "gate_tbl.rds"))
    
    gate_tbl %>%
      dplyr::mutate(cut = cyt) %>%
      dplyr::select(marker, everything()) %>%
      dplyr::filter(gate_combn == "min_clust")
  }) 
}) %>% 
  setNames(names(pop_to_path_base_gate))

get_faust_phenotype_cyt_expr <- function(dir_faust, pop_list, chnl, 
                                         gate_tbl
                                         ) {
  dir_faust <- file.path(
    dirname(here::here()),
    "DataTidyACSCyTOFFAUST",
    "datalarge",
    "visc",
    "faustData"
  )

  # take faust sampleData
  fcs_vec_faust <- list.dirs(
    file.path(dir_faust, "sampleData"),
    recursive = FALSE,
    full.names = FALSE
  )
  
  orig_to_tidy_fcs <- setNames(
    DataTidyACSCyTOFFAUST::clean_fcs_for_matching(fcs_vec_faust), 
    fcs_vec_faust
  )

  dir_sample <- file.path(
    dir_faust, 
    "sampleData", 
    fcs_vec_faust[1]
  )
  
  ann_tbl <- readRDS(
    file.path(
      dir_sample, 
      "resMat.rds"
    )
  ) %>%
    tibble::as_tibble()
  
  ann_tbl_faust <- readr::read_csv(
    file.path(
      dir_sample, 
      "faustAnnotation.csv"
    ), 
    col_names = FALSE, 
    col_types = "c"
  )
  
  ex_tbl <- readRDS(
    file.path(
      dir_sample, 
      "exprsMat.rds"
    )
  ) %>%
    tibble::as_tibble()

  # for each faust sample:
  # combine with FAUST phenotypic annotation
  pop <- "cd4"
  
  ind_vec <- rep(TRUE, nrow(ex_tbl))
  pop_defn <- pop_list[[pop]]

  ind_list <- purrr::map(seq_along(pop_defn), function(i) {
    pttrn <- paste0(names(pop_defn)[i], "~", pop_defn[i])
    grepl(pttrn, ann_tbl_faust[[1]])
  }) 

  ind_vec <- purrr::reduce(
    .x = ind_list, 
    .f = function(x, y) x & y
  )
  
  ex <- ex_tbl[ind_vec, ]
  
  gate_tbl_pop <- pop_to_gate_tbl[[pop]]
  
  ann_list_cyt <- purrr::map(cyt) 
  
  
  # get indices mapped to sample ids (but not stims)
  # ---------------------------
  dir_gs <- normalizePath(
    file.path(dir_faust, 
              "gsData",
              paste0("gs_cytof_acs_", pop)))
  gs <- flowWorkspace::load_gs(dir_gs)
  
  pop_gate <- "root"
  ind_in_batch_lab_vec <- c("1" = "ebv", "2" = "mtb", "3" = "p1",
                            "4" = "p4", "5" = "uns")
  
  fcs_stim_tbl <- purrr::map_df(seq_along(gs), function(i) {
    fn <- flowCore::description(
      flowWorkspace::gh_pop_get_data(gs[[i]])
    )$FILENAME
    last_slash_loc_tbl <- stringr::str_locate_all(fn, "/")[[1]]
    last_slash_loc <- last_slash_loc_tbl[nrow(last_slash_loc_tbl), "end"][[1]]
    fn <- stringr::str_sub(fn, last_slash_loc + 1) %>%
      stringr::str_remove(".fcs")
    second_underscore_loc <- stringr::str_locate_all(fn, "_")[[1]][2, "end"][[1]]
    fcs <- stringr::str_sub(fn, end = second_underscore_loc - 1)
    stim <- stringr::str_split(fn, "_")[[1]]
    stim <- stim[[length(stim)]]
    tibble::tibble(fcs = fcs, stim = stim)
  })
  
  fcs_vec_gs <- fcs_stim_tbl$fcs
  
  ind_batch_list <- outliergatev1:::.get_ind_batch_list(
    data = gs, batch_size = 5,
    cytof_fcs_to_clin_map = NULL,
    fcs = fcs_vec_gs
    )
  
  .get_batch_for_ind <- function(ind) {
    names(ind_batch_list)[purrr::map_lgl(ind_batch_list, function(x) ind %in% x)]
  }
  .get_batch_for_ind(1)
  ind_in_batch <- ifelse(ind %% 5 == 0, 5, ind - (ind %/% 5) * 5)
  
  stim <- ind_in_batch_lab_vec[as.character(ind_in_batch)]
  
  fcs_vec_faust[1]
  
  match_tbl <- DataTidyACSCyTOFFAUST::fcs_to_sampleid_map %>%
    dplyr::filter(MatchFCSName == DataTidyACSCyTOFFAUST::clean_fcs_for_matching(fcs_vec_faust[1])) %>%
    dplyr::mutate(SampleID = paste0(SubjectID, "_", VisitType))
  sample_id <- match_tbl$SampleID 
  stim <- match_tbl$Stim
  gate_tbl_ind <- gate_tbl_pop %>%
    dplyr::filter(.data$ind == .env$ind) 

  
  
  
  # identify cells that belong to each main cell type

  # gate them using pre-defined cytokine gates

  # get out:
  # - count matrices (total of cell type, total of phenotype, total of each cyt combn)
  # - FCS files (of what? cyt+ cells? yes, I think so).

  # final count matrix:
  # each row is one sample for one stim, for a given phenotype showing the count

  # then apply COMPASS, get per cyt_combn prob of a response
  # get overall prob of a response
  n_cell <- nrow(ex)
  
  pos_tbl <- purrr::map(cyt_vec, function(cyt) setNames(tibble::tibble(x = rep(0, n_cell)), cyt)) %>%
    dplyr::bind_cols()
  
  inc_list <- purrr::map(cyt_vec, function(cyt) {
    inc_vec <- ex[[cut_curr]] > gate_tbl[gate_tbl$cut == cut_curr & gate_tbl$ind == ind, 'gate', drop = TRUE]
    pos_tbl[[cut_curr]] <- inc_vec %>% as.numeric()
  })
  
  inc_vec <- purrr::reduce(
    .x = inc_list, 
    .f = function(x, y) x & y
  )
  

  
  for(cut_curr in cut){
    
  }
  
  ind_list <- map(seq_along(cut), function(i) 0:1)
  combn_mat <- expand.grid(ind_list)
  all_neg_ind_vec <- map_lgl(1:nrow(combn_mat), function(i){
    map_lgl(1:ncol(combn_mat), function(j) !combn_mat[i, j, drop = TRUE]) %>%
      all
  })
  combn_mat <- combn_mat[!all_neg_ind_vec,]
  combn_tbl <- combn_mat %>%
    tibble::as_tibble() %>%
    setNames(cut) %>%
    dplyr::mutate(count = NA_integer_)
  for(i in 1:nrow(combn_tbl)){
    ind_vec <- rep(TRUE, nrow(pos_tbl))
    for(j in 1:ncol(combn_mat)){
      ind_vec <- ind_vec & pos_tbl[[j]] == combn_tbl[i,j, drop = TRUE]
    }
    combn_tbl[i,"count"] <- sum(ind_vec)
  }


  }
