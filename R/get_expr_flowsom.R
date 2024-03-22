#' @rdname get_expr_tbl
#' @title Get expression dat with dim-red, prob and FlowSOM dat
#' @inheritParams get_dir_fcs_cyt
#' @param ind integer vector. Indices of FCS files to load. Number of FCS files to load. If \code{NULL},
#' then all files are loaded. Default is \code{NULL}.
#'
#' @return A tibble with columns corresponding to marker names, with
#' columns specifying the FCS name and sample info.
#'
#' @export
get_expr_flowsom <- function(dir_gating,
                             pop,
                             chnl,
                             exc,
                             ind = NULL) {
  # prep
  # ---------------

  # directory housing FCS files
  dir_fcs_cyt <- get_dir_fcs_cyt(
    dir_gating = dir_gating,
    pop = pop,
    chnl = chnl,
    exc = exc
  )

  # prepare data
  # -------------------------

  # load
  fcs_vec <- list.files(dir_fcs_cyt)
  if (length(fcs_vec) == 0) {
    return(NULL)
  }
  if (!is.null(ind)) {
    fcs_vec <- fcs_vec[ind]
  }
  ex_tbl <- purrr::map(fcs_vec, function(x) {
    ex <- try(flowCore::exprs(flowCore::read.FCS(file.path(dir_fcs_cyt, x))) |>
      tibble::as_tibble() |>
      dplyr::mutate(fcs = x))
    if (identical(class(ex), "try-error")) {
      return(NULL)
    }
    ex
  }) |>
    purrr::compact() |>
    dplyr::bind_rows()

  # flowframe
  fr <- flowCore::read.FCS(file.path(dir_fcs_cyt, fcs_vec[1]))

  # format
  chnl_lab_vec_fcs <- c(UtilsCytoRSV::chnl_lab(fr), c("fcs" = "fcs"))
  marker_lab_vec_fcs <- c(UtilsCytoRSV::marker_lab(fr), c("fcs" = "fcs"))

  # remove columnns without descriptors and CD45
  ex_tbl <- ex_tbl[, colnames(ex_tbl) %in% names(chnl_lab_vec_fcs)]
  ex_tbl <- ex_tbl[, !colnames(ex_tbl) %in% (names(chnl_lab_vec_fcs)[chnl_lab_vec_fcs == "CD45"])]


  # rename columns
  colnames(ex_tbl) <- chnl_lab_vec_fcs[colnames(ex_tbl)]

  # transform columns
  ex_tbl <- ex_tbl |>
    dplyr::select(-fcs) |>
    dplyr::mutate_all(function(x) asinh(x / 5)) |>
    dplyr::mutate(fcs = ex_tbl$fcs)

  # get sample info from FCS name

  fcs_vec <- ex_tbl$fcs

  sample_info_tbl <- purrr::map_df(unique(fcs_vec), function(fcs) {
    split_vec <- stringr::str_split(stringr::str_sub(fcs, end = -5), "_")[[1]]
    tibble::tibble(
      SubjectID = split_vec[1],
      VisitType = split_vec[2],
      Progressor = split_vec[4],
      timeToTB = split_vec[3],
      stim = split_vec[5],
      fcs = fcs
    )
  })

  ex_tbl <- ex_tbl |>
    dplyr::mutate(fcs = fcs_vec) |>
    dplyr::left_join(sample_info_tbl, by = "fcs")

  ex_tbl
}


#' @rdname get_expr_tbl
ex_tbl_cyt <- get_expr_flowsom
