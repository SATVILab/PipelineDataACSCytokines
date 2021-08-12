#' @title Get directory in which FCS files are stored
#'
#' @param dir_gating character. Directory in which
#' gating results are stored, i.e. the last
#' part should be "gating".
#' @param pop 'cd4', 'cd8' or 'tcrgd'. Cell population.
#' @param chnl character vector. Channels that should be
#' gated together.
#' @param exc list of one list. Each element list should specify
#' the cytokines to be excluded.
#'
#' @export
#'
#' @examples
#' # no excluded cyt combns
#' get_dir_fcs_cyt(
#'   dir_gating = dir_cytofacs_gt,
#'   pop = 'cd4',
#'   chnl = c('Ho165Di', 'Gd158Di', 'Nd146Di', 'Dy164Di'),
#'   exc = list(NULL))
#' # exclude single TNF
#' get_dir_fcs_cyt(
#'   dir_gating = dir_cytofacs,
#'   pop = 'cd4',
#'   chnl = c('Ho165Di', 'Gd158Di', 'Nd146Di', 'Dy164Di'),
#'   exc = list(list("Nd146Di")))
get_dir_fcs_cyt <- function(dir_gating,
                            pop,
                            chnl,
                            exc){

  # prep
  # -------------

  dir_gs <- switch(
    pop,
    "cd4" = 'gs_cytof_acs_cd4',
    "cd8" = 'gs_cytof_acs_cd8',
    "tcrgd" = 'gs_cytof_acs_tcrgd'
  )
  dir_chnl <- paste0(chnl, collapse = "_")
  dir_exc <- purrr::map_chr(exc, function(exc_curr){
    if(is.null(exc_curr)) return("exc-none")
    paste0("exc-", paste0(exc_curr, collapse = "_&_"))
  })

  # out
  # -----------------

  file.path(
    dir_gating,
    dir_gs,
    "root",
    "ebvmtbp1p4uns",
    dir_chnl,
    dir_exc,
    "locb0.15_min_clust",
    "fcs"
  )

}
