#' @title Save objects in a list
#'
#' @description Save multiple objects of multiple classes
#' in multiple formats with less error-prone naming
#' method of using a named list.
#'
#' @param ... name-value pairs. Saved.
#' @param obj_list named list. If not named, then names
#' will be provided.
#' @param ftype "pdf", "png" and/or "rds". Save format.
#' @param height,width,scale numeric. Passed to \code{cowplot::ggsave}.
#' @param dir_save directory. To save files to.
#'
#' @return \code{invisible(TRUE)}.
#'
#' @export
save_list <- function(...,
                      obj_list,
                      ftype,
                      height = 6,
                      width = 6,
                      scale = 1,
                      dir_save) {

  if (!dir.exists(dir_save)) {
    dir.create(dir_save, recursive = TRUE)
  }

  p_dots <- list(...)
  obj_list <- p_dots %>%
    append(obj_list)

  if (class(obj_list)[1] == "gg") {
    obj_list <- list("p" = obj_list)
  }

  ftype_vec_allowed <- c("rds", "png", "pdf")
  ftype_vec_allowed <- ftype_vec_allowed %>%
    c(paste0(".", ftype_vec_allowed))
  ftype <- ftype[ftype %in% ftype_vec_allowed]
  if (any(grepl("rds", ftype))) {
    fn_rds <- names(obj_list)[i]
    fn_rds <- switch(as.character(grepl(".rds$")),
      "TRUE" = fn_rds,
      "FALSE" = paste0(fn_rds, ".rds")
    )
    saveRDS(
      object = obj_list[[i]],
      file = file.path(
        dir_save,
        names(obj_list)[i]
      )
    )
  }

  ftype <- stringr::str_remove(ftype, "^.")
  ftype_gd_vec <- ftype[ftype %in% c("pdf", "png")]
  for (gd_ind in seq_along(ftype_gd_vec)) {
    gd <- ftype_gd_vec[gd_ind]
    for (i in seq_along(obj_list)) {
      if (!class(obj_list[[i]])[1] == "gg") next
      cowplot::ggsave2(
        filename = file.path(
          dir_save,
          paste0(
            names(obj_list)[i],
            gd
          )
        ),
        plot = obj_list[[i]],
        width = weight[i],
        height = height[i],
        units = "cm",
        scale = scale
      )
    }
  }

  return(invisible(TRUE))
}


.ensure_named <- function(obj) {
  if (!is.vector(obj) && !is.list(obj)) {
    stop(paste0(
      "obj must be a vector or a list"
    ))
  }
  nm_vec <- names(obj)
  unnamed_vec_ind <- which(nm_vec == "")
  nm_vec[unnamed_vec_ind] <- paste0("v", seq_len(length(unnamed_vec_ind)))
  names(obj) <- nm_vec
  obj
}

.remove_slashes <- function(obj) {
  if (!is.vector(obj) && !is.list(obj)) {
    stop(paste0(
      "obj must be a vector or a list"
    ))
  }
  nm_vec <- names(obj)
  nm_vec <- purrr::map_chr(nm_vec, function(x) {
    stringr::str_replace_all(x, "/", "_")
  })
  nm_vec <- purrr::map_chr(
    nm_vec, function(x) stringr::str_replace_all(x, "\\\\", "_")
  )
  names(obj) <- nm_vec
  obj
}

.ensure_unique <- function(obj) {
  if (!is.vector(obj) && !is.list(obj)) {
    stop("obj must be a vector or a list")
  }
  nm_vec <- names(obj)
  uni_vec <- unique(nm_vec)
  nm_vec_rep <- nm_vec
  for (i in seq_along(uni_vec)) {
    nm <- uni_vec[i]
    ind_vec <- which(nm_vec == nm)
    if (length(ind_vec) > 1) {
      nm_rep <- paste0(nm, "_", seq_len(length(ind_vec)))
      nm_vec_rep[ind_vec] <- nm_rep
    }
  }
  names(obj) <- nm_vec_rep
  obj
}

.fix_names <- function(obj) {
  obj %>%
    .ensure_named() %>%
    .remove_slashes() %>%
    .ensure_unique()
}
