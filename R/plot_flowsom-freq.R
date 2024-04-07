#' @title Plot cluster frequencies per cyt+ and all cells, by stim and prog
#'
#' @param data,data_no_uns.
#' Whatever works with \code{.data} parameter of \code{saver::save_rds_eval}.
#' @param dir_save character. Directory to save to.
#'
#' @export
plot_freq <- function(data,
                      data_no_uns,
                      dir_save,
                      fill_col,
                      fill_lab) {
  if (dir.exists(dir_save)) {
    unlink(dir_save, recursive = TRUE)
  }
  if (!dir.exists(dir_save)) {
    dir.create(
      dir_save,
      recursive = TRUE
    )
  }

  # parameters
  # ------------------------
  facet_scale_vec <- c("free", "fixed", "asn")
  x_dodge_vec <- c("none", "stim")
  trans_vec_y <- c("identity", "asn")
  y_vec <- c(
    "freq_tot_bs", "freq_cyt_stim",
    "freq_cyt_bs", "freq_tot_stim"
  )

  arg_tbl <- expand.grid(
    facet_scale = facet_scale_vec,
    x_dodge = x_dodge_vec,
    trans_y = trans_vec_y,
    y = y_vec
  ) |>
    dplyr::tibble::as_tibble() |>
    dplyr::mutate_all(as.character)

  arg_tbl <- arg_tbl |>
    dplyr::filter(!(trans_y == "asn" & facet_scale == "fixed"))

  # plot
  # ------------------------

  purrr::walk(seq_len(nrow(arg_tbl)), function(i) {
    arg_row <- arg_tbl[i, ]
    axis_lab <- c(
      "x" = "Progressor status",
      "y" = switch(arg_row$y,
        "freq_cyt_stim" = "Frequency of cytokine-positive cells",
        "freq_cyt_bs" =
          "Background-subtracted frequency of\n cytokine-positive cells",
        "freq_tot" = "Frequency of all cells",
        "freq_tot_bs" = "Background-subtracted frequency of\n all cells"
      )
    )
    x_dodge <- switch(as.character(arg_row$x_dodge == "none"),
      "TRUE" = NULL,
      "FALSE" = arg_row$x_dodge
    )

    fill_col <- switch(as.character(arg_row$x_dodge == "none"),
      "TRUE" = fill_col[["prog"]],
      "FALSE" = fill_col[["stim"]]
    )

    fill_lab <- switch(as.character(arg_row$x_dodge == "none"),
      "TRUE" = fill_lab[["prog"]],
      "FALSE" = fill_lab[["stim"]]
    )
    .plot_freq(
      .data = switch(as.character(arg_row$y == "freq_tot" &&
        arg_row$x_dodge != "none"),
      "TRUE" = data_no_uns,
      "FALSE" = data
      ),
      x = "Progressor",
      y = arg_row$y,
      dir_save = dir_save,
      x_dodge = x_dodge,
      fill_col = fill_col,
      fill_lab = fill_lab,
      axis_lab = axis_lab,
      font_size = 14,
      facet_scale = arg_row$facet_scale,
      trans_y = arg_row$trans_y,
      width = NULL,
      height = NULL
    )
  })

  invisible(TRUE)
}


#' @title Plot cluster
.plot_freq <- function(.data,
                       x,
                       y,
                       dir_save,
                       x_dodge = NULL,
                       fill_col = NULL,
                       axis_lab = NULL,
                       font_size = 14,
                       facet_scale = c("free", "asn", "fixed"),
                       fill_lab,
                       trans_y = "identity",
                       n_col = NULL,
                       height,
                       width) {
  # create base plot to be faceted later in a variety
  # of ways
  path_rds <- saver::save_rds_eval(
    fn_or_call = .plot_freq_rds,
    .data = .data,
    filename = paste0("p-", y, "-", facet_scale),
    dir_save = dir_save,
    return_obj = FALSE,
    eval = FALSE,
    test = FALSE,
    x = x,
    y = y,
    x_dodge = x_dodge,
    fill_col = fill_col,
    fill_lab = fill_lab,
    axis_lab = axis_lab,
    font_size = font_size,
    facet_scale = facet_scale,
    trans_y = trans_y
  )

  p <- saver::load_rds_eval(path_rds)

  n_clust <- length(unique(p$data$clust))

  if (is.null(height)) {
    height <- 5 * ceiling(n_clust / 4)
  }
  if (is.null(width)) {
    width <- min(19.5, n_clust * 19.5 / 3)
  }

  purrr::walk(c("pdf", "png"), function(gd) {
    cowplot::ggsave2(
      filename = stringr::str_remove(path_rds, ".rds$") |>
        paste0(".", gd),
      plot = p,
      height = height,
      width = width,
      units = "cm"
    )
  })

  invisible(TRUE)
}

.plot_freq_rds <- function(.data,
                           x,
                           y,
                           x_dodge = x_dodge,
                           fill_col,
                           fill_lab,
                           axis_lab,
                           font_size,
                           facet_scale,
                           trans_y,
                           n_col = NULL,
                           height = NULL,
                           width = NULL) {
  .fill <- ifelse(is.null(x_dodge), x, x_dodge)
  n_clust <- length(unique(.data$clust))
  if (is.null(n_col)) {
    n_col <- pmin(n_clust, 4)
  }
  .data[[y]] <- pmax(0, .data[[y]])
  ggplot(
    data = .data,
    mapping = aes(x = .data[[x]], y = .data[[y]])
  ) +
    cowplot::theme_cowplot(font_size = font_size) +
    geom_hline(yintercept = 0) +
    cowplot::background_grid(major = "y") +
    switch(as.character(is.null(x_dodge)),
      "TRUE" = geom_boxplot(
        aes(fill = .data[[.fill]]),
        outlier.size = -1
      ),
      "FALSE" = geom_boxplot(
        aes(fill = .data[[.fill]]),
        outlier.size = -1,
        position = "dodge"
      )
    ) +
    switch(as.character(is.null(x_dodge)),
      "TRUE" = geom_point(size = 0.5),
      "FALSE" = geom_point(
        size = 0.5,
        group = .data[[.fill]],
        position = "dodge"
      )
    ) +
    scale_fill_manual(
      values = fill_col,
      labels = fill_lab
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank()
    ) +
    labs(
      x = axis_lab[1],
      y = axis_lab[2]
    ) +
    facet_wrap(~clust,
      scales = switch(facet_scale == "fixed",
        "fixed",
        "free"
      ),
      ncol = n_col
    ) +
    scale_y_continuous(trans = trans_y)
}
