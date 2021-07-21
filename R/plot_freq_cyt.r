#' @title Plot cluster 
.plot_freq_cyt <- function(data,
                           x,
                           y,
                           mult_stim
                           ) {

  # create base plot to be faceted later in a variety
  # of ways

  p_expr <- quote(
    ggplot(
      data = .data,
      mapping = aes(x = Progressor, y = freq_cyt_stim)
    ) +
      cowplot::theme_cowplot() +
      geom_hline(yintercept = 0) +
      cowplot::background_grid(major = "y") +
      geom_boxplot(aes(fill = Progressor), outlier.size = -1) +
      geom_point(size = 0.5) +
      scale_fill_manual(values = c(
        "Progressor" = "orange",
        "Control" = "dodgerblue"
      )) +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()
      ) +
      labs(
        x = "Progressor status",
        y = "Frequency of cytokine-positive cells"
      )
  )

  p_freq_cyt <- saver::save_rds_eval(
    fn_or_call = p_expr,
    .data = data,
    filename = "p_freq_cyt_base.rds",
    dir_save = dir_save,
    return_obj = TRUE,
    eval_fun = TRUE,
    test = FALSE
  )

  file.remove(
    file.path(
      dir_save,
      "p_freq_cyt_base.rds"
    )
  )

  p_list <- list(
    p_freq_cyt_free = p_freq_cyt +
      facet_wrap(~clust, scales = "free"),
    p_freq_cyt_free_asn = p_freq_cyt_free +
      scale_y_continuous(trans = "asn"),
    p_freq_cyt_fixed = p_freq_cyt +
      facet_wrap(~clust, scales = "fixed")
  )

  save_list(
    obj_list = p_list,
    ftype = c("rds", "pdf", "png"),
    width = 20,
    height = length(stability_lab_vec) * 1.25,
  )

}

#' @title Plot all stimulations together
.plot_freq_cyt_fill <- function(data,
                                fill) {

  p <- ggplot(
    data = data,
    mapping = aes(x = Progressor, y = freq_cyt_stim)
    ) +
    cowplot::theme_cowplot() +
    geom_hline(yintercept = 0) +
    cowplot::background_grid(major = "y")

  p <- p + switch(fill,
    "prog" = geom_boxplot(aes(fill = Progressor), outlier.size = -1) +
      scale_fill_manual(values = c(
        "Progressor" = "orange",
        "Control" = "dodgerblue"
      )),
    "stim" = geom_boxplot(
      aes(fill = stim),
      outlier.size = -1,
      position = "dodge") +
      scale_fill_manual(
        values = c(
          "uns" = "#a6cee3",
          "ebv" = "#1f78b4",
          "mtbaux" = "#b2df8a",
          "mtb" = "#b2df8a",
          "p1" = "#33a02c"
        ),
        labels = stim_lab_vec_full_break
      ) +
      theme(legend.position = "bottom")
  )

  p <- p +
    geom_point(size = 0.5) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank()
    ) +
    labs(
      x = "Progressor status",
      y = "Frequency of cytokine-positive cells"
    )
  p
}


