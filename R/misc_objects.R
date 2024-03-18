#' @title FAUST project directory
#'
#' @export
# dir_data <- file.path(here::here(), "..", "..", "Data")

#' @title Cytokine labelling function
#'
#' @description Function to label cytokines
#'
#' @export
label_cytokines <- function(cyt) {
  purrr::map(cyt, function(cyt_ind) {
    switch(cyt_ind,
      "CCR4-Beads" = "CCR4",
      "beads" = "Beads",
      "IFNg-beads" = ,
      "IFNg" = bquote(paste(plain(paste("IFN")), gamma)),
      "Cd11c" = "CD11c",
      "HLA-DR-beads" = ,
      "HLA-DR" = "HLADR",
      "MIP1b" = bquote(paste(plain(paste("MIP1")), beta)),
      "Perf-beads" = ,
      "Perf" = "Perforin",
      "TCRgd-CD19" = bquote(paste(plain(paste("TCR")), gamma, delta, plain(paste("-CD19")))),
      "TNFa" = "TNF",
      "Va7.2" = bquote(paste(plain(paste("V")), alpha, plain(paste("7.2")))),
      cyt_ind
    )
  })
}
