% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_post_probs.R
\name{get_post_probs}
\alias{get_post_probs}
\title{Get posterior probabilities of response for each sample for each individual stim and per-stim list of responders}
\usage{
get_post_probs(
  stats_combn_tbl,
  stim,
  exc,
  uns_chr = "uns",
  conv_mtbaux_to_mtb = TRUE,
  gn = NULL,
  chnl_lab_cyt_custom = NULL,
  chnl_lab_cyt_fcs,
  post_prob_min
)
}
\arguments{
\item{stats_combn_tbl}{df.}

\item{stim}{'all_u', 'all', 'p1', 'p4', 'ebv' or 'mtbaux'. Stimulations for which posterior probabilities should
be calculated}

\item{exc}{character vector. Contains channel names of cytokine combination(s) separated by underscores ("_") 
that are not included as part of the summed response of cytokine-positive cells.}

\item{uns_chr}{character. Name for unstim in \code{stim_vec_u}. Default is "uns".}

\item{conv_mtbaux_to_mtb}{logical. If \code{TRUE}, then when looping over stimulations to 
calculate posterior probabilities, "mtbaux" is converted to "mtb" in the loop index.}

\item{gn}{character. Name of gate in \code{stats_combn_tbl}. If \code{NULL}, then set to "locb0.15_min_clust".
Default is \code{NULL}.}

\item{chnl_lab_cyt_custom}{named character vector. Converts cytokine chnl names to custom marker names. 
If \code{NULL}, then preset labelling vector is used. Default is \code{NULL}.}

\item{chnl_lab_cyt_fcs}{named character vector. Converts cytokine chnl names to FCS-derived marker names. No default, and
must be supplied.}

\item{post_prob_min}{0 to 1. Minimum posterior probability to count as a responder.}

\item{stim_vec_u}{character vector. Consists of names of individual stimulations. If \code{NULL}, then 
set to default value. Default is \code{NULL}.}
}
\value{
List with one named element for each stimulation, each containing
a list with named elements 'post_prob' and 'responders', which are the per-sample posterior probabilities
and the responding samples, respectively.
}
\description{
Uses MIMOSA.
}
