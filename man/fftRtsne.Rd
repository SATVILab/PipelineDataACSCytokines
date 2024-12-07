% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fast_tsne.R
\name{fftRtsne}
\alias{fftRtsne}
\title{Compute FIt-SNE of a dataset}
\usage{
fftRtsne(
  X,
  dims = 2,
  perplexity = 30,
  theta = 0.5,
  max_iter = 750,
  fft_not_bh = TRUE,
  ann_not_vptree = TRUE,
  stop_early_exag_iter = 250,
  exaggeration_factor = 12,
  no_momentum_during_exag = FALSE,
  start_late_exag_iter = -1,
  late_exag_coeff = 1,
  mom_switch_iter = 250,
  momentum = 0.5,
  final_momentum = 0.8,
  learning_rate = "auto",
  n_trees = 50,
  search_k = -1,
  rand_seed = -1,
  nterms = 3,
  intervals_per_integer = 1,
  min_num_intervals = 50,
  K = -1,
  sigma = -30,
  initialization = "pca",
  max_step_norm = 5,
  data_path = NULL,
  result_path = NULL,
  load_affinities = NULL,
  fast_tsne_path = NULL,
  nthreads = 0,
  perplexity_list = NULL,
  get_costs = FALSE,
  df = 1,
  FAST_TSNE_SCRIPT_DIR = dirname(.libPaths()[2])
)
}
\arguments{
\item{dims}{dimensionality of the embedding. Default 2.}

\item{perplexity}{perplexity is used to determine the
bandwidth of the Gaussian kernel in the input
space.  Default 30.}

\item{theta}{Set to 0 for exact.  If non-zero, then will use either
Barnes Hut or FIt-SNE based on nbody_algo.  If Barnes Hut, then
this determins the accuracy of BH approximation.
Default 0.5.}

\item{max_iter}{Number of iterations of t-SNE to run.
Default 750.}

\item{fft_not_bh}{if theta is nonzero, this determines whether to
use FIt-SNE or Barnes Hut approximation. Default is FIt-SNE.
set to be True for FIt-SNE}

\item{ann_not_vptree}{use vp-trees (as in bhtsne) or approximate nearest neighbors (default).
set to be True for approximate nearest neighbors}

\item{stop_early_exag_iter}{When to switch off early exaggeration.
Default 250.}

\item{exaggeration_factor}{coefficient for early exaggeration
(>1). Default 12.}

\item{no_momentum_during_exag}{Set to 0 to use momentum
and other optimization tricks. 1 to do plain,vanilla
gradient descent (useful for testing large exaggeration
coefficients)}

\item{start_late_exag_iter}{When to start late exaggeration. 'auto' means
that late exaggeration is not used, unless late_exag_coeff>0. In that
case, start_late_exag_iter is set to stop_early_exag_iter. Otherwise,
set to equal the iteration at which late exaggeration should begin.
  Default 'auto'}

\item{late_exag_coeff}{Late exaggeration coefficient.
Set to -1 to not use late exaggeration.
 Default -1}

\item{learning_rate}{Set to desired learning rate or 'auto', which
sets learning rate to N/exaggeration_factor where N is the sample size, or to 200 if
     N/exaggeration_factor < 200.
         Default 'auto'}

\item{nterms}{If using FIt-SNE, this is the number of
               interpolation points per sub-interval
intervals_per_integer - See min_num_intervals}

\item{min_num_intervals}{Let maxloc = ceil(max(max(X)))
and minloc = floor(min(min(X))). i.e. the points are in
a [minloc]^no_dims by [maxloc]^no_dims interval/square.
The number of intervals in each dimension is either
min_num_intervals or ceil((maxloc -
minloc)/intervals_per_integer), whichever is
larger. min_num_intervals must be an integer >0,
and intervals_per_integer must be >0. Default:
min_num_intervals=50, intervals_per_integer =
1}

\item{K}{Number of nearest neighbours to get when using fixed sigma
Default -30 (None)}

\item{sigma}{Fixed sigma value to use when perplexity==-1
Default -1 (None)}

\item{initialization}{'pca', 'random', or N x no_dims array to intialize the solution.
Default: 'pca'}

\item{max_step_norm}{Maximum distance that a point is allowed to move on
one iteration. Larger steps are clipped to this value. This prevents
possible instabilities during gradient descent.  Set to -1 to switch it
off. Default: 5}

\item{load_affinities}{If 1, input similarities are loaded from a file and not computed
If 2, input similarities are saved into a file.
If 0, affinities are neither saved nor loaded}

\item{perplexity_list}{if perplexity==0 then perplexity combination will
be used with values taken from perplexity_list. Default: NULL}

\item{df}{Degree of freedom of t-distribution, must be greater than 0.
Values smaller than 1 correspond to heavier tails, which can often 
resolve substructure in the embedding. See Kobak et al. (2019) for
details. Default is 1.0}

\item{FAST_TSNE_SCRIPT_DIR}{character. 
Directory to where FITsne binaries have been saved to.
Default is \code{dirname(.libPaths()[2])}.}
}
\description{
Compute FIt-SNE of a dataset
}
\details{
Note that original version is at \code{base::system.file("extdata", "fast_tsne.R", package = 'PipelineDataACSCytokines')}. 
The only difference is that \code{FAST_TSNE_SCRIPT_DIR} has been made a function parameter with a default (see parameter description for details).
}
