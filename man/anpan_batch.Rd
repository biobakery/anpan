% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modelling.R
\name{anpan_batch}
\alias{anpan_batch}
\title{Apply anpan to a many bugs}
\usage{
anpan_batch(
  bug_dir,
  meta_file,
  out_dir,
  genomes_dir = NULL,
  prefiltered_dir = NULL,
  model_type = "fastglm",
  covariates = c("age", "gender"),
  outcome = "crc",
  omit_na = FALSE,
  filtering_method = "kmeans",
  discretize_inputs = TRUE,
  minmax_thresh = NULL,
  discard_poorly_covered_samples = TRUE,
  skip_large = TRUE,
  save_fit = TRUE,
  annotation_file = NULL,
  save_filter_stats = TRUE,
  verbose = TRUE,
  plot_result = TRUE,
  plot_ext = "pdf",
  q_threshold = 0.1,
  beta_threshold = 1,
  n_top = 50,
  width = 10,
  height = 8,
  ...
)
}
\arguments{
\item{bug_dir}{a directory of gene family files}

\item{meta_file}{path to a metadata tsv}

\item{out_dir}{path to the desired output directory}

\item{genomes_dir}{an optional directory of genome files}

\item{prefiltered_dir}{an optional directory to pre-filtered data from an earlier run to skip the
filtering step}

\item{model_type}{either "horseshoe" or "fastglm"}

\item{covariates}{character vector of covariates to include in the model}

\item{outcome}{character string of the outcome variable}

\item{omit_na}{logical indicating whether to omit incomplete cases of the metadata}

\item{filtering_method}{method to use for filtering samples. Either "kmeans" or "none"}

\item{discretize_inputs}{logical indicating whether to discretize the input abundance
measurements (0/nonzero --> FALSE/TRUE) before passing them to the modelling function}

\item{minmax_thresh}{genes must have at least this many (or N - this many)
non-zero observations or else be discarded. NULL defaults to \code{floor(0.005*nrow(metadata))}.}

\item{discard_poorly_covered_samples}{logical indicating whether to discard samples where the
genes of a bug are poorly covered}

\item{skip_large}{logical indicating whether to skip bugs with over 5k genes. Only used when
model_type = "horseshoe".}

\item{save_fit}{logical indicating whether to save horseshoe fit objects. Only used when
model_type = "horseshoe".}

\item{annotation_file}{a path to a file giving annotations for each gene}

\item{save_filter_stats}{logical indicating whether to save filter statistics}

\item{plot_result}{logical indicating whether or not to plot the results}

\item{plot_ext}{extension to use for plots}

\item{q_threshold}{FDR threshold to use for inclusion in the plot.}

\item{beta_threshold}{Regression coefficient threshold to use for inclusion in the plot. Set to 0
to include everything.}

\item{n_top}{number of top elements to show from the results}

\item{width}{width of saved plot in inches}

\item{height}{height of saved plot in inches}

\item{...}{arguments to pass to [cmdstanr::sample()] if applicable}
}
\value{
a data.table of model statistics for each bug:gene combination
}
\description{
This function calls anpan() on each gene family file in the \code{bug_dir} directory
  and makes a composite data + results plot for each.
}
\details{
\code{bug_dir} should be a directory of gene (or SNV or pathway) abundance files, one
  for each bug.

  \code{annotation} file must have two columns named "gene" and "annotation"

  See \code{?anpan()} for the format / usage if providing genome files. If provided, genomes_dir
  must contain ONLY the genome files themselves.
}
\seealso{
[anpan()]
}
