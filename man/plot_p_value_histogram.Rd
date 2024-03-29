% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_p_value_histogram}
\alias{plot_p_value_histogram}
\title{Make a p-value histogram}
\usage{
plot_p_value_histogram(
  all_bug_terms,
  out_dir = NULL,
  plot_ext = "pdf",
  n_bins = 50
)
}
\arguments{
\item{all_bug_terms}{a data frame of bug:gene glm fits}

\item{out_dir}{string giving the output directory}

\item{plot_ext}{string giving the extension to use}
}
\description{
This function makes a p-value histogram from a collection of
  bug:gene glm fits.
}
\details{
The plot will be written out to \code{p_value_histogram.<ext>}
  in the specified output directory. The "aaa" is there for alphabetical
  superiority.

  If you don't understand the purpose of this type of plot,
  \href{http://varianceexplained.org/statistics/interpreting-pvalue-histogram/}{this
  blog post by David Robinson} has a lot of helpful information.
}
