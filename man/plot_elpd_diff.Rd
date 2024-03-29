% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_elpd_diff}
\alias{plot_elpd_diff}
\title{Plot the ELPD difference interval}
\usage{
plot_elpd_diff(
  anpan_pglmm_res,
  probs = c(0.5, 0.95),
  color_category = NULL,
  verbose = TRUE
)
}
\arguments{
\item{anpan_pglmm_res}{a result from \code{anpan_pglmm} or \code{anpan_pglmm_batch}}

\item{probs}{the probability to cover in the inner and outer intervals}

\item{color_category}{an optional string giving the name of the column in the input for a
categorical variable used to color the intervals}

\item{verbose}{verbose}
}
\description{
Plot the ELPD difference interval
}
\details{
The validity of the intervals shown in the plot hinges on the normality approximation of
  the loo model comparison. See the [Cross validation
  FAQ](https://mc-stan.org/loo/articles/online-only/faq.html#se_diff) for more details.

  For batch results, you can set the \code{input_file} column to a factor to alter the vertical
  sorting of input files. By default it sorts according to ELPD difference.

  If you only want to highlight a subset of intervals with colors, set the \code{color_category}
  variable to NA for all other entries.
}
