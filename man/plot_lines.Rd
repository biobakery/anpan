% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_lines}
\alias{plot_lines}
\title{Make a filter-labelled line plot}
\usage{
plot_lines(
  bug_file = NULL,
  meta_file = NULL,
  covariates,
  outcome,
  genomes_stats,
  omit_na = FALSE,
  fgf = NULL,
  bug_name = NULL,
  plot_ext = "pdf",
  plot_dir = NULL,
  subset_line = 200
)
}
\arguments{
\item{bug_file}{a path to a gene family file}

\item{meta_file}{a path to the corresponding metadata file}

\item{covariates}{covariates to account for (as a vector of strings)}

\item{outcome}{the name of the outcome variable}

\item{omit_na}{logical indicating whether to omit incomplete cases of the metadata}

\item{fgf}{a filtered gene family data frame}

\item{bug_name}{name of the bug}

\item{plot_ext}{extension to use for plots}

\item{subset_line}{integer gives the number of points to take along the lines}
}
\description{
Make a filter-labelled line plot
}
\details{
The required input is either \itemize{ \item{the gene family file and the metadata file}
  \item{OR a pre-filtered gene family file}} \code{subset_line} is used to make saving plots
  faster. Plotting thousands of lines each with tens of thousands of points along them is too
  much visual detail and makes saving the plot very slow. Set \code{subset_line} to 0 to turn off
  subsetting.
}
