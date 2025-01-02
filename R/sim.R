#' Draw a sample from a multivariate normal distribution
#'
#' @description Draw a sample from a multivariate normal distribution with
#'   specified mean and covariance.
#'
#' @param mean mean vector, defaults to 0 vector
#' @param Sigma covariance matrix
#' @param L optional cholesky factor of \code{Sigma}
#'
#' @details One of \code{Sigma} or \code{L} must be provided. Pre-computing L
#'   with \code{L = t(chol(Sigma))} will make calling this function repeatedly
#'   much faster.
#' @returns A vector giving one draw from the specified multivariate normal
#'   distribution.
#' @export
draw_mvnorm = function(mean = NULL, Sigma = NULL,  L = NULL) {

  is.null(L) && is.null(Sigma) && cli::cli_abort('One of Sigma or L must be provided.')

  L = L %||% t(base::chol(Sigma))

  d = nrow(L)

  mean = mean %||% rep(0, d)

  mean + drop(L %*% matrix(rnorm(d), ncol = 1))
}
