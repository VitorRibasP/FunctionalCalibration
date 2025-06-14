#' @title Weight Estimation
#'
#' @description
#' Estimates the weights associated with the functional coefficients
#' \eqn{\alpha_l(x)} using the  using Ordinary Least Squares.
#'
#' The problem can be formulated as:
#'
#' \deqn{A(x) = \displaystyle \sum_{l=1}^{L} y_l \alpha_l(x)}
#'
#' where \eqn{A(x)} is the aggregated function evaluated at each point \eqn{x},
#' \eqn{\alpha_l(x)} are the functional coefficients, and \eqn{y_l} are the weights to be estimated.
#'
#' @param data A numeric vector representing one sample of the aggregated function \eqn{A(x)}, evaluated at a grid of points \eqn{x}.
#' @param alpha A numeric matrix where each column represents the values of a function \eqn{\alpha_l(x)} evaluated at the same grid of points as \code{data}.
#'
#' @return The function returns a vector with the estimated weights obtained using Ordinary Least Squares.
#'
#' @examples
#' weight_estimation(simulated_data$data[,1], simulated_data$alphas)
#'
#' @export
weight_estimation <- function(data, alpha) {
  if (class(alpha) == "data.frame") {
    alphas <- as.matrix(alpha)
  }
  coefs <- coef(lm(data ~ . - 1, cbind(data, alpha)))
  if (!(is.null(colnames(alpha)))) {
    names(coefs) <- colnames(alpha)
  }
  return(
    coefs
  )
}
