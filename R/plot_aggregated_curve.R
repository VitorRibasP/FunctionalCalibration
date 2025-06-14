#' @title Aggregated Curve Plot
#'
#' @description
#'
#' Generates the plot of the aggregated curve based on the functional coefficients and their corresponding weights.
#' The aggregated curve is computed as:
#'
#' \deqn{A(x) = \displaystyle \sum_{l=1}^{L} y_l \alpha_l(x)}
#'
#' @param alpha A numeric matrix where each column represents the values of a function \eqn{\alpha_l(x)} evaluated at each point in \eqn{x}.
#' @param weights A numeric vector with the weight values corresponding to each function \eqn{\alpha_l(x)}.

#' @param title A string specifying the title of the plot.
#' @param x A numeric vector of values at which the function is evaluated. If \code{NULL}, the default is the sequence \code{1:nrow(alpha)}.
#'
#' @return The function returns the plot of the aggregated function.
#'
#' @examples
#' plot_aggregated_curve(simulated_data$alphas, c(0.7, 0.3))
#' plot_aggregated_curve(simulated_data$alphas, c(0.7, 0.3), "Aggregated Curve Example", simulated_data$x)
#'
#' @export
plot_aggregated_curve <- function(alpha, weights, title = NULL, x = NULL) {
  if (is.null(x)) {
    x <- 1:nrow(alpha)
  }
  if (class(alpha) == "data.frame") {
    alpha <- as.matrix(alpha)
  }
  plot(x, alpha%*%weights, type = "l", xlab = "x", ylab = "y", main = title, lwd = 2)
  Plot <- recordPlot()
  return(Plot)
}
