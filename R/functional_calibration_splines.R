#' @title Functional Data Calibration with Splines
#'
#' @description
#' This function performs functional calibration based on the following model:
#'
#' \deqn{A_i(x_m) = \displaystyle \sum_{l=1}^{L} y_{il} \alpha_l(x_m) + e_i(x_m), \quad i = 1,...,I, \quad m = 1,...,M = 2^J}
#'
#' where the functions \eqn{\alpha_l(x)} are estimated using spline basis functions.
#'
#' In matrix notation, the model is represented as:
#'
#' \deqn{A = \alpha y + e}
#'
#' @param data A matrix \eqn{M} x \eqn{I} where each column represents one sample of the aggregated function — the matrix \eqn{A} in the model.
#' @param weights A matrix \eqn{L} x \eqn{I} representing the weight values associated with each sample — the matrix \eqn{y} in the model.
#' @param x A numeric vector of values at which the function is evaluated.
#' @param n_functions Number of spline basis functions to be used for estimating \eqn{\alpha_l(x)}.
#'
#' @return
#'
#' The function returns a list containing two objects.
#'
#' \describe{
#'   \item{\code{alpha}}{A matrix with the estimated functional coefficients \eqn{\alpha}.}
#'   \item{\code{Plots}}{A list of plot objects, each representing the corresponding function \eqn{\alpha_l(x)}.}
#'}
#'
#' @references Saraiva, M. A., & Dias, R. (2009). Analise não-parametrica de dados funcionais: uma aplicação a quimiometria (Doctoral dissertation, Master’s thesis, Universidade Estadual de Campinas, Campinas).
#'
#' @examples
#' functional_calibration_splines(simulated_data$data, simulated_data$weights, simulated_data$x)
#' functional_calibration_splines(simulated_data$data, simulated_data$weights, simulated_data$x, 12)
#'
#' @importFrom splines bs
#' @importFrom grDevices recordPlot
#'
#' @export
functional_calibration_splines <- function(data, weights, x,
                                           n_functions = 10) {

  if (ncol(data) != ncol(weights)) {
    stop("Error: data dimension is not compatible with weights dimension.")
  }

  if (inherits(data, "data.frame")) {
    data <- as.matrix(data)
  }

  if (inherits(weights, "data.frame")) {
    weights <- as.matrix(weights)
  }

  X <- as.vector(data)
  d <- matrix(NA, nrow = ncol(data)*length(x), ncol = nrow(weights))
  for (i in 1:nrow(weights)) {
    d[,i] <- rep(weights[i,], each = length(x))
  }

  B <- bs(x, degree = 3, df = n_functions)
  D <- d*do.call(rbind, replicate(ncol(data),
                                  do.call(cbind, replicate(nrow(weights), B[,1], simplify = FALSE))
                                  , simplify = FALSE))

  for (i in 2:n_functions) {
    D <- cbind(D, d*do.call(rbind, replicate(ncol(data),
                                             do.call(cbind, replicate(nrow(weights), B[,i], simplify = FALSE))
                                             , simplify = FALSE)))
  }

  theta <- solve(t(D)%*%D)%*%t(D)%*%X

  alpha <- matrix(NA, nrow = length(x), ncol = nrow(weights))
  for (i in 1:ncol(alpha)) {
    thetas <- theta[seq(i, nrow(theta), by = ncol(alpha))]
    alpha[,i] <- B%*%thetas
  }
  if(is.null(rownames(weights))) {
    titles <- paste0("alpha", 1:nrow(weights))
  } else {
    titles <- rownames(weights)
  }

  Plots <- list()

  for (i in 1:ncol(alpha)) {
    plot(x, alpha[,i], type = "l", lwd = 2,
         main = paste0("Estimated ", titles[i]), xlab = "x", ylab = "y")
    Plots[[i]] <- recordPlot()
  }

  return(list("alpha" = alpha, "Plots" = Plots))
}
