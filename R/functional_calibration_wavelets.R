#' @title Logistic Density
#'
#' @description
#' Computes the function:
#'
#' \deqn{g(\theta; \tau) = \displaystyle \frac{\exp\left(-\frac{\theta}{\tau}\right)}{\tau \left(1 + \exp\left(-\frac{\theta}{\tau}\right)\right)^2}}
#'
#' @param theta Numeric value of \eqn{\theta}.
#' @param tau Numeric value of \eqn{\tau}.
#'
#' @return A numeric value representing the result of the function \eqn{g(\theta; \tau)} for the specified inputs.
#'
#' @examples
#' Logistic_Density(0.3, 0.5)
#'
#'@keywords internal
Logistic_Density <- function(theta, tau) {
  return(exp(-theta/tau)/(tau*(1 + exp(-theta/tau))^2))
}
#' @title Bayesian Shrinkage
#'
#' @description
#' A Bayesian shrinkage method applied to empirical coefficients \eqn{d}, aiming to denoise them.
#'
#' The shrinkage function is defined as:
#' \deqn{\delta(d) = \displaystyle \frac{(1 - p) \int_{\mathbb{R}} (\sigma u + d) \, g(\sigma u + d; \tau) \, \phi(u) \, du}{\frac{p}{\sigma} \phi\left( \frac{d}{\sigma} \right) + (1 - p) \int_{\mathbb{R}} g(\sigma u + d; \tau) \, \phi(u) \, du}}
#'
#' where \eqn{\phi(x)} is the probability density function of the standard normal distribution,
#' and \eqn{g(\theta; \tau)} is the logistic density function.
#'
#' @param d Numeric value of the empirical coefficient to be denoised.
#' @param tau Numeric value of \eqn{\tau}.
#' @param p Numeric value of \eqn{p}.
#' @param sigma Numeric value of \eqn{\sigma}.
#' @param MC A logical evaluating to \code{TRUE} or \code{FALSE} indicating if the integrals will be approximated using Monte Carlo.
#'
#' @return A numeric value representing the result of the Bayesian shrinkage applied to the empirical coefficient \eqn{d}.
#'
#' @examples
#' Bayesian_Shrinkage(0.3, 5, 0.95, 0.1, FALSE)
#'
#'@keywords internal
Bayesian_Shrinkage <- function(d, tau, p, sigma, MC = FALSE) {
  if (!MC) {
    N <- rnorm(500)
    return(
      ((1-p)*mean((sigma*N + d) * Logistic_Density(sigma*N + d, tau)))/
        (p/sigma*dnorm(d/sigma) + (1-p)*mean(Logistic_Density(sigma*N + d, tau)))
    )
  } else {
    return(
      ((1-p)*integrate(function(u) (sigma*u + d) * Logistic_Density(sigma*u + d, tau) * dnorm(u), -Inf, Inf)$value)/
        (p/sigma*dnorm(d/sigma) + (1-p)*integrate(function(u) Logistic_Density(sigma*u + d, tau) * dnorm(u), -Inf, Inf)$value)
    )
  }
}
#' @title Functional Data Calibration with Wavelets
#'
#' @description
#' This function performs functional calibration based on the following model:
#'
#' \deqn{A_i(x_m) = \displaystyle \sum_{l=1}^{L} y_{il} \alpha_l(x_m) + e_i(x_m), \quad i = 1,...,I, \quad m = 1,...,M = 2^J}
#'
#' where the functions \eqn{\alpha_l(x)} are estimated using wavelet decomposition.
#'
#' In matrix notation, the model is represented as:
#'
#' \deqn{A = \alpha y + e}
#'
#' @param data A matrix \eqn{M} x \eqn{I} where each column represents one sample of the aggregated function — the matrix \eqn{A} in the model.
#' @param weights A matrix \eqn{L} x \eqn{I} representing the weight values associated with each sample — the matrix \eqn{y} in the model.
#' @param wavelet A string indicating the wavelet family to be used in the Discrete Wavelet Transform (DWT).
#' @param method A string specifying the shrinkage method applied to the empirical wavelet coefficients. Options are: "bayesian", "universal", "sure", "probability", or "cv".
#' @param tau A numeric value for the \eqn{\tau} parameter in the Bayesian shrinkage. If \code{NULL}, it is estimated from the data.
#' @param p A numeric value for the \eqn{p} parameter in the Bayesian shrinkage. If \code{NULL}, it is estimated from the data.
#' @param sigma A numeric value for the \eqn{\sigma} parameter in the Bayesian shrinkage. If \code{NULL}, it is estimated from the data.
#' @param MC A logical evaluating to \code{TRUE} or \code{FALSE} indicating if the integrals in the Bayesian shrinkage are approximated using Monte Carlo simulation.
#' @param type A string indicating whether the thresholding should be "soft" or "hard" (applies only when the method is not "bayesian").
#' @param singular A logical evaluating to \code{TRUE} or \code{FALSE} indicating if it adds a small constant (1e-10) to the diagonal of \eqn{yy^T} to stabilize the matrix inversion.
#' @param x A numeric vector of values at which the function is evaluated. If \code{NULL}, the default is the sequence \code{1:nrow(data)}.
#'
#' @return
#'
#' The function returns a list containing two objects:
#'
#' \describe{
#'   \item{\code{alpha}}{A matrix with the estimated functional coefficients \eqn{\alpha}.}
#'   \item{\code{Plots}}{A list of plot objects, each representing the corresponding function \eqn{\alpha_l(x)}.}
#'}
#'
#' @references dos Santos Sousa, A. R. (2024). A wavelet-based method in aggregated functional data analysis. Monte Carlo Methods and Applications, 30(1), 19-30.
#'
#' @examples
#' functional_calibration_wavelets(simulated_data$data, simulated_data$weights)
#' functional_calibration_wavelets(simulated_data$data, simulated_data$weights, tau = 5, p = 0.95, sigma = 0.1, x = simulated_data$x)
#' functional_calibration_wavelets(simulated_data$data, simulated_data$weights, method = "universal")
#'
#' @importFrom wavethresh wd threshold wr
#'
#' @export
functional_calibration_wavelets <- function(data, weights, wavelet = "DaubExPhase",
                                            method = "bayesian", tau = 1,
                                            p = NULL, sigma = NULL, MC = FALSE,
                                            type = "soft", singular = FALSE,
                                            x = NULL) {

  if (ncol(data) != ncol(weights)) {
    stop("Error: data dimension is not compatible with weights dimension.")
  }

  if (bitwAnd(nrow(data), nrow(data) - 1) != 0) {
    stop("Error: dataset length is not a power of 2.")
  }

  if (!(method %in% c("bayesian", "universal", "sure", "probability", "cv"))) {
    stop("Error: the methods of shrinkage avaiable are: bayesian, universal, sure, probability and cv.")
  }

  if (singular) {
    e <- 1e-10
  } else {
    e <- 0
  }

  if (class(data) == "data.frame") {
    data <- as.matrix(data)
  }

  if (class(weights) == "data.frame") {
    weights <- as.matrix(weights)
  }

  if (method == "bayesian") {

    D <- matrix(NA, nrow = nrow(data), ncol = ncol(data))

    for (i in 1:ncol(data)) {
      DWT <- wd(data[,i], family = wavelet)
      D[,i] <- c(DWT$D, DWT$C[length(DWT$C)])
    }

    if (is.null(sigma)) {
      sigma <- numeric(ncol(data))
      for (i in 1:ncol(data)) {
        sigma[i] <- median(abs(D[,i][1:(2^(log2(nrow(data)) - 1))]))/0.6745
      }
    } else {
      sigma <- rep(sigma, ncol(data))
    }

    if (is.null(p)) {
      p <- 1 - 1/((log2(nrow(data)) - 1):0 + 1)^2
    } else {
      p <- rep(p, log2(nrow(data)))
    }
    delta_D <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    for (i in 1:ncol(data)) {
      k <- 1
      for (j in 1:log2(nrow(data))) {
        for (l in 1:2^(log2(nrow(data)) - j)) {
          delta_D[k, i] <- Bayesian_Shrinkage(D[k, i], tau, p[j], sigma[i], MC)
          k <- k + 1 }
      }
      delta_D[k, i] <- Bayesian_Shrinkage(D[k, i], tau, p[length(p)], sigma[i], MC)
    }

  } else {
    delta_D <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
    for (i in 1:ncol(data)) {
      DWT <- wd(data[,i], family = wavelet)
      DWT <- threshold(DWT, policy = method, type = type)
      delta_D[,i] <- c(DWT$D, DWT$C[length(DWT$C)])
    }
  }
  Gama <- delta_D%*%t(weights)%*%solve(weights%*%t(weights) + e*diag(nrow(weights)))

  alpha <- matrix(NA, nrow = nrow(Gama), ncol = ncol(Gama))
  for (i in 1:ncol(Gama)) {
    DWT$D <- Gama[-nrow(Gama),i]
    DWT$C[length(DWT$C)] <- Gama[nrow(Gama),i]
    alpha[,i] <- wr(DWT)
  }

  if (is.null(x)) {
    x <- 1:nrow(alpha)
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
