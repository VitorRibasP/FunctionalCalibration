#' @title Simulated Data Correlated
#' @description
#' This is a simulated dataset designed to illustrate the functionalities of the package.
#' It contains 100 samples of aggregated data generated from two functions, \eqn{\alpha_1(x)} and \eqn{\alpha_2(x)}, with added noise AR(1) process with autoregressive parameter 0.5 and Gaussian process noise distributed as \eqn{N(0, 0.0075)}.
#'
#' The functions used in the simulation are:
#'
#' \deqn{ \alpha_1(x) = \sum_{j=1}^{11} h_j \left(1 + \left|\frac{x - x_j}{w_j}\right|\right)^{-4} \quad \alpha_2(x) = \sqrt{x(1-x)}\sin\left(\frac{2\pi(1 + 0.05)}{x + 0.05}\right)}
#'
#' where
#' \eqn{x = (0.1,\,0.13,\,0.15,\,0.23,\,0.25,\,0.40,\,0.44,\,0.65,\,0.76,\,0.78,\,0.81)},
#'
#' \eqn{h = (4,\,5,\,3,\,4,\,5,\,4.2,\,2.1,\,4.3,\,3.1,\,5.1,\,4.2)},
#'
#' \eqn{w = (0.005,\,0.005,\,0.006,\,0.01,\,0.01,\,0.03,\,0.01,\,0.01,\,0.005,\,0.008,\,0.005)}.
#'
#'
#' The simulations were performed over an equally spaced grid of 1024 points in the interval [0, 1].
#' These functions were linearly combined using random concentrations to generate the samples, with the addition of AR(1) noise.
#'
#' @return
#'
#' \describe{
#'   \item{\code{data}}{A data frame with 1024 rows and 100 columns. \cr
#'                          Each column represents one sample of the aggregated functions with AR(1) noise.}
#'   \item{\code{weigths}}{A data frame with 2 rows and 100 columns. \cr
#'                                Each column contains the random concentrations used to aggregate the two functions in each sample.}
#'   \item{\code{x}}{A numeric vector of length 1024. \cr
#'                   The grid of x-values used in the simulation, equally spaced from 0 to 1.}
#'   \item{\code{alphas}}{A data frame with 1024 rows and 2 columns. \cr
#'                            The true values of the functions \eqn{\alpha_1(x)} and \eqn{\alpha_2(x)} evaluated over the x grid.}
#' }
"simulated_data_cor"
#' @title Simulated Data Splines
#' @description
#' This is a simulated dataset designed to illustrate the functionalities of the package.
#' It contains 100 samples of aggregated data generated from two functions, \eqn{\alpha_1(x)} and \eqn{\alpha_2(x)}, with added Gaussian noise \eqn{N(0, 0.01)}.
#'
#' The functions used in the simulation are:
#'
#' \deqn{\alpha_1(x) = \sin(5x) e^{-x^2} \quad \alpha_2(x) = \frac{1.2 \log(1 + 9x)}{\log(10)} e^{-x}}
#'
#' The simulations were performed over an equally spaced grid of 1024 points in the interval [0, 1].
#' These functions were linearly combined using random concentrations to generate the samples, with the addition of Gaussian noise.
#'
#' @return
#'
#' \describe{
#'   \item{\code{data}}{A data frame with 1024 rows and 100 columns. \cr
#'                          Each column represents one sample of the aggregated functions with Gaussian noise.}
#'   \item{\code{weigths}}{A data frame with 2 rows and 100 columns. \cr
#'                                Each column contains the random concentrations used to aggregate the two functions in each sample.}
#'   \item{\code{x}}{A numeric vector of length 1024. \cr
#'                   The grid of x-values used in the simulation, equally spaced from 0 to 1.}
#'   \item{\code{alphas}}{A data frame with 1024 rows and 2 columns. \cr
#'                            The true values of the functions \eqn{\alpha_1(x)} and \eqn{\alpha_2(x)} evaluated over the x grid.}
#' }
"simulated_data_spl"
#' @title Simulated Data Wavelets
#' @description
#' This is a simulated dataset designed to illustrate the functionalities of the package.
#' It contains 100 samples of aggregated data generated from two functions, \eqn{\alpha_1(x)} and \eqn{\alpha_2(x)}, with added Gaussian noise \eqn{N(0, 0.01)}.
#'
#' The functions used in the simulation are:
#'
#' \deqn{ \alpha_1(x) = \sum_{j=1}^{11} h_j \left(1 + \left|\frac{x - x_j}{w_j}\right|\right)^{-4} \quad \alpha_2(x) = \sqrt{x(1-x)}\sin\left(\frac{2\pi(1 + 0.05)}{x + 0.05}\right)}
#'
#' where
#' \eqn{x = (0.1,\,0.13,\,0.15,\,0.23,\,0.25,\,0.40,\,0.44,\,0.65,\,0.76,\,0.78,\,0.81)},
#'
#' \eqn{h = (4,\,5,\,3,\,4,\,5,\,4.2,\,2.1,\,4.3,\,3.1,\,5.1,\,4.2)},
#'
#' \eqn{w = (0.005,\,0.005,\,0.006,\,0.01,\,0.01,\,0.03,\,0.01,\,0.01,\,0.005,\,0.008,\,0.005)}.
#'
#' The simulations were performed over an equally spaced grid of 1024 points in the interval [0, 1].
#' These functions were linearly combined using random concentrations to generate the samples, with the addition of Gaussian noise.
#'
#' @return
#'
#' \describe{
#'   \item{\code{data}}{A data frame with 1024 rows and 100 columns. \cr
#'                          Each column represents one sample of the aggregated functions with Gaussian noise.}
#'   \item{\code{weigths}}{A data frame with 2 rows and 100 columns. \cr
#'                                Each column contains the random concentrations used to aggregate the two functions in each sample.}
#'   \item{\code{x}}{A numeric vector of length 1024. \cr
#'                   The grid of x-values used in the simulation, equally spaced from 0 to 1.}
#'   \item{\code{alphas}}{A data frame with 1024 rows and 2 columns. \cr
#'                            The true values of the functions \eqn{\alpha_1(x)} and \eqn{\alpha_2(x)} evaluated over the x grid.}
#' }
"simulated_data_wav"
