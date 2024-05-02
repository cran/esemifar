#' Data-driven Local Polynomial for the Trend's Derivatives in Equidistant Time
#' Series
#'
#'  This function runs through an iterative process in order to find the
#'  optimal bandwidth for the nonparametric estimation of the first or second
#'  derivative of the trend in an equidistant time series (with long-memory
#'  errors) and subsequently employs the obtained bandwidth via local
#'  polynomial regression.
#'
#'@param y a numeric vector that contains the time series ordered from past to
#'present.
#'@param d an integer \code{1} or \code{2} that defines the order of
#'derivative; the default is \code{d = 1}.
#'@param pmin an integer value \eqn{>= 0} that defines the minimum
#'autoregressive order to calculate the BIC-criterion for; is set to \code{0}
#'by default; decimal numbers will be rounded off to integers.
#'@param pmax an integer value \eqn{>= 0} that defines the maximum
#'autoregressive order to calculate the BIC-criterion for; is set to \code{0}
#'by default; decimal numbers will be rounded off to integers.
#'@param qmin an integer value \eqn{>= 0} that defines the minimum
#'moving-average order to calculate the BIC-criterion for; is set to \code{0}
#'by default; decimal numbers will be rounded off to integers.
#'@param qmax an integer value \eqn{>= 0} that defines the maximum
#'moving-average order to calculate the BIC-criterion for; is set to \code{0}
#'by default; decimal numbers will be rounded off to integers.
#'@param mu an integer \code{0}, ..., \code{3} that represents the smoothness
#'parameter of the kernel weighting function and thus defines the kernel
#'function that will be used within the local polynomial regression; is set to
#'\code{1} by default.
#'@param mu.p an integer \code{0}, ..., \code{3} that represents the smoothness
#'parameter of the kernel weighting function for the iterative process to
#'obtain initial estimates for \eqn{c_f}, \eqn{d} and \eqn{b_0}; is set to
#'\code{1} by default.
#'
#'\tabular{cl}{
#'\strong{Number} \tab \strong{Kernel}\cr
#'\code{0} \tab Uniform Kernel\cr
#'\code{1} \tab Epanechnikov Kernel\cr
#'\code{2} \tab Bisquare Kernel\cr
#'\code{3} \tab Triweight Kernel
#'}
#'@param pp an integer \code{1} (local linear regression) or \code{3} (local
#'cubic regression) that indicates the order of polynomial upon which
#'\eqn{c_f}, \eqn{d} and \eqn{b_0} will be calculated by
#'\code{\link{tsmoothlm}}; the default is \code{pp = 1}.
#'@param bStart.p a numeric object that indicates the starting value of the
#'bandwidth for the iterative process to obtain initial estimates for \eqn{c_f},
#'\eqn{d} and \eqn{b_0}; should be \eqn{> 0}; is set to \code{0.15} by default.
#'@param InfR.p a character object that represents the inflation
#'rate in the form \eqn{h_d = h^a} of the bandwidth for the iterative process
#'to obtain initial estimates for \eqn{c_f},
#'\eqn{d} and \eqn{b_0}; is set to \code{"Opt"} by default.
#'
#'@details
#'
#'The trend is estimated based on the additive
#'nonparametric regression model for an equidistant time series
#'\deqn{y_t = m(x_t) + \epsilon_t,}
#'where \eqn{y_t} is the observed time series, \eqn{x_t} is the rescaled time
#'on the interval \eqn{[0, 1]}, \eqn{m(x_t)} is a smooth and deterministic
#'trend function and \eqn{\epsilon_t} are stationary errors with
#'\eqn{E(\epsilon_t) = 0} and is assumed to follow a FARIMA(\eqn{p, d, q})
#'model (see also Beran and Feng, 2002).
#'
#'The iterative-plug-in (IPI) algorithm, which numerically minimizes the
#'Asymptotic Mean Squared Error (AMISE), is based on the proposal of Beran
#'and Feng (2002a).
#'
#'The variance factor \eqn{c_f}, the long memory parameter \eqn{d} and the
#'starting bandwidth \eqn{b0} are first obtained from a pilot-estimation of
#'the time series' nonparametric trend (\eqn{\nu = 0}) with polynomial order
#'\eqn{p_p}. The estimate is then plugged into the iterative procedure for
#'estimating the first or second derivative (\eqn{\nu = 1} or \eqn{\nu = 2}).
#'For further details on the asymptotic theory or the algorithm, we refer the
#'user to Letmathe, Beran and Feng (2021).
#'
#'The function itself is applicable in the following way: Based on a data input
#'\code{y}, an order of polynomial \code{pp} for the variance factor estimation
#'procedure, a starting value for the relative bandwidth \code{bStart.p} in the
#'variance factor estimation procedure and a kernel function defined by the
#'smoothness parameter \code{mu}, an optimal bandwidth is numerically calculated
#'for the trend's derivative of order \code{d}. In fact, aside from the input
#'vector \code{y}, every argument has a default setting that can be adjusted for
#'the individual case. However, it is recommended to initially use the default
#'values for the estimation of the
#'first derivative and adjust the argument \code{d} to \code{d = 2} for the
#'estimation of the second derivative.
#'The initial bandwidth does not affect the resulting optimal bandwidth in
#'theory. However in practice, local minima of the AMISE can influence the
#'results. For more specific information on the input arguments consult the
#'section \emph{Arguments}.
#'
#'After the bandwidth estimation, the nonparametric derivative of the series
#'is calculated with respect to the obtained optimal bandwidth by means of a
#'local polynomial regression. The output object is then a list that contains,
#'among other components, the original time series, the estimates of the
#'derivative and the estimated optimal bandwidth.
#'
#'The default print method for this function delivers key numbers such as
#'the iteration steps and the generated optimal bandwidth rounded to the fourth
#'decimal. The exact numbers and results such as the estimated nonparametric
#'trend series are saved within the output object and can be addressed via the
#'\code{$} sign.
#'
#'
#'@return The function returns a list with different components:
#'
#'\describe{
#'\item{b0}{the optimal bandwidth chosen by the IPI-algorithm.}
#'\item{bStart.p}{the starting bandwidth for the nonparametric trend estimation
#'that leads to the initial estimates; input argument.}
#'\item{cf0}{the estimated variance factor.}
#'\item{InfR.p}{the inflation rate setting.}
#'\item{iterations}{the bandwidths of the single iterations steps}
#'\item{mu.p}{the smoothness parameter of the second order kernel; input
#'argument.}
#'\item{n}{the number of observations.}
#'\item{niterations}{the total number of iterations until convergence.}
#'\item{orig}{the original input series; input argument.}
#'\item{p}{the order of polynomial for the local polynomial
#'regression used within derivative estimation procedure.}
#'\item{pp}{the order of polynomial for the local polynomial
#'regression used in the pilot estimation; input argument.}
#'\item{v}{the considered order of the trend's derivative; input argument
#'\code{d}.}
#'\item{ws}{the weighting system matrix used within the local polynomial
#'regression; this matrix is a condensed version of a complete weighting system
#'matrix; in each row of \code{ws}, the weights for conducting the smoothing
#'procedure at a specific observation time point can be found; the first
#'\eqn{[nb + 0.5]} rows, where \eqn{n} corresponds to the number of
#'observations, \eqn{b} is the bandwidth considered for smoothing and
#'\eqn{[.]} denotes the integer part, contain the weights at the
#'\eqn{[nb + 0.5]} left-hand boundary points; the weights in row
#'\eqn{[nb + 0.5] + 1} are representative for the estimation at all
#'interior points and the remaining rows contain the weights for the right-hand
#'boundary points; each row has exactly \eqn{2[nb + 0.5] + 1} elements,
#'more specifically the weights for observations of the nearest
#'\eqn{2[nb + 0.5] + 1} time points; moreover, the weights are normalized,
#'i.e. the weights are obtained under consideration of the time points
#'\eqn{x_t = t/n}, where \eqn{t = 1, 2, ..., n}.}
#'\item{ye}{the nonparametric estimates of the derivative.}
#'}
#'
#'@md
#'
#'@export
#'
#'@references
#' Letmathe, S., Beran, J. and Feng, Y. (2021). An extended exponential SEMIFAR
#' model with application in R. Discussion Paper. Paderborn University.
#'
#'
#'@author
#'\itemize{
#'\item Yuanhua Feng (Department of Economics, Paderborn University), \cr
#'Author of the Algorithms \cr
#'Website: \url{https://wiwi.uni-paderborn.de/en/dep4/feng/}
#'\item Sebastian Letmathe (Scientific Employee) (Department of Economics,
#'Paderborn
#'University), \cr
#'Package Creator and Maintainer
#'\item Dominik Schulz (Scientific Employee) (Department of Economics, Paderborn
#'University), \cr
#'Author
#'}
#'
#'@examples
#'\donttest{
#'# Logarithm of test data
#'test_data <- gdpG7
#'y <- log(test_data$gdp)
#'n <- length(y)
#'t <- seq(from = 1962, to = 2020, length.out = n)
#'
#'# Applied dsmooth function for the trend's first derivative
#'result_d <- dsmoothlm(y, d = 1, pp = 1, pmax = 1, qmax = 1, InfR.p = "Opt")
#'estim <- result_d$ye
#'
#'# Plot of the results
#'plot(t, estim, xlab = "Year", ylab = "First derivative", type = "l",
#'  main = paste0("Estimated first derivative of the trend for log-quarterly ",
#'  "G7-GDP, Q1 1962 - Q4 2019"), cex.axis = 0.8, cex.main = 0.8,
#'  cex.lab = 0.8, bty = "n")
#'
#'# Print result
#'result_d
#'}
#'# The main function "dsmoothlm"------------------------------------------

# d = 1 or 2 indicate iterative plug-in for m' or m'', respectively
dsmoothlm <- function(y,
                      d = c(1, 2),
                      pmin = c(0, 1, 2, 3, 4, 5),
                      pmax = c(0, 1, 2, 3, 4, 5),
                      qmin = c(0, 1, 2, 3, 4, 5),
                      qmax = c(0, 1, 2, 3, 4, 5),
                      mu = c(0, 1, 2, 3),
                      mu.p = c(0, 1, 2, 3),
                      pp = c(1, 3),
                      bStart.p = 0.15,
                      InfR.p = c("Opt", "Nai", "Var")) {

  if (length(y) <= 1 || !all(!is.na(y)) || !is.numeric(y)) {
    stop("The argument 'y' must be a numeric vector with length > 1 and ",
         "without NAs.")
  }

  if (!(length(d) %in% c(1, 2)) || !all(!is.na(d)) || !is.numeric(d)) {
    stop("The argument 'd' must be a single integer value (1 or 2).")
  }
  d <- floor(d)

  if (!(length(mu) %in% c(1, 4)) || !all(!is.na(mu)) || !is.numeric(mu)) {
    stop("The argument 'mu' must be a single integer value (0, 1, 2 or 3).")
  }
  mu <- floor(mu)

  if (!(length(mu.p) %in% c(1, 4)) || !all(!is.na(mu.p)) || !is.numeric(mu.p)) {
    stop("The argument 'mu' must be a single integer value (0, 1, 2 or 3).")
  }
  mu.p <- floor(mu.p)

  if (!(length(pp) %in% c(1, 2)) || !all(!is.na(pp)) || !is.numeric(pp)) {
    stop("The argument 'pp' must be a single integer value (either 1 or 3).")
  }
  pp <- floor(pp)

  if (length(InfR.p) != 1 || !(InfR.p %in% c("Opt", "Nai", "Var"))) {
    stop("Input of argument 'InfR.p' incorrect. Input not recognized.")
  }
  if (length(bStart.p) != 1 || is.na(bStart.p) || !is.numeric(bStart.p) ||
      (bStart.p <= 0)) {
    stop("The argument 'bStart.p' must be a single non-NA double value with ",
         "bStart.p > 0.")
  }

  # if (!length(method_error) %in% c(1, 2) || !all(!is.na(method_error)) ||
  #     !is.character(method_error)) {
  #   stop("The argument 'method_error' must be a single non-NA character value.")
  # }

  if (all(d == c(1, 2))) d <- 1
  if (all(mu == c(0, 1, 2, 3))) mu <- 1
  if (all(mu.p == c(0, 1, 2, 3))) mu.p <- 1
  if (all(pp == c(1, 3))) pp <- 1
  if (all(pmin == c(0, 1, 2, 3, 4, 5))) pmin <- 0
  if (all(pmax == c(0, 1, 2, 3, 4, 5))) pmax <- 0
  if (all(qmin == c(0, 1, 2, 3, 4, 5))) qmin <- 0
  if (all(qmax == c(0, 1, 2, 3, 4, 5))) pmax <- 0
  if (all(InfR.p == c("Opt", "Nai", "Var"))) InfR.p <- "Opt"
  # if (all(method_error == c("fracdiff", "arfima"))) method_error <- "fracdiff"


  if (!(d %in% c(1, 2))) {
    stop("Input of argument 'd' incorrect. It must be set to either 1 or 2.")
  }
  if (!length(pmin) %in% c(1, 5) || !all(!is.na(pmin)) || !is.numeric(pmin)) {
    stop("The argument 'pmin' must be a single integer value (1, 2, 3, 4 or 5).")
  }
  pmin <- floor(pmin)

  if (!length(pmax) %in% c(1, 5) || !all(!is.na(pmax)) || !is.numeric(pmax)) {
    stop("The argument 'pmax' must be a single integer value (1, 2, 3, 4 or 5).")
  }
  pmax <- floor(pmax)

  if (!length(qmin) %in% c(1, 5) || !all(!is.na(qmin)) || !is.numeric(qmin)) {
    stop("The argument 'pmin' must be a single integer value (1, 2, 3, 4 or 5).")
  }
  qmin <- floor(qmin)

  if (!length(qmax) %in% c(1, 5) || !all(!is.na(qmax)) || !is.numeric(qmax)) {
    stop("The argument 'pmin' must be a single integer value (1, 2, 3, 4 or 5).")
  }
  qmax <- floor(qmax)

  if (!(pp %in% c(1, 3))) {
    stop("Input of argument 'pp' incorrect. It must be set to either 1 or 3.")
  }

  # if (length(method_error) != 1 || !(method_error %in% c("fracdiff", "arfima"))) {
  #   stop("Input of argument 'method_error' incorrect. Method not recognized.")
  # }

  result.p <- tsmoothlmCalc(y, qmin = qmin, qmax = qmax, pmin = pmin,
                            pmax = pmax, p = pp, mu = mu.p, InfR = InfR.p,
                            bStart.p, bb = 1, method = "lpr")

  cf0 <- result.p$cf0
  d.BIC <- result.p$d.BIC
  p.BIC <- result.p$p.BIC
  q.BIC <- result.p$q.BIC
  bStart <- result.p$b0
  # Using the knn idea, bb=1, or not
  bb <- 1  # default
  if (d == 1) {InfR <- "Nai"} else {InfR <- "Var"}  # default
  cb <- 0.05

  # Input parameters
  n <- length(y)
  p <- d + 1
  k <- p + 1
  pd <- p + 2
  runc <- 1
  n1 <- trunc(n * cb)

  # New method for the kernel constants with p = 1 or 3------------------------

  # Kernel constants-----------------------------------------------------------
  m <- 1000000  # for the numerical integral
  u <- (-m:m) / (m + 0.5)
  # For d = 1, the four 3rd-order kernels for m' (Table 5.7, Mueller, 1988)
  Vc <- 2 * gamma(1 - 2 * d.BIC) * sin(pi * d.BIC) # Constant in the variance
  wkp <- lookup$lookup_3kerns[(mu + 1), d][[1]](u)

  Rp <- Vc * lookup$lookup_3kerns_kdf[(mu + 1), d][[1]](d.BIC)
  mukp <- sum((u ** k) * wkp) / m

  # Two constants in the bandwidth
  c1 <- factorial(k) ** 2 * (2 * d + 1) / (2 * (k - d))  # for d = 1 or 2
  c2 <- (1 - 2 * cb) * (1 - 2 * d.BIC) * Rp / (mukp) ** 2

  steps <- rep(NA, 40)

  bd2_func <- lookup$InfR2_lookup[d, InfR][[1]]

  # The main iteration---------------------------------------------------------

  noi <- 40  # maximal number of iterations

  for (i in 1:noi) {
    if (runc == 1) {
      if (i > 1) {bold1 <- bold}
      # Look up EIM inflation rates from the lookup table
      if (i == 1) {
        bold <- bStart ; bd <- bStart
      } else {
          bold <- bopt ; bd <- bd2_func(bold, d.BIC)
          }
      if (bd >= 0.49) {bd <- 0.49}

      yed <- smoots::gsmooth(y, k, pd, mu, bd, bb)$ye
      I2 <- sum(yed[max(1, n1):((1 - cb) * n)]^2) / (n - 2 * n1)
	    c3 <- cf0 / I2

      if (d == 1) {
        bopt <- (c1 * c2 * c3) ** (1 / (7 - 2 * d.BIC)) * n ** ((2 * d.BIC - 1)
                                                                / (7 - 2 * d.BIC))
        if (bopt < n ** (-7 / 9)) {bopt <- n ** (-7 / 9)}
      }
      if (d == 2) {
        bopt <- (c1 * c2 * c3) ** (1 / (9 - 2 * d.BIC)) * n ** ((2 * d.BIC - 1)
                                                                / (9 - 2 * d.BIC))
        if (bopt < n ** (-9 / 11)) {bopt <- n ** (-9 / 11)}
      }
      if (bopt > 0.49) {bopt <- 0.49}
      steps[i] <- bopt

      if (i > 2 && abs(bold - bopt) / bopt < 1 / n) {runc <- 0}
      if (i > 3 && abs(bold1 - bopt) / bopt < 1 / n){
        bopt <- (bold + bopt) / 2
        runc <- 0
      }
    }
  }

  # Smooth with the selected bandwidth-----------------------------------------

  if (d == 1 && bopt < n ** (-7 / 9)) {bopt <- n ** (-7 / 9)}
  if (d == 2 && bopt < n ** (-9 / 11)) {bopt <- n ** (-9 / 11)}
  if (bopt > 0.49) {bopt <- 0.49}

  est.opt <- smoots::gsmooth(y, d, p, mu, bopt, bb)
  ye <- est.opt$ye
  ws <- est.opt$ws
  # Final results
  results <- list(ye = ye, b0 = bopt, cf0 = cf0, p.BIC = p.BIC, q.BIC = q.BIC,
       d.BIC = d.BIC, v = d, n = n, niterations = length(steps[!is.na(steps)]),
       orig = y, iterations = steps[!is.na(steps)], pp = pp,
       mu.p = mu.p, InfR.p = InfR.p, bStart = bStart, bb = bb, cb = cb,
       ws = ws, p = d + 1)
  class(results) <- "esemifar"
  attr(results, "function") <- "dsmoothlm"

  drop(results)
}
# End of the function
