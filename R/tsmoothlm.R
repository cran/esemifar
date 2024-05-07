#' Advanced Data-driven Nonparametric Regression for the Trend in Equidistant
#' Time Series
#'
#' This function runs an iterative plug-in algorithm to find the optimal
#' bandwidth for the estimation of the nonparametric trend in equidistant
#' time series (with long-memory errors) and then employs the resulting
#' bandwidth via either local polynomial or kernel regression.
#'
#'@param y a numeric vector that contains the time series ordered from past to
#'present.
#'@param p an integer \code{1} (local linear regression) or \code{3} (local
#'cubic regression); represents the order of polynomial within the local
#'polynomial regression (see also the 'Details' section); is set to \code{1} by
#'default; is automatically set to \code{1} if \code{method = "kr"}.
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
#'
#'\tabular{cl}{
#'\strong{Number} \tab \strong{Kernel}\cr
#'\code{0} \tab Uniform Kernel\cr
#'\code{1} \tab Epanechnikov Kernel\cr
#'\code{2} \tab Bisquare Kernel\cr
#'\code{3} \tab Triweight Kernel
#'}
#'@param InfR a character object that represents the inflation
#'rate in the form \eqn{h_d = h^a} for the bandwidth in the estimation of
#'\eqn{I[m^{(k)}]}{I[m^(k)]} (see also the 'Details' section); is set to
#'\code{"Opt"} by default.
#'
#'\tabular{cl}{
#'\strong{Inflation rate} \tab \strong{Description}\cr
#'\code{"Opt"} \tab Optimal inflation rate \eqn{a_{p,O}}{a_[p,O]}
#'(\eqn{(5-2d)/(7-2d)} for \code{p = 1}; \eqn{(9-2d)/(11-2d)} for
#'\code{p = 3})\cr
#'\code{"Nai"} \tab Naive inflation rate \eqn{a_{p,N}}{a_[p,N]}
#'(\eqn{(5-2d)/(9-2d)} for \code{p = 1}; \eqn{(9-2d)/(13-2d)} for
#'\code{p = 3})\cr
#'\code{"Var"} \tab Stable inflation rate \eqn{a_{p,S}}{a_[p,S]} (\eqn{1/2} for
#'\code{p = 1} and \code{p = 3})
#'}
#'@param bStart a numeric object that indicates the starting value of the
#'bandwidth for the iterative process; should be \eqn{> 0}; is set to
#'\code{0.15} by default.
#'@param bb can be set to \code{0} or \code{1}; the parameter controlling the
#'bandwidth used at the boundary; is set to \code{1} by default.
#'
#'\tabular{cl}{
#'\strong{Number (\code{bb})} \tab \strong{Estimation procedure at boundary
#'points}\cr
#'\code{0} \tab Fixed bandwidth on one side with possible large
#'bandwidth on the other side at the boundary\cr
#'\code{1} \tab The \eqn{k}-nearest neighbor method will be used
#'}
#'@param cb a numeric value that indicates the percentage of omitted
#'observations on each side of the observation period for the automated
#'bandwidth selection; is set to \code{0.05} by default.
#'@param method the final smoothing approach; \code{"lpr"} represents the local
#'polynomial regression, whereas \code{"kr"} implements a kernel regression;
#'is set to \code{"lpr"} by default.
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
#'model (see also Beran and Feng, 2002a, Beran and Feng, 2002b and Beran
#'and Feng, 2002c).
#'
#'The iterative-plug-in (IPI) algorithm, which numerically minimizes the
#'Asymptotic Mean Squared Error (AMISE), is based on the proposal of Beran
#'and Feng (2002a).
#'
#'
#'The function calculates suitable estimates for \eqn{c_f}, the variance
#'factor, and \eqn{I[m^{(k)}]}{I[m^(k)]} over different iterations. In each
#'iteration, a bandwidth is obtained in accordance with the AMISE that once
#'more serves as an input for the following iteration. The process repeats
#'until either convergence or the 40th iteration is reached. For further
#'details on the asymptotic theory or the algorithm, please see Letmathe et
#'al., 2023.
#'
#'To apply the function, the following arguments are needed: a data input
#'\code{y}, an order of polynomial \code{p}, a kernel weighting function
#'defined by the smoothness parameter \code{mu}, an inflation rate setting
#'\code{InfR} (see also Beran and Feng, 2002b), a starting value for the
#'relative bandwidth \code{bStart}, a
#'boundary method \code{bb}, a boundary cut-off percentage \code{cb} and a
#'final smoothing method \code{method}. In fact, aside from the input vector
#'\code{y}, every argument has a default setting that can be adjusted for the
#'individual case. Theoretically, the initial bandwidth does not affect the
#'selected optimal bandwidth. However, in practice local minima of the AMISE
#'might exist and influence the selected bandwidth. Therefore, the default
#'setting is \code{bStart = 0.15}. In the rare
#'case of a clearly unsuitable optimal bandwidth, a starting bandwidth that
#'differs from the default value is a first possible approach to obtain a
#'better result. Other argument adjustments can be tried as well. For more
#'specific information on the input arguments consult the section
#'\emph{Arguments}.
#'
#'When applying the function, an optimal bandwidth is obtained based on a
#'strongly modified version of the IPI algorithm of Beran and Feng (2002a). In
#'a second step, the nonparametric trend of the series is calculated with
#'respect to the chosen bandwidth and the selected regression method (\code{lpf}
#'or \code{kr}). Please note that \code{method = "lpf"} is strongly recommended
#'by the authors. Moreover, it is notable that \code{p} is automatically set to
#'\code{1} for \code{method = "kr"}. The output object is then a list that
#'contains, among other components, the original time series, the estimated
#'trend values and the series without the trend.
#'
#'The default print method for this function delivers only key numbers such as
#'the iteration steps and the generated optimal bandwidth rounded to the fourth
#'decimal. The exact numbers and results such as the estimated nonparametric
#'trend series are saved within the output object and can be addressed via the
#'\code{$} sign.
#'
#'
#'@return The function returns a list with different components:
#'
#'\describe{
#'\item{FARIMA.BIC}{the Bayesian Information Criterion of the optimal
#'FARIMA(\eqn{p,d,q}) model.}
#'\item{cb}{the percentage of omitted observations on each side of the
#'observation period.}
#'\item{b0}{the optimal bandwidth chosen by the IPI-algorithm.}
#'\item{bb}{the boundary bandwidth method used within the IPI; always equal to
#'1.}
#'\item{bStart}{the starting value of the (relative) bandwidth; input
#'argument.}
#'\item{cf0}{the estimated variance factor; in contrast to the definitions
#'given in the \emph{Details} section, this object actually contains an
#'estimated value of \eqn{2\pi c_f}, i.e. it corresponds to the estimated sum
#'of autocovariances.}
#'\item{d.BIC}{the long-memory parameter of the optimal FARIMA(\eqn{p,d,q})
#'model.}
#'\item{FARMA.BIC}{the model fit of the selected FARIMA(\eqn{p,d,q} model.}
#'\item{I2}{the estimated value of \eqn{I[m^{(k)}]}{I[m^(k)]}.}
#'\item{InfR}{the setting for the inflation rate according to the chosen
#'algorithm.}
#'\item{iterations}{the bandwidths of the single iterations steps}
#'\item{mu}{the smoothness parameter of the second order kernel; input
#'argument.}
#'\item{n}{the number of observations.}
#'\item{niterations}{the total number of iterations until convergence.}
#'\item{orig}{the original input series; input argument.}
#'\item{p.BIC}{the order p of the optimal FARIMA(\eqn{p,d,q}) model.}
#'\item{p}{the order of polynomial used in the IPI-algorithm; also used for the
#'final smoothing, if \code{method = "lpr"}; input argument.}
#'\item{q.BIC}{the order \eqn{q} of the optimal FARIMA(\eqn{p,d,q})
#'model.}
#'\item{res}{the estimated residual series.}
#'\item{v}{the considered order of derivative of the trend; is always zero for
#'this function.}
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
#'\item{ye}{the nonparametric estimates of the trend.}
#'}
#'
#'@export
#'
#'@references
#' Beran, J. and Y. Feng (2002a). Iterative plug-in algorithms for SEMIFAR
#' models - definition, convergence, and asymptotic properties. Journal of
#' Computational and Graphical Statistics 11(3), 690-713.
#'
#' Beran, J. and Feng, Y. (2002b). Local polynomial fitting with long-memory,
#' short-memory and antipersistent errors. Annals of the Institute of
#' Statistical Mathematics, 54(2), 291-311.
#'
#' Beran, J. and Feng, Y. (2002c). SEMIFAR models - a semiparametric approach
#' to modelling trends, longrange dependence and nonstationarity. Computational
#' Statistics & Data Analysis 40(2), 393-419.
#'
#' Letmathe, S., Beran, J. and Feng, Y. (2023). An extended exponential SEMIFAR
#' model with application in R. Communications in Statistics - Theory and Methods:
#' 1-13.
#'
#'
#'@author
#'\itemize{
#'\item Yuanhua Feng (Department of Economics, Paderborn University), \cr
#'Author of the Algorithms \cr
#'Website: \url{https://wiwi.uni-paderborn.de/en/dep4/feng/}
#'\item Sebastian Letmathe (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'Package Creator and Maintainer
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'Author
#'}
#'
#'@importFrom Rcpp cppFunction
#'@useDynLib esemifar
#'
#'@examples
#'
#'\donttest{
#'### Example 1: G7-GDP ###
#'
#'# Logarithm of test data
#'# -> the logarithm of the data is assumed to follow the additive model
#'test_data <- gdpG7
#'y <- log(test_data$gdp)
#'n <- length(y)
#'
#'# Applied tsmooth function for the trend
#'result <- tsmoothlm(y, p = 1, pmax = 1, qmax = 1, InfR = "Opt")
#'trend1 <- result$ye
#'
#'# Plot of the results
#'t <- seq(from = 1962, to = 2020, length.out = n)
#'plot(t, y, type = "l", xlab = "Year", ylab = "log(G7-GDP)", bty = "n",
#'  lwd = 1, lty = 3,
#'  main = "Estimated trend for log-quarterly G7-GDP, Q1 1962 - Q4 2019")
#'points(t, trend1, type = "l", col = "red", lwd = 1)
#'title(sub = expression(italic("Figure 1")), col.sub = "gray47",
#'  cex.sub = 0.6, adj = 0)
#'result
#'}
#'
#'
# The main function------------------------------------------------------------

tsmoothlm <- function(y,
                      pmin = c(0, 1, 2, 3, 4, 5),
                      pmax = c(0, 1, 2, 3, 4, 5),
                      qmin = c(0, 1, 2, 3, 4, 5),
                      qmax = c(0, 1, 2, 3, 4, 5),
                      p = c(1, 3),
                      mu = c(0, 1, 2, 3),
                      InfR = c("Opt", "Nai", "Var"),
                      bStart = 0.15,
                      bb = c(0, 1),
                      cb = 0.05,
                      method = c("lpr", "kr"))  {

  # Input if no inputs were made to specific arguments

  if (length(y) <= 1 || !all(!is.na(y)) || !is.numeric(y)) {
    stop("The argument 'y' must be a numeric vector with length > 1 and ",
         "without NAs.")
  }

  if (!length(p) %in% c(1, 2) || !all(!is.na(p)) || !is.numeric(p)) {
    stop("The argument 'p' must be a single integer value (either 1 or 3).")
  }
  p <- floor(p)

  if (!length(mu) %in% c(1, 4) || !all(!is.na(mu)) || !is.numeric(mu)) {
    stop("The argument 'mu' must be a single integer value (0, 1, 2 or 3).")
  }
  mu <- floor(mu)

  if (!length(InfR) %in% c(1, 3)) {
    stop("The argument 'InfR' must be a single non-NA character value.")
  }

  if (length(bStart) != 1 || is.na(bStart) || !is.numeric(bStart) ||
      (bStart <= 0)) {
    stop("The argument 'bStart' must be a single non-NA double value with ",
         "bStart > 0.")
  }


  if (!length(bb) %in% c(1, 2) || !all(!is.na(bb)) || !is.numeric(bb)) {
    stop("The argument 'bb' must be a single non-NA integer (either 0 or 1).")
  }
  bb <- floor(bb)

  if (length(cb) != 1 || is.na(cb) || !is.numeric(cb) ||
      (cb < 0 || cb >= 0.5)) {
    stop("The argument 'cb' must be a single non-NA double value with ",
         "0 <= cb < 0.5.")
  }

  if (!length(method) %in% c(1, 2) || !all(!is.na(method)) ||
      !is.character(method)) {
    stop("The argument 'method' must be a single non-NA character value.")
  }

  # if (!length(method_error) %in% c(1, 2) || !all(!is.na(method_error)) ||
  #     !is.character(method_error)) {
  #   stop("The argument 'method_error' must be a single non-NA character value.")
  # }


  if (all(mu == c(0, 1, 2, 3))) mu <- 1
  if (all(InfR == c("Opt", "Nai", "Var"))) InfR <- "Opt"
  if (all(bb == c(0, 1))) bb <- 1
  if (all(method == c("lpr", "kr"))) method <- "lpr"
  # if (all(method_error == c("fracdiff", "arfima"))) method_error <- "fracdiff"
  if (all(p == c(1, 3)) || method == "kr") p <- 1
  if (all(pmin == c(0, 1, 2, 3, 4, 5))) pmin <- 0
  if (all(pmax == c(0, 1, 2, 3, 4, 5))) pmax <- 0
  if (all(qmin == c(0, 1, 2, 3, 4, 5))) qmin <- 0
  if (all(qmax == c(0, 1, 2, 3, 4, 5))) qmax <- 0

  # Stop, if incorrect inputs were given
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
  if (length(p) != 1 || !(p %in% c(1, 3))) {
    stop("Input of argument 'p' incorrect. It must be either 1 or 3.")
  }
  if (length(mu) != 1 || !(mu %in% 0:4)) {
    stop("Input of argument 'mu' incorrect. It must be an integer from 0 to 4.")
  }
  if (length(InfR) != 1 || !(InfR %in% c("Opt", "Nai", "Var"))) {
    stop("Input of argument 'InfR' incorrect. Input not recognized.")
  }
  if(length(bb) != 1 || !(bb %in% c(0, 1))) {
    stop("Input of argument 'bb' incorrect. It must be either 0 or 1.")
  }
  if (length(method) != 1 || !(method %in% c("lpr", "kr"))) {
    stop("Input of argument 'method' incorrect. Method not recognized.")
  }
  # if (length(method_error) != 1 || !(method_error %in% c("fracdiff", "arfima"))) {
  #   stop("Input of argument 'method_error' incorrect. Method not recognized.")
  # }

  results <- tsmoothlmCalc(y = y, p = p, pmin = pmin, pmax = pmax, qmin = qmin,
                           qmax = qmax, mu = mu, InfR = InfR, bStart = bStart,
                           bb = bb, cb = cb, method = method)

  class(results) <- "esemifar"
  attr(results, "function") <- "tsmoothlm"
  attr(results, "type") <- "tsmoothlm"
  results
}
# End of the function
