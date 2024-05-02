#'AR Representation of an ARMA Model
#'
#'Output has representation with positive signs (on the right-hand side of the equation); inputs are both with positive signs (on right-hand side of equation).
#'
#'@param ar the AR-coefficient series ordered by lag.
#'@param ma the MA-coefficient series ordered by lag.
#'@param max_i the maximum index up until which to return the coefficient series.
#'
#'@details
#'Consider the ARMA model
#'\deqn{X_t = ar_1 X_{t-1} + ... + ar_p X_{t-p}+ma_1 e_{t-1}+...+ma_q e_{t-q}+e_t,}
#'where \eqn{e_t} are the innovations. \eqn{ar_i}, \eqn{i=1, ..., p}, are the AR-coefficients to pass to the
#'argument \code{ar}, \eqn{ma_j}, \eqn{j = 1, ..., q}, are the MA-coefficients
#'to pass to the argument \code{ma}.The function then returns the coefficients
#'from the corresponding infinite-order AR-representation
#'\deqn{-e_t = c_0 X_t + c_1 X_{t-1}+c_2 X_{t-2} + c_3 X_{t-3} + ...,}
#'where \eqn{c_l}, \eqn{l = 0, 1, 2, ...}, are the coefficients. Following this
#'notation, \eqn{c_0 = -1} by definition.
#'
#'@export
#'
#'@return
#'A numeric vector is returned.
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'Author
#'}
#'
#'@examples
#'arma_to_ar(ar = 0.75, ma = 0.5, max_i = 100)
#'

arma_to_ar <- function(ar = numeric(0), ma = numeric(0), max_i = 1000) {
  c(ARinftySHORT(ar = ar, ma = ma, max_i = max_i))
}

#'MA Representation of an ARMA Model
#'
#'Output has representation with positive signs (on the right-hand side of the equation); inputs are both with positive signs (on right-hand side of equation).
#'
#'@param ar the AR-coefficient series ordered by lag.
#'@param ma the MA-coefficient series ordered by lag.
#'@param max_i the maximum index up until which to return the coefficient series.
#'
#'@details
#'Consider the ARMA model
#'\deqn{X_t = ar_1 X_{t-1} + ... + ar_p X_{t-p}+ma_1 e_{t-1}+...+ma_q e_{t-q}+e_t,}
#'where \eqn{e_t} are the innovations. \eqn{ar_i}, \eqn{i=1, ..., p}, are the AR-coefficients to pass to the
#'argument \code{ar}, \eqn{ma_j}, \eqn{j = 1, ..., q}, are the MA-coefficients
#'to pass to the argument \code{ma}.The function then returns the coefficients
#'from the corresponding infinite-order MA-representation
#'\deqn{X_t = c_0 e_t + c_1 e_{t-1}+c_2 e_{t-2} + c_3 e_{t-3} + ...,}
#'where \eqn{c_l}, \eqn{l = 0, 1, 2, ...}, are the coefficients. Following this
#'notation, \eqn{c_0 = 1} by definition.
#'
#'@export
#'
#'@return
#'A numeric vector is returned.
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'Author
#'}
#'
#'@examples
#'arma_to_ma(ar = 0.75, ma = 0.5, max_i = 100)
#'

arma_to_ma <- function(ar = numeric(0), ma = numeric(0), max_i = 1000) {
  c(MAinftySHORT(ar = ar, ma = ma, max_i = max_i))
}

#'Filter Coefficients of the Fractional Differencing Operator
#'
#'Output is with positive signs on the left-hand side of the equation.
#'
#'@param d the fractional differencing coefficient.
#'@param max_i the maximum index up until which to return the coefficient series.
#'
#'@details
#'Consider the FARIMA model
#'\deqn{(1-B)^d Y_t = ar_1 X_{t-1} + ... + ar_p X_{t-p}+ma_1 e_{t-1}+...+ma_q e_{t-q}+e_t,}
#'where \eqn{e_t} are the innovations and where \eqn{X_t=(1-B)^d Y_t}.
#'\eqn{d} is the fractional differencing
#'coefficient.
#'
#'The fractional differencing operator \eqn{(1-B)^d} can alternatively be expressed
#'as an infinite coefficient series, so that
#'\deqn{(1-B)^d=\sum_{l=0}^{\infty}b_l B^k,}
#'where \eqn{B} is the backshift operator and where \eqn{b_l}, \eqn{l=0,1,2,...},
#'are the coefficients. Note that \eqn{b_0=1} by definition.
#'
#'The function returns the series of coefficients \eqn{\{b_l, l =0,1,2,...\}}.
#'
#'@export
#'
#'@return
#'A numeric vector is returned.
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'Author
#'}
#'
#'@examples
#'d_to_coef(d = 0.3, max_i = 100)
#'

d_to_coef <- function(d, max_i = 1000) {
  k <- 1:max_i
  coef_out <- c(1, (-1)^k * choose(d, k))
  coef_out
}

#'AR Representation of a FARIMA Model
#'
#'Output has representation with positive signs (on the right-hand side of the equation); inputs are both with positive signs (on right-hand side of equation).
#'
#'@param ar the AR-coefficient series ordered by lag.
#'@param ma the MA-coefficient series ordered by lag.
#'@param d the fractional differencing coefficient.
#'@param max_i the maximum index up until which to return the coefficient series.
#'
#'@details
#'Consider the FARIMA model
#'\deqn{(1-B)^d Y_t = ar_1 X_{t-1} + ... + ar_p X_{t-p}+ma_1 e_{t-1}+...+ma_q e_{t-q}+e_t,}
#'where \eqn{e_t} are the innovations and where \eqn{X_t=(1-B)^d Y_t}.
#'\eqn{ar_i}, \eqn{i=1, ..., p}, are the AR-coefficients to pass to the
#'argument \code{ar}, \eqn{ma_j}, \eqn{j = 1, ..., q}, are the MA-coefficients
#'to pass to the argument \code{ma}. \eqn{d} is the fractional differencing
#'coefficient. The function then returns the coefficients
#'from the corresponding infinite-order AR-representation
#'\deqn{-e_t = c_0 Y_t + c_1 Y_{t-1}+c_2 Y_{t-2} + c_3 Y_{t-3} + ...,}
#'where \eqn{c_l}, \eqn{l = 0, 1, 2, ...}, are the coefficients. Following this
#'notation, \eqn{c_0 = -1} by definition.
#'
#'@export
#'
#'@return
#'A numeric vector is returned.
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'Author
#'}
#'
#'@examples
#'farima_to_ar(ar = 0.75, ma = 0.5, d = 0.3, max_i = 100)
#'

farima_to_ar <- function(ar = numeric(0), ma = numeric(0), d = 0, max_i = 1000) {

  d_op <- d_to_coef(d = d, max_i = max_i)

  ar_op <- c(arma_to_ar(ar = ar, ma = ma, max_i = max_i))

  c(ARinftyLONG(ar = ar_op, d_coef = d_op))

}

#' MA Representation of a FARIMA Model
#'
#'Output has representation with positive signs (on the right-hand side of the equation); inputs are both with positive signs (on right-hand side of equation).
#'
#'@param ar the AR-coefficient series ordered by lag.
#'@param ma the MA-coefficient series ordered by lag.
#'@param d the fractional differencing coefficient.
#'@param max_i the maximum index up until which to return the coefficient series.
#'
#'@details
#'Consider the FARIMA model
#'\deqn{(1-B)^d Y_t = ar_1 X_{t-1} + ... + ar_p X_{t-p}+ma_1 e_{t-1}+...+ma_q e_{t-q}+e_t,}
#'where \eqn{e_t} are the innovations and where \eqn{X_t=(1-B)^d Y_t}.
#'\eqn{ar_i}, \eqn{i=1, ..., p}, are the AR-coefficients to pass to the
#'argument \code{ar}, \eqn{ma_j}, \eqn{j = 1, ..., q}, are the MA-coefficients
#'to pass to the argument \code{ma}. \eqn{d} is the fractional differencing coefficient.
#'The function then returns the coefficients
#'from the corresponding infinite-order AR-representation
#'\deqn{Y_t = c_0 e_t + c_1 e_{t-1}+c_2 e_{t-2} + c_3 e_{t-3} + ...,}
#'where \eqn{c_l}, \eqn{l = 0, 1, 2, ...}, are the coefficients. Following this
#'notation, \eqn{c_0 = 1} by definition.
#'
#'@export
#'
#'@return
#'A numeric vector is returned.
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'Author
#'}
#'
#'@examples
#'farima_to_ma(ar = 0.75, ma = 0.5, d = 0.3, max_i = 100)
#'

farima_to_ma <- function(ar = numeric(0), ma = numeric(0), d = 0, max_i = 1000) {

  d_op <- d_to_coef(d = d, max_i = max_i)

  ma_op <- c(arma_to_ma(ar = ar, ma = ma, max_i = max_i))

  c(MAinftyLONG(ma = ma_op, d_coef = d_op))

}
