#' esemifar: A package for data-driven nonparametric estimation of the trend and
#' its derivatives in equidistant time series.
#'
#' The \code{esemifar} package provides different applicable functions for the
#' estimation of the trend or its derivatives in equidistant time series.
#' The main functions include an automated bandwidth selection method for time
#' series with long-memory errors.
#'
#'@section Functions (version 1.0.0):
#'The \code{esemifar} functions are either meant for calculating nonparametric
#'estimates of the trend of a time series or its derivatives.
#'
#'
#'\code{dsmoothlm} is a function that calculates the derivatives of the
#'trend after obtaining the optimal bandwidth by an iterative plug-in
#'algorithm.
#'
#'\code{tsmoothlm} is the central function of the package. It allows
#'the user to conduct a local polynomial regression of the trend based on
#'an optimal bandwidth that is obtained by an iterative plug-in algorithm.
#'Inflation rate (and other factors) can be manually
#'and individually adjusted as arguments in the function
#'(see also: \code{\link{tsmoothlm}}).
#'
#'\code{critMatlm} is a quick tool for the calculation of information criteria
#'for FARIMA(\eqn{p,d,q}) models with different order combinations \eqn{p} and
#'\eqn{q}. The function returns a matrix with the obtained values of the
#'selected criterion for the different combinations of \eqn{p} and \eqn{q}
#'(see also: \code{\link{critMatlm}}).
#'

#'@section Datasets:
#'The package includes two datasets: \code{airLDN} (see also:
#'\code{\link{airLDN}}) with daily observations of individual air pollutants
#'from 2014 to 2020 and \code{gdpG7} (see also: \code{\link{gdpG7}}) that has
#'data concerning the quarterly G7 GDP between Q1 1962 and Q4 2019.
#'
#'@section License:
#'The package is distributed under the General Public License v3
#'([GPL-3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))).
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
#' Letmathe, S., Beran, J. and Feng, Y. (2021). An extended exponential SEMIFAR
#' model with application in R. Discussion Paper. Paderborn University.
#'
#'@author
#'\itemize{
#'\item Yuanhua Feng (Department of Economics, Paderborn University), \cr
#'Author of the Algorithms \cr
#'Website: \url{https://wiwi.uni-paderborn.de/en/dep4/feng/}
#'\item Sebastian Letmathe (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'Package Creator and Maintainer
#'}
#'
#' @docType package
#' @name esemifar
NULL
