% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.esemifar_fc.R
\name{plot.esemifar_fc}
\alias{plot.esemifar_fc}
\title{Plot Method for Class \code{"esemifar_fc"}}
\usage{
\method{plot}{esemifar_fc}(x, y = NULL, t = NULL, ...)
}
\arguments{
\item{x}{an object of class \code{"esemifar_fc"}, for example generated by a
call to \code{\link{predict.esemifar}}.}

\item{y}{currently without use; for compatibility only.}

\item{t}{a numeric vector with series of time points for observations and
forecasts.}

\item{...}{further arguments of \code{\link[stats]{plot.ts}} to adjust
for example the axis limits via \code{xlim} and \code{ylim}.}
}
\value{
This method returns \code{NULL}.
}
\description{
Create basic R plots for forecasting objects of class \code{"esemifar_fc"}.
}
\details{
This is a plot method to visualize the forecasting results for an ESEMIFAR
model. Common plot arguments can be implemented to change the
appearance.
}
\examples{
lgdp <- log(esemifar::gdpG7$gdp)
est <- tsmoothlm(lgdp, pmax = 1, qmax = 1)
# Under normality
fc <- predict(est, n.ahead = 10, method = "norm", expo = TRUE)
plot(fc)

}
\author{
\itemize{
\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
University), \cr
Author and Package Creator
}
}
