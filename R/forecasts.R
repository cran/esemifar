
#'Create an ESEMIFAR Model Object Suitable for Forecasting
#'
#'An alternative method to create an ESEMIFAR estimation object
#'stitched together from a nonparametric and a parametric part.
#'
#'@param nonpar_model an estimation object returned by \code{\link{gsmooth}}.
#'@param par_model an estimation object returned by \code{\link[fracdiff]{fracdiff}} fitted
#'to the residuals of \code{nonpar_model}.
#'
#'@details
#'The main function \code{\link{tsmoothlm}} already returns a fully estimated
#'ESEMIFAR model. In some instances, alternative specifications of the
#'nonparametric and parametric model parts are needed, for which
#'\code{\link{tsmoothlm}} with its automated estimation algorithm does not
#'provide sufficient flexibility. Therefore, this function allows to stitch
#'together a nonparametric model part returned by \code{\link{gsmooth}} and
#'a FARIMA part for the residuals obtained via \code{\link[fracdiff]{fracdiff}}.
#'The resulting object can then be used for forecasting.
#'
#'@importFrom stats qnorm
#'@importFrom utils tail
#'
#'@export
#'
#'@return
#'The function returns a list of class \code{"esemifar"} with elements
#'\code{nonpar_model} and \code{par_model}.
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'Author
#'}
#'
#'@examples
#'lgdp <- log(esemifar::gdpG7$gdp)
#'nonpar <- gsmooth(lgdp, b = 0.15)
#'res <- nonpar$res
#'par <- fracdiff::fracdiff(res, nar = 1, nma = 1)
#'model <- esemifar(nonpar_model = nonpar, par_model = par)
#'model
#'

esemifar <- function(nonpar_model, par_model) {
  out <- list(nonpar_model = nonpar_model, par_model = par_model)
  class(out) <- "esemifar"
  attr(out, "function") <- "tsmoothlm"
  attr(out, "type") <- "semifar"
  attr(out, "method") <- "other"
  out
}

generate_future.Farima <- function(object, h, obs) {

  res <- object$residuals
  et <- sample(res - mean(res), size = h, replace = TRUE)

  ar <- farima_to_ar(ar = object$ar, ma = -object$ma, d = object$d, max_i = length(obs) + h - 1)[-1]

  meanObs <- mean(obs)

  out <- c(FARIMAfutureObs(obs = obs, ar_inf = ar, et = et, meanObs = meanObs))

  out

}

point_forecast.Farima <- function(object, h, obs, mean_v) {

  ar <- farima_to_ar(ar = object$ar, ma = -object$ma, d = object$d, max_i = length(obs) + h - 1)[-1]

  out <- c(FARIMAforecastAR(obs = obs, ar_inf = ar, meanObs = mean_v, m = h))

  out

}

adv_boot_creator <- function(ar, ma, d, mean_v, et, obs, h, est_old) {

  p <- length(ar)
  q <- length(ma)
  force(d)
  force(et)
  force(h)
  force(mean_v)
  n <- length(obs)
  force(est_old)
  et_s <- et - mean(et)
  rm(et)

  function() {

    et_new <- sample(et_s, size = 5000 + n + q, replace = TRUE)
    innov_burnin <- utils::head(et_new, 5000)
    innov_n <- et_new[5001:(5000 + n + q)]

    x_new <- suppressWarnings(
      fracdiff::fracdiff.sim(n = n, ar = ar, ma = -ma, d = d,
        innov = innov_n, start.innov = innov_burnin, n.start = 5000)$series +
        mean_v
    )

    mean_v_new <- mean(x_new)

    est_new <- suppressWarnings(fracdiff::fracdiff(x_new, nar = p, nma = q))

    pfc <- point_forecast.Farima(est_new, h = h, obs = obs, mean_v_new)

    obs_new <- generate_future.Farima(est_old, h = h, obs = obs)

    obs_new - pfc

  }
}

#'ESEMIFAR Prediction Method
#'
#'Point and interval forecasts (under the normality assumption or via a
#'bootstrap) for fitted ESEMIFAR models.
#'
#'@param object an object returned by either \code{\link{tsmoothlm}} or
#'\code{\link{esemifar}}.
#'@param n.ahead a single numeric value that represents the forecasting horizon.
#'@param alpha a numeric vector with confidence levels for the forecasting
#'intervals; the default \code{c(0.95, 0.99)} represents 95-percent and
#'99-percent forecasting interval bounds that will be computed.
#'@param method whether to obtain the forecasting intervals under the
#'normality assumption (\code{"norm"}) or via a bootstrap (\code{"boot"}).
#'@param bootMethod only for \code{method = "boot"}: whether to simulate
#'future paths only (\code{"simple"}) or whether to
#'re-estimate the FARIMA model for the re-sampled series and to then obtain
#'simulated predictive roots (\code{"advanced"}).
#'@param npaths only for \code{method = "boot"}: the number of bootstrap
#'iterations.
#'@param quant.type only for \code{method = "boot"}: the quantile type as
#'in the argument \code{type} of the function \code{\link[stats]{quantile}}.
#'@param boot_progress only for \code{method = "boot"}: whether to show a
#'progress bar in the console.
#'@param expo whether to exponentiate all results at the end.
#'@param trend_extrap how to extrapolate the estimated trend into the future:
#'linearly (\code{"linear"}) or constantly (\code{"constant"}).
#'@param future only for \code{method = "boot"}: use parallel programming
#'for the bootstrap via the \code{future} framework?
#'@param num_cores only for \code{method = "boot"} and \code{future = TRUE}:
#'how many cores to use in the parallel programming.
#'@param ... no purpose; for compatibility only.
#'
#'@details
#'Produce point and interval forecasts based on ESEMIFAR models. Throughout,
#'the infinite-order AR-representation of the parametric FARIMA part is considered
#'to produce point forecasts and future paths of the series. The trend is usually
#'extrapolated linearly (or constantly as an alternative).
#'
#'@export
#'
#'@return
#'The function returns a list of class \code{"esemifar"} with elements
#'\code{nonpar_model} and \code{par_model}.
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Scientific Employee) (Department of Economics,
#'Paderborn University), \cr
#'Author
#'}
#'
#'@return
#'A list with various elements is returned.
#'\describe{
#'\item{obs}{the observed series.}
#'\item{mean}{the point forecasts.}
#'\item{lower}{the lower bounds of the forecasting intervals.}
#'\item{upper}{the upper bounds of the forecasting intervals.}
#'\item{model}{the fitted ESEMIFAR model object.}
#'\item{level}{the confidence levels for the forecasting intervals.}
#'}
#'
#'@examples
#'lgdp <- log(esemifar::gdpG7$gdp)
#'est <- tsmoothlm(lgdp, pmax = 1, qmax = 1)
#'# Under normality
#'fc <- predict(est, n.ahead = 10, method = "norm", expo = TRUE)
#'fc$mean
#'fc$lower
#'fc$upper
#'


predict.esemifar <- function(object, n.ahead = 5, alpha = c(0.95, 0.99),
                             method = c("norm", "boot"), bootMethod = c("simple", "advanced"),
                             npaths = 5000, quant.type = 8, boot_progress = TRUE,
                             expo = FALSE,
                             trend_extrap = c("linear", "constant"),
                             future = TRUE,
                             num_cores = future::availableCores() - 1,...) {


  stopifnot("Object cannot be used for prediction." = (inherits(object, "esemifar") & attr(object, "function") == "tsmoothlm"))

  method <- match.arg(method)
  bootMethod <- match.arg(bootMethod)
  trend_extrap <- match.arg(trend_extrap)

  stopifnot(
    "n.ahead must be a positive integer" = (length(n.ahead) == 1 & is.numeric(n.ahead) & n.ahead >= 1),
    "alpha must be a vector with values between 0 and 1" = (length(alpha) > 0 & all(alpha > 0 & alpha < 1)),
    'method must be from either "norm" or "boot"' = (length(method) %in% c(1, 2) & all(method %in% c("norm", "boot"))),
    'bootMethod must be from either "simple" or "advanced"' = (length(bootMethod) %in% c(1, 2) & all(bootMethod %in% c("simple", "advanced"))),
    "npaths must be a single integer value > 0" = (length(npaths) == 1 & npaths > 0),
    "quant.type must be a valid value for the type argument in the quantile function of stats" = (length(quant.type) == 1 & quant.type %in% c(1:9)),
    "expo must be logical" = (length(expo) == 1 & is.logical(expo)),
    'trend_extrap must be either from "linear" or "constant"' = (length(trend_extrap) %in% c(1, 2) & all(trend_extrap %in% c("linear", "constant"))),
    "future must be logical" = (length(future) == 1 & is.logical(future)),
    "num_cores must be an integer > 0" = (length(num_cores) == 1 & num_cores > 0),
    "boot_progress must be logical" = (length(boot_progress) == 1 & is.logical(boot_progress))
  )

  if (attr(object, "type") == "semifar") {
    trend <- object$nonpar_model$ye
    farima <- object$par_model
    resids <- object$nonpar_model$res
    x <- object$nonpar_model$orig
  } else if (attr(object, "type") == "tsmoothlm") {
    trend <- object$ye
    farima <- object$FARMA.BIC
    resids <- object$res
    x <- object$orig
  }
  mean_obs <- mean(resids)

  trend_step <- diff(tail(trend, 2))
  trend.fc <- tail(trend, 1) + (1:n.ahead) * trend_step

  p <- length(farima$ar)
  q <- length(farima$ma)

  if (p > 0) {
    ar <- farima$ar
  } else if (p == 0) {
    ar <- 0
  }

  if (q > 0) {
    ma <- -farima$ma
  } else if (q == 0) {
    ma <- 0
  }
  d <- farima$d
  n <- length(trend)
  m <- n + n.ahead - 1
  h <- n.ahead

  farima.fc <- point_forecast.Farima(object = farima, h = h, resids, mean_obs)
  y.fc <- c(trend.fc) + c(farima.fc)

#-----------------------------------------------------------

  alpha.names <- paste0(100 * alpha, "%")

  if (method == "boot" && bootMethod == "simple") {

    if (future) {

      oldplan <- future::plan()
      future::plan(future::multisession(), workers = num_cores)
      on.exit(future::plan(oldplan), add = TRUE, after = TRUE)

    }

    paths <- matrix(unlist(furrr::future_map(
      1:npaths,
      function(.x, object, h, resids) {
        c(generate_future.Farima(object, h, resids))
      },
      object = farima, h = h, resids = resids,
      .progress = boot_progress,
      .options = furrr::furrr_options(seed = TRUE)
    )), nrow = h, ncol = npaths, byrow = FALSE)

    if (boot_progress) {
      cat("", fill = TRUE)
    }

    fc_low <- t(apply(
      paths,
      MARGIN = 1,
      FUN = stats::quantile,
      probs = 0.5 - (alpha) / 2,
      type = quant.type
    ))

    fc_up <- t(apply(
      paths,
      MARGIN = 1,
      FUN = stats::quantile,
      probs = 0.5 + (alpha) / 2,
      type = quant.type
    ))

  } else if (method == "boot" && bootMethod == "advanced") {

    et <- farima$residuals
    sim_fun <- adv_boot_creator(ar, ma, d, mean_obs, et, resids, h, farima)

    if (future) {

      oldplan <- future::plan()
      future::plan(future::multisession(), workers = num_cores)
      on.exit(future::plan(oldplan), add = TRUE, after = TRUE)

    }

    predRoots <- matrix(unlist(furrr::future_map(
      1:npaths,
      ~ sim_fun(),
      .progress = boot_progress,
      .options = furrr::furrr_options(seed = TRUE)
    )), ncol = h, nrow = npaths, byrow = TRUE)

    if (boot_progress) {
      cat("", fill = TRUE)
    }

    level_low <- 0.5 - alpha / 2
    level_up <- 1 - level_low

    diffs <- t(apply(predRoots, MARGIN = 2, FUN = stats::quantile, type = quant.type, probs = c(level_low, level_up)))

    for (i in 1:h) {
      diffs[i, ] <- diffs[i, ] + farima.fc[[i]]
    }

    fc_low <- diffs[, 1:length(alpha)]
    fc_up <- diffs[, (length(alpha) + 1):(2 * length(alpha))]

  } else if (method == "norm") {       # normal case

    mainf.coef <- farima_to_ma(ar = ar, ma = ma, d = d, max_i = h - 1)

    alpha1 <- 50 - (alpha * 100) / 2
    alpha2 <- 50 + (alpha * 100) / 2
    alpha <- c(alpha1, alpha2)
    alpha <- alpha / 100
    quant.norm <- qnorm(alpha)

    sig2 <- farima$sigma^2

    sd.fcast <- sqrt(sig2 * cumsum(mainf.coef^2))

    l <- length(alpha) / 2

    fc_low <- matrix(unlist(lapply(
      quant.norm[1:l],
      function(.x, sd.fcast, farima.fc) {
        .x * sd.fcast + farima.fc
      }, sd.fcast = sd.fcast, farima.fc = c(farima.fc)
    )), nrow = h, ncol = l)

    fc_up <- matrix(unlist(lapply(
      quant.norm[(l + 1):(2 * l)],
      function(.x, sd.fcast, farima.fc) {
        .x * sd.fcast + farima.fc
      }, sd.fcast = sd.fcast, farima.fc = c(farima.fc)
    )), nrow = h, ncol = l)


  }

  colnames(fc_low) <- alpha.names
  colnames(fc_up) <- alpha.names

  y.fc_low <- c(trend.fc) + fc_low
  y.fc_up <- c(trend.fc) + fc_up
  fitted <- trend + farima$fitted

  if (expo) {
    y.fc <- exp(y.fc)
    y.fc_low <- exp(y.fc_low)
    y.fc_up <- exp(y.fc_up)
    fitted <- exp(fitted)
    x <- exp(x)
  }

  out <- list(
    obs = x,
    mean = y.fc,
    lower = y.fc_low,
    upper = y.fc_up,
    model = object,
    level = alpha
  )

  class(out) <- "esemifar_fc"

  out

}
