cf0.FARMA.est <- function(Xt, pmin, pmax, qmin, qmax) {
  p.farma <- rep(pmin:pmax, times = qmax - qmin + 1)
  q.farma <- rep(qmin:qmax, each = pmax - pmin + 1)

  # error_fun <- switch(
  #   method_error,
  #   "fracdiff" = fracdiff_preset,
  #   "arfima" = arfima_preset
  # )

  bic <- mapply(function(x, p0, q0) {
    farma.nobs <- length(x)
    # farma.est <- error_fun(x, nar = p0, nma = q0)
    farma.est <- fracdiff_preset(x, nar = p0, nma = q0)
    stats::AIC(farma.est, k = log(farma.nobs))
  }, p0 = p.farma, q0 = q.farma,
  MoreArgs = list(x = Xt))
  opt <- which(bic == min(bic))
  p.BIC <- p.farma[opt]
  q.BIC <- q.farma[opt]
  # FARMA.BIC <- error_fun(Xt, nar = p.BIC, nma = q.BIC)
  FARMA.BIC <- fracdiff_preset(Xt, nar = p.BIC, nma = q.BIC)

#-------------------------------

  if (p.BIC != 0) {
    # ar <- ar_select(FARMA.BIC, method_error, p.BIC)
    ar <- ar_select(FARMA.BIC, "fracdiff", p.BIC)
    sc.AR <- sum(ar)
  } else {
    ar <- 0
    sc.AR <- 0
  }
  if (q.BIC != 0) {
    # ma <- ma_select(FARMA.BIC, method_error, p.BIC, q.BIC)
    ma <- ma_select(FARMA.BIC, "fracdiff", p.BIC, q.BIC)
    sc.MA <- sum(ma)
  } else {
    ma <- 0
    sc.MA <- 0
  }
  # d.BIC <- d_select(FARMA.BIC, method_error, p.BIC, q.BIC)
  d.BIC <- d_select(FARMA.BIC, "fracdiff", p.BIC, q.BIC)

  if((p.BIC == 0) & (q.BIC == 0)) {
    # Better like this using the designated function for computing the
    # fractionally differenced series
    r <- fracdiff::diffseries(Xt, d.BIC)
    n <- length(r)
    sigma2 <- sum(r**2) / n
  } else {
    # sigma2 <- sigma2_select(FARMA.BIC, method_error)
    sigma2 <- sigma2_select(FARMA.BIC, "fracdiff")
  }

  cf0.FARMA <- ((1 - sc.MA) / (1 - sc.AR)) ^ 2 * sigma2 / (2 * pi)
  results <- list(cf0.FARMA = cf0.FARMA, FARMA.BIC = FARMA.BIC,
                  p.BIC = p.BIC, q.BIC = q.BIC, d.BIC = d.BIC)
  class(results) <- "esemifar"
  attr(results, "function") <- "cf0.FARMA.est"
  results
}
