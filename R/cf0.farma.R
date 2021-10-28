cf0.FARMA.est <- function(Xt, pmin, pmax, qmin, qmax) {
  p.farma <- rep(pmin:pmax, times = qmax - qmin + 1)
  q.farma <- rep(qmin:qmax, each = pmax - pmin + 1)
  bic <- mapply(function(x, p0, q0) {
    farma.nobs <- length(x)
    farma.est <- suppressWarnings(fracdiff::fracdiff(x, nar = p0, nma = q0,
                                                     drange = c(0, 0.5)))
    -2 * farma.est$log.likelihood + log(farma.nobs) * (p0 + q0)
  }, p0 = p.farma, q0 = q.farma,
  MoreArgs = list(x = Xt))
  opt <- which(bic == min(bic))
  p.BIC <- p.farma[opt]
  q.BIC <- q.farma[opt]
  FARMA.BIC <- suppressWarnings(fracdiff::fracdiff(Xt, nar = p.BIC, nma = q.BIC,
                                                  drange = c(0, 0.5)))
  if (p.BIC != 0) {
    ar <- FARMA.BIC$ar
    sc.AR <- sum(ar)
  } else {
    ar <- 0
    sc.AR <- 0
  }
  if (q.BIC != 0) {
    ma <- FARMA.BIC$ma
    sc.MA <- sum(ma)
  } else {
    ma <- 0
    sc.MA <- 0
  }
  d.BIC <- FARMA.BIC$d
  # r <- Xt
  # n <- length(r)
  # b <- dfrac(n, d.BIC)$b
  # r <- c(r * 0, r)
  # r <- stats::filter(r, b, sides = 1)[(n + 1):(2 * n)]
  r <- fracdiff::fracdiff(Xt, drange = c(d.BIC, d.BIC))$residuals # alternativ

  if((p.BIC == 0) & (q.BIC == 0))
  {
    n <- length(r)
    sigma2 <- sum(r**2) / n
  }
  else
  {
    coefest <- c(ar, -ma)[which(c(ar, -ma) != 0)]
    #sigma2 <- stats::arima(r, order = c(p.BIC, 0, q.BIC))$sigma2
    sigma2 <- stats::arima(r, order = c(p.BIC, 0, q.BIC), fixed = c(coefest),
                           include.mean = FALSE)$sigma2
  }

  cf0.FARMA <- ((1 - sc.MA) / (1 - sc.AR)) ^ 2 * sigma2 / (2 * pi)
  results <- list(cf0.FARMA = cf0.FARMA, FARMA.BIC = FARMA.BIC,
                  p.BIC = p.BIC, q.BIC = q.BIC, d.BIC = d.BIC)
  class(results) <- "esemifar"
  attr(results, "function") <- "cf0.FARMA.est"
  results
}
