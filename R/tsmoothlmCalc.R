
tsmoothlmCalc <- function(y,
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

  # Input parameters
  n <- length(y)
  k <- p + 1
  pd <- p + 2
  runc <- 1
  n1 <- trunc(n * cb)

  steps <- rep(NA, 40)
  bd_func <- lookup$InfR_lookup[as.character(p), InfR][[1]]


  # Some constants
  m <- 1000000  # for the numerical integral
  u <- (-m:m) / (m + 0.5)

  if (p == 1) {
    wkp <- (1 - u^2)^(mu)
    lfun_intern <- lookup$p1p3_lookup[mu + 1, 1][[1]]
    bopt_calc <- function(c1, c2, c3, d.BIC, n) {
      bopt <- (c1 * c2 * c3)^(1 / (5 - 2 * d.BIC)) * n^((2 * d.BIC - 1)
                                                        / (5 - 2 * d.BIC))
      if (bopt < n^(-5 / 7)) {bopt <- n^(-5 / 7)}
      bopt
    }
  } else if (p == 3) {
    wkp <- lookup$p3_lookup[mu + 1][[1]](u)
    lfun_intern <- lookup$p1p3_lookup[mu + 1, 2][[1]]
    bopt_calc <- function(c1, c2, c3, d.BIC, n) {
      bopt <- (c1 * c2 * c3)^(1 / (9 - 2 * d.BIC)) * n^((2 * d.BIC - 1)
                                                        / (9 - 2*d.BIC))
      if (bopt < n^(-9 / 11)) {bopt <- n^(-9 / 11)}
      bopt
    }
  }
  Rp_calc <- function(Vc, d.BIC) {
    Vc * lfun_intern(d.BIC)
  }
  mukp <- sum((u^k) * wkp) / m ### beta
  c1 <- factorial(k)^2 / (2 * k)  # a constant in the bandwidth

  # The main iteration------------------------------------------------------

  noi <- 40  # maximal number of iterations

  for (i in 1:noi) {
    if (runc == 1) {
      if (i > 1) {bold1 <- bold}
      if (i == 1) {bold <- bStart} else {bold <- bopt}

      ye <- smoots::gsmooth(y, 0, p, mu, bold, bb)$ye

      # The optimal bandwidth-----------------------------------------------
      # Estimating the variance factor
      yd <- y - ye

      cf0.FARMA <- cf0.FARMA.est(yd, pmin, pmax, qmin, qmax)
      cf0 <- cf0.FARMA$cf0.FARMA
      p.BIC <- cf0.FARMA$p.BIC
      q.BIC <- cf0.FARMA$q.BIC
      d.BIC <- cf0.FARMA$d.BIC

      #
      # Look up the EIM inflation rate in the internal list 'lookup'

      bd <- bd_func(bold, d.BIC)

      if (bd >= 0.49) {bd <- 0.49}

      yed <- smoots::gsmooth(y, k, pd, mu, bd, bb)$ye

      I2 <- sum(yed[max(1, n1):((1 - cb) * n)]^2) / (n - 2 * n1)
      Vc <- 2 * gamma(1 - 2 * d.BIC) * sin(pi * d.BIC) # Const. in the variance
      Rp <- Rp_calc(Vc = Vc, d.BIC = d.BIC)

      # Three constants in the bandwidth
      c2 <- (1 - 2 * cb) * (1 - 2 * d.BIC) * Rp / (mukp)^2
      c3 <- cf0 / I2

      bopt <- bopt_calc(c1 = c1, c2 = c2, c3 = c3, d.BIC = d.BIC, n = n)
      if (bopt > 0.49) {bopt <- 0.49}
      steps[i] <- bopt
      if (i > 2 && abs(bold - bopt) / bopt < 1 / n) {runc = 0}
      if (i > 3 && abs(bold1 - bopt) / bopt < 1 / n) {
        bopt <- (bold + bopt) / 2
        runc <- 0
      }
    }
  }

  # Smooth with the selected bandwidth--------------------------------------
  if (p == 1 && bopt < n^(-5 / 7)) {bopt <- n^(-5 / 7)}
  if (p == 3 && bopt < n^(-9 / 11)) {bopt <- n^(-9 / 11)}
  if (bopt > 0.49) {bopt <- 0.49}

  if (method == "lpr") {

    est.opt <- smoots::gsmooth(y, 0, p, mu, bopt, bb)
    ye <- est.opt$ye
    ws <- est.opt$ws
    res <- y - ye

    cf0.FARMA <- cf0.FARMA.est(res, pmin, pmax, qmin, qmax)
    cf0 <- cf0.FARMA$cf0.FARMA
    p.BIC <- cf0.FARMA$p.BIC
    q.BIC <- cf0.FARMA$q.BIC
    d.BIC <- cf0.FARMA$d.BIC
    FARMA.BIC <- cf0.FARMA$FARMA.BIC

    results <- list(ye = ye, b0 = bopt, cf0 = cf0, I2 = I2,
                    cf0 = cf0, p.BIC = p.BIC, q.BIC = q.BIC, d.BIC = d.BIC,
                    FARMA.BIC = FARMA.BIC,
                    n = n, niterations = length(steps[!is.na(steps)]),
                    res = res, orig = y, iterations = steps[!is.na(steps)],
                    p = p, mu = mu, InfR = InfR, bStart = bStart,
                    bb = bb, cb = cb, ws = ws)
    attr(results, "method") <- "lpr"
  } else if (method == "kr") {
    est.opt <- smoots::knsmooth(y, mu, bopt, bb)
    res <- est.opt$res
    cf0.FARMA <- cf0.FARMA.est(res, pmin, pmax, qmin, qmax)
    cf0 <- cf0.FARMA$cf0.FARMA
    p.BIC <- cf0.FARMA$p.BIC
    q.BIC <- cf0.FARMA$q.BIC
    d.BIC <- cf0.FARMA$d.BIC
    FARMA.BIC <- cf0.FARMA$FARMA.BIC
    results <- list(ye = est.opt$ye, b0 = bopt, cf0 = cf0, I2 = I2,
                    cf0 = cf0, p.BIC = p.BIC, q.BIC = q.BIC, d.BIC = d.BIC,
                    FARMA.BIC = FARMA.BIC,
                    n = n, niterations = length(steps[!is.na(steps)]),
                    res = res, orig = y, iterations = steps[!is.na(steps)],
                    p = p, mu = mu, InfR = InfR, bStart = bStart,
                    bb = bb, cb = cb)
    attr(results, "method") <- "kr"
  }
  results
}
# End of the function
