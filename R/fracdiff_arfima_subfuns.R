fracdiff_preset <- function(x, nar, nma) {
  suppressWarnings(
    fracdiff::fracdiff(x, nar = nar, nma = nma,
      drange = c(0, 0.5))
  )
}

# arfima_preset <- function(x, nar, nma) {
#   suppressWarnings(
#     arfima::arfima(z = x, order = c(nar, 0, nma),
#                    quiet = TRUE)
#   )
# }

ar_select <- function(fit, method_error, p.BIC) {

  switch(
    method_error,
    "fracdiff" = fit$ar,
    "arfima" = utils::head(summary(fit)$coef[[1]], p.BIC)
  )

}

ma_select <- function(fit, method_error, p.BIC, q.BIC) {

  switch(
    method_error,
    "fracdiff" = fit$ma,
    "arfima" = summary(fit)$coef[[1]][(p.BIC + 1):(p.BIC + q.BIC)]
  )

}

d_select <- function(fit, method_error, p.BIC, q.BIC) {

  switch(
    method_error,
    "fracdiff" = fit$d,
    "arfima" = summary(fit)$coef[[1]][[p.BIC + q.BIC + 1]]
  )

}

sigma2_select <- function(fit, method_error) {

  switch(
    method_error,
    "fracdiff" = fit$sigma^2,
    "arfima" = summary(fit)$sigma2[[1]]
  )

}
