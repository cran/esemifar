# esemifar 2.0.1
- a bug in the predict-method was fixed that did not implement the constant
  extrapolation of the trend for `trend_extrap = "constant"`.
- a minor bug in the handling of the residual mean when forecasting was corrected.  
- the README was updated to show examples for the new forecasting approach.
- the main literature in the manual was updated to the by now officially
  published paper.
- a link to the package 'RcppArmadillo' was removed in the manual, as it 
  caused a NOTE on some operating systems.


# esemifar 2.0.0
- predict-method added for output of 'tsmoothlm()' that allows for point
  and interval forecasts (under normality or via bootstrap) based on ESEMIFAR
  models.
- plot functionality updated to allow for plot selection via function argument in
  addition to the interactive console selection.
- functions added to obtain coefficients of different representations of ARMA
  and FARIMA models:
  infinite-order AR representation ('arma_to_ar()' / 'farima_to_ar()'), 
  infinite-order MA representation ('arma_to_ma()' / 'farima_to_ma()'), and
  the infinite coefficient series according to a fractional differencing
  operator ('d_to_coef()').

# esemifar 1.0.1

- implemented a S3 method which extracts fitted values from an 'esemifar' class
  object.
- implemented a S3 method which extracts residuals from an 'esemifar' class
  object.
