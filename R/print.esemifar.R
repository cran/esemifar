#' Print Method for the Package 'esemifar'
#'
#'This function regulates how objects created by the package \code{esemifar} are
#'printed.
#'
#'@param x an input object of class \code{esemifar}.
#'@param ... included for compatibility; additional arguments will however
#'not affect the output.
#'
#'@export
#'
#'@return
#'None
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Scientific employee) (Department of Economics, Paderborn
#'University), \cr
#'}


# Print function for the R package 'esemifar'-------------------------------
print.esemifar <- function(x, ...) {
  if (attr(x, "function") == "tsmoothlm") {
    if (attr(x, "method") == "lpr") {
      cat("-------------------------------------------------", fill = TRUE)
      cat("| Results of the nonparametric trend estimation |", fill = TRUE)
      cat("-------------------------------------------------", fill = TRUE)
      cat("Method: Local Polynomial Regression", fill = TRUE)
      result_vector <- c(as.character(x$n), x$niterations,
                         sprintf("%.4f", x$b0))
      result_dataframe <- data.frame(result_vector)
      rnames_dataframe <- c("Number of observations:",
                            "Iterations until convergence:",
                            "Optimal bandwidth by IPI:")
      colnames(result_dataframe) <- ""
      rownames(result_dataframe) <- rnames_dataframe
      print.data.frame(result_dataframe)
      cat("", fill = TRUE)

      cat("Iterative plug-in algorithm:", fill = TRUE)
      cat("----------------------------")
      ipi_vec <- c(x$bStart, x$p, x$mu, x$InfR,
                   x$bb, x$cb)
      ipi_df <- data.frame(ipi_vec)
      rnames_ipi <- c("Bandwidth starting value:", "Order of polynomial:",
                      "Smoothness parameter:",
                      "Inflation rate:", "Boundary method:",
                      "Boundary cut-off:")
      colnames(ipi_df) <- ""
      rownames(ipi_df) <- rnames_ipi
      print.data.frame(ipi_df)

      cat("", fill = TRUE)
      cat("Components of the object ($):", fill = TRUE)
      cat("-----------------------------")
      abbreviations <- c("ye", "orig", "res", "ws", "b0", "cf0")
      abbr <- data.frame(abbreviations)
      colnames(abbr) <- ""
      rownames(abbr) <- c("Estimates:", "Original series:", "Residuals:",
                          "Weights:", "Optimal bandwidth:",
                          "Estimated variance factor:")
      print.data.frame(abbr, right = FALSE)
      cat(" ", fill = TRUE)
      cat("Iterations:", fill = TRUE)
      cat("-----------", fill = TRUE)
      if (x$niterations < 10) {
        it.names <- paste0("i = ", 1:x$niterations)
      } else {
        it.names <- paste0("i = ", sprintf("%2.f", 1:x$niterations))
      }
      print.data.frame(data.frame(bandwidth = sprintf("%.4f", x$iterations),
                                  row.names = it.names))

    } else if (attr(x, "method") == "kr") {
      cat("-------------------------------------------------", fill = TRUE)
      cat("| Results of the nonparametric trend estimation |", fill = TRUE)
      cat("-------------------------------------------------", fill = TRUE)
      cat("Method: Kernel Regression", fill = TRUE)
      result_vector <- c(as.character(x$n), x$niterations,
                         sprintf("%.4f", x$b0))
      result_dataframe <- data.frame(result_vector)
      rnames_dataframe <- c("Number of observations:",
                            "Iterations until convergence:",
                            "Optimal bandwidth by IPI:")
      colnames(result_dataframe) <- ""
      rownames(result_dataframe) <- rnames_dataframe
      print.data.frame(result_dataframe)
      cat("", fill = TRUE)

      cat("Iterative plug-in algorithm:", fill = TRUE)
      cat("----------------------------")
      ipi_vec <- c(x$bStart, x$p, x$mu, x$InfR,
                   x$bb, x$cb)
      ipi_df <- data.frame(ipi_vec)
      rnames_ipi <- c("Bandwidth starting value:", "Order of polynomial:",
                      "Smoothness parameter:",
                      "Inflation rate:", "Boundary method:",
                      "Boundary cut-off:")
      colnames(ipi_df) <- ""
      rownames(ipi_df) <- rnames_ipi
      print.data.frame(ipi_df, right = FALSE)

      cat("", fill = TRUE)
      cat("Components of the object ($):", fill = TRUE)
      cat("-----------------------------")
      abbreviations <- c("ye", "orig", "res", "b0", "cf0")
      abbr <- data.frame(abbreviations)
      colnames(abbr) <- ""
      rownames(abbr) <- c("Estimates:", "Original series:", "Residuals:",
                          "Optimal bandwidth:", "Estimated variance factor:")
      print.data.frame(abbr, right = FALSE)
      cat(" ", fill = TRUE)
      cat("Iterations:", fill = TRUE)
      cat("-----------", fill = TRUE)
      if (x$niterations < 10) {
        it.names <- paste0("i = ", 1:x$niterations)
      } else {
        it.names <- paste0("i = ", sprintf("%2.f", 1:x$niterations))
      }
      print.data.frame(data.frame(bandwidth = sprintf("%.4f", x$iterations),
                                  row.names = it.names))

    } else {
      print.default(x)
    }
  } else if(attr(x, "function") == "dsmoothlm") {
    cat("------------------------------------------------------",
        fill = TRUE)
    cat("| Results of the nonparametric derivative estimation |",
        fill = TRUE)
    cat("------------------------------------------------------",
        fill = TRUE)
    cat("Method: Local Polynomial Regression", fill = TRUE)
    result_vector <- c(as.character(x[["v"]]), x$n, x$niterations,
                       sprintf("%.4f", x$b0))
    result_dataframe <- data.frame(result_vector)
    rnames_dataframe <- c("Order of derivative:", "Number of observations:",
                          "Iterations until convergence:",
                          "Optimal bandwidth by IPI:")
    colnames(result_dataframe) <- ""
    rownames(result_dataframe) <- rnames_dataframe
    print.data.frame(result_dataframe)
    cat("", fill = TRUE)
    cat("Iterative plug-in algorithm:", fill = TRUE)
    cat("----------------------------")

    ipi_vec <- c(round(x$bStart, 4), x$pp, x$p, x$mu.p, x$InfR.p)
    ipi_df <- data.frame(ipi_vec)
    rnames_ipi <- c("Bandwidth starting value (pilot-IPI):",
                    "Order of polynomial (pilot-IPI):",
                    "Order of polynomial (IPI):",
                    "Smoothness parameter (pilot-IPI):",
                    "Inflation rate (pilot-IPI):")
    colnames(ipi_df) <- ""
    rownames(ipi_df) <- rnames_ipi
    print.data.frame(ipi_df)

    cat("", fill = TRUE)
    cat("Components of the object ($):", fill = TRUE)
    cat("-----------------------------")
    abbreviations <- c("ye", "orig", "ws", "b0", "cf0")
    abbr <- data.frame(abbreviations)
    colnames(abbr) <- ""
    rownames(abbr) <- c("Estimates:", "Original series:", "Weights:",
                        "Optimal bandwidth:", "Estimated variance factor:")
    print.data.frame(abbr, right = FALSE)
    cat(" ", fill = TRUE)
    cat("Iterations:", fill = TRUE)
    cat("-----------", fill = TRUE)
    if (x$niterations < 10) {
      it.names <- paste0("i = ", 1:x$niterations)
    } else {
      it.names <- paste0("i = ", sprintf("%2.f", 1:x$niterations))
    }
    print.data.frame(data.frame(bandwidth = sprintf("%.4f", x$iterations),
                                row.names = it.names))
  }
}

