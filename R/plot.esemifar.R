#' Plot Method for the Package 'esemifar'
#'
#'This function regulates how objects created by the package \code{esemifar} are
#'plotted.
#'
#'@param x an input object of class \code{esemifar}.
#'@param t an optional vector with time points that will be considered for
#'the x-axis within the plot; is set to NULL by default and uses a vector
#'\code{1:length(x$ye)} for time points.
#'@param ... additional arguments of the standard plot method.
#'
#'@export
#'
#'@importFrom graphics lines
#'
#'@return
#'None
#'
#'@author
#'\itemize{
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'}
#'

plot.esemifar <- function(x, t = NULL, ...) {

  dots <- list(...)
  dots[["x"]] <- t
  dots[["y"]] <- x
  if (attr(dots$y, "function") == "tsmoothlm") {
    cat("Plot choices for esemifar object:", fill = TRUE)
    choices <- c(1, 2, 3, 4)
    choice_names <- c("Original series:", "Trend series:", "Residual series:",
                      "Original with trend series:")
    choices_df <- data.frame(choices)
    colnames(choices_df) <- ""
    rownames(choices_df) <- choice_names
    print.data.frame(choices_df)
    plot_choice <- readline(prompt="Please enter the corresponding number: ")
    plot_choice <- as.numeric(plot_choice)

    if (plot_choice == 1) {
      if (is.null(dots[["main"]])) {
        dots[["main"]] <- "Original series"
      }
      if (is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- "Time"
      }
      if (is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- "y"
      }
      if (is.null(dots[["type"]])) {
        dots[["type"]] <- "l"
      }
      dots[["y"]] <- dots[["y"]]$orig
      n <- length(dots[["y"]])
      if (is.null(dots[["x"]])) {
        dots[["x"]] <- 1:n
      }

      if (is.null(dots[["xlim"]])) {
        dots[["xlim"]] <- c(dots[["x"]][[1]], dots[["x"]][[n]])
      }
      if (is.null(dots[["ylim"]])) {
        x.low <- sum(dots[["x"]] < dots[["xlim"]][[1]]) + 1
        x.up <- sum(dots[["x"]] <= dots[["xlim"]][[2]])
        dots[["ylim"]] <- c(min(dots[["y"]][x.low:x.up]),
                            max(dots[["y"]][x.low:x.up]))
      }

      do.call(what = graphics::plot, args = dots)

    } else if (plot_choice == 2) {
      if (is.null(dots[["main"]])) {
        dots[["main"]] <- "Estimated trend series"
      }
      if (is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- "Time"
      }
      if (is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- "ye"
      }
      if (is.null(dots[["type"]])) {
        dots[["type"]] <- "l"
      }
      dots[["y"]] <- dots[["y"]]$ye
      n <- length(dots[["y"]])
      if (is.null(dots[["x"]])) {
        dots[["x"]] <- 1:n
      }

      if (is.null(dots[["xlim"]])) {
        dots[["xlim"]] <- c(dots[["x"]][[1]], dots[["x"]][[n]])
      }
      if (is.null(dots[["ylim"]])) {
        x.low <- sum(dots[["x"]] < dots[["xlim"]][[1]]) + 1
        x.up <- sum(dots[["x"]] <= dots[["xlim"]][[2]])
        dots[["ylim"]] <- c(min(dots[["y"]][x.low:x.up]),
                            max(dots[["y"]][x.low:x.up]))
      }

      do.call(what = graphics::plot, args = dots)

    } else if (plot_choice == 3) {
      if (is.null(dots[["main"]])) {
        dots[["main"]] <- "Residual series"
      }
      if (is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- "Time"
      }
      if (is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- "resid"
      }
      if (is.null(dots[["type"]])) {
        dots[["type"]] <- "l"
      }
      dots[["y"]] <- dots[["y"]]$res
      n <- length(dots[["y"]])
      if (is.null(dots[["x"]])) {
        dots[["x"]] <- 1:n
      }

      if (is.null(dots[["xlim"]])) {
        dots[["xlim"]] <- c(dots[["x"]][[1]], dots[["x"]][[n]])
      }
      if (is.null(dots[["ylim"]])) {
        x.low <- sum(dots[["x"]] < dots[["xlim"]][[1]]) + 1
        x.up <- sum(dots[["x"]] <= dots[["xlim"]][[2]])
        dots[["ylim"]] <- c(min(dots[["y"]][x.low:x.up]),
                            max(dots[["y"]][x.low:x.up]))
      }

      do.call(what = graphics::plot, args = dots)

    } else if (plot_choice == 4) {
      if (is.null(dots[["col"]])) dots[["col"]] <- "black"
      color.n <- color.name(dots[["col"]])
      if (is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- "Time"
      }
      if (is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- "y & ye"
      }
      if (is.null(dots[["type"]])) {
        dots[["type"]] <- "l"
      }
      if (is.null(dots[["main"]])) {
        dots[["main"]] <- paste0("Original (", color.n,
                                 ") with trend (red) series")
      }
      ye <- dots[["y"]]$ye
      dots[["y"]] <- dots[["y"]]$orig
      n <- length(dots[["y"]])
      if (is.null(dots[["x"]])) {
        dots[["x"]] <- 1:n
      }
      if (is.null(dots[["xlim"]])) {
        dots[["xlim"]] <- c(dots[["x"]][[1]], dots[["x"]][[n]])
      }
      if (is.null(dots[["ylim"]])) {
        x.low <- sum(dots[["x"]] < dots[["xlim"]][[1]]) + 1
        x.up <- sum(dots[["x"]] <= dots[["xlim"]][[2]])
        dots[["ylim"]] <- c(min(dots[["y"]][x.low:x.up], ye[x.low:x.up]),
                            max(dots[["y"]][x.low:x.up], ye[x.low:x.up]))
      }

      do.call(what = graphics::plot, args = dots)
      lines(dots[["x"]], ye, col = "red")

    }
  } else if (attr(dots[["y"]], "function") == "dsmoothlm") {

    cat("Plot choices for smoots object:", fill = TRUE)
    choices <- c(1, 2)
    choice_names <- c("Original series:", "Derivative series:")
    choices_df <- data.frame(choices)
    colnames(choices_df) <- ""
    rownames(choices_df) <- choice_names
    print.data.frame(choices_df)
    plot_choice <- readline(prompt="Please enter the corresponding number: ")
    plot_choice <- as.numeric(plot_choice)

    if (plot_choice == 1) {
      if (is.null(dots[["main"]])) {
        dots[["main"]] <- "Original series"
      }
      if (is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- "Time"
      }
      if (is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- "y"
      }
      if (is.null(dots[["type"]])) {
        dots[["type"]] <- "l"
      }
      dots[["y"]] <- dots[["y"]]$orig
      n <- length(dots[["y"]])
      if (is.null(dots[["x"]])) {
        dots[["x"]] <- 1:n
      }

      if (is.null(dots[["xlim"]])) {
        dots[["xlim"]] <- c(dots[["x"]][[1]], dots[["x"]][[n]])
      }
      if (is.null(dots[["ylim"]])) {
        x.low <- sum(dots[["x"]] < dots[["xlim"]][[1]]) + 1
        x.up <- sum(dots[["x"]] <= dots[["xlim"]][[2]])
        dots[["ylim"]] <- c(min(dots[["y"]][x.low:x.up]),
                            max(dots[["y"]][x.low:x.up]))
      }

      do.call(what = graphics::plot, args = dots)

    } else if (plot_choice == 2) {
      if (is.null(dots[["main"]])) {
        if (attr(dots[["y"]], "function") == "dsmoothlm") {
          dots[["main"]] <- paste0("Estimated derivative series, d = ",
                                   dots[["y"]]$v)
        } else if (attr(dots[["y"]], "function") == "gsmooth") {
          dots[["main"]] <- paste0("Estimated derivative series, d = ",
                                   dots[["y"]]$v)
        }
      }
      if (is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- "Time"
      }
      if (is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- "ye"
      }
      if (is.null(dots[["type"]])) {
        dots[["type"]] <- "l"
      }
      model <- dots[["y"]]
      dots[["y"]] <- dots[["y"]]$ye
      n <- length(dots[["y"]])
      if (is.null(dots[["x"]])) {
        dots[["x"]] <- 1:n
      }

      if (is.null(dots[["xlim"]])) {
        dots[["xlim"]] <- c(dots[["x"]][[1]], dots[["x"]][[n]])
      }
      if (is.null(dots[["ylim"]])) {
        x.low <- sum(dots[["x"]] < dots[["xlim"]][[1]]) + 1
        x.up <- sum(dots[["x"]] <= dots[["xlim"]][[2]])
        dots[["ylim"]] <- c(min(dots[["y"]][x.low:x.up]),
                            max(dots[["y"]][x.low:x.up]))
      }

      do.call(what = graphics::plot, args = dots)

    }
  }
}
