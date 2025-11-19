#' Get Predictions from Stan Model Fit
#'
#' Extracts predicted values from a Stan model fit.
#'
#' @param stan_fit A \code{stanfit} object returned by \code{get_stanfit}.
#' @return A summary of the predictions and bar charts of each prediction.
#'
#' @examples
#' \donttest{
#'   # Generate toy data
#'   y_data <- data_simu(n = 60, alpha = 0.5, rho = 0.3, theta = c(5),
#'                       mod_type = "zi", distri = "poi")
#'
#'   # Fit a small Stan model (may take > 5s on first compile)
#'   stan_fit <- get_stanfit(mod_type = "zi", distri = "poi", y = y_data)
#'
#'   # Get predicted values from the Stan model fit
#'   get_pred(stan_fit = stan_fit)
#' }
#'
#' @import rstan
#' @importFrom ggplot2 ggplot aes geom_bar theme_minimal labs theme element_text
#' @importFrom gridExtra grid.arrange
#' @importFrom knitr kable
#' @importFrom stats median IQR
#' @importFrom utils globalVariables
#' @export



get_pred <- function(stan_fit) {

  if (!inherits(stan_fit, "stanfit")) {
    stop("The parameter 'stan_fit' must be a valid 'stanfit' object
         returned by Stan.")
  }

  y_pred <- data.frame(rstan::extract(stan_fit, pars = "y_pred"))

  y_pred_summary <- data.frame(
    Mode = apply(y_pred, 2, get_mode),
    Median = apply(y_pred, 2, median),
    IQR = apply(y_pred, 2, IQR),
    Min = apply(y_pred, 2, min),
    Max = apply(y_pred, 2, max)
  )

  num_cols <- ncol(y_pred)
  plot_list <- list()

  for (i in 1:num_cols) {
    column_data <- as.data.frame(table(y_pred[, i]))
    colnames(column_data) <- c("Value", "Frequency")

    p <- ggplot2::ggplot(column_data, aes(x = as.factor(Value),
                                          y = Frequency)) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(
        title = paste("Bar Chart for Prediction", i),
        x = "Predicted Value",
        y = "Frequency"
      )
    plot_list[[i]] <- p
  }

  return(list(summary = y_pred_summary, plots = plot_list))
}
