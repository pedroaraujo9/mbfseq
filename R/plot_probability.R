#' Plot Posterior Probability Trajectories
#'
#' Visualize the posterior mean and credible intervals of group probabilities over time from a posterior summary object.
#'
#' @param posterior_summary A list as returned by \code{posterior_summary}, containing \code{prob_summ_df}.
#'
#' @return A \code{ggplot} object displaying the posterior probability trajectories and credible intervals.
#'
#' @details
#' This function creates a faceted line plot of posterior mean probabilities over time for each group, with shaded ribbons representing credible intervals. Each facet corresponds to a group, and colors indicate categories or latent classes.
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs facet_grid
#' @examples
#' \dontrun{
#' # Assume 'fit' is a fitted model object and posterior_summary(fit) has been run
#' summ <- posterior_summary(fit)
#' plot_probability(summ)
#' }
#' @export
plot_probability = function(posterior_summary) {
  posterior_summary$prob_summ_df |>
    dplyr::mutate(w = paste0("Group ", w)) |>
    ggplot2::ggplot(ggplot2::aes(x = time, y = mean, color = cat)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = li, ymax = ui, fill = cat, y = mean, x = time),
                         alpha = 0.3, inherit.aes = FALSE) +
    ggplot2::labs(x = "Time", y = "Probability", fill = "Z", color = "Z") +
    ggplot2::facet_grid(. ~ w)
}
