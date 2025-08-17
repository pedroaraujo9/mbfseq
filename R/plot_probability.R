 #' Plot Posterior Stage Probabilities
 #'
 #' Visualize posterior stage probabilities over time for each group, with credible intervals as ribbons.
 #'
 #' @param posterior_summary Output from posterior_summary containing summarized probabilities and intervals.
 #'
 #' @return A ggplot2 object showing posterior probabilities and credible intervals for each group over time.
 #'
 #' @details
 #' Plots the mean and credible interval (ribbon) for each category and group, using ggplot2 and dplyr.
 #'
 #' @examples
 #' \dontrun{
 #' # Plot posterior probabilities from a summary
 #' plot_probability(posterior_summary)
 #' }
 #'
 #' @importFrom dplyr mutate
 #' @importFrom ggplot2 ggplot aes geom_line geom_ribbon labs facet_grid
 #' @export
plot_probability = function(posterior_summary) {
  posterior_summary$prob_stage_summ |>
    dplyr::mutate(w = paste0("Group ", w)) |>
    ggplot2::ggplot(ggplot2::aes(x = time, y = mean, color = cat)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = li, ymax = ui, fill = cat, y = mean, x = time),
                         alpha = 0.3, inherit.aes = FALSE) +
    ggplot2::labs(x = "Time", y = "Probability", fill = "Z", color = "Z") +
    ggplot2::facet_grid(. ~ w)
}
