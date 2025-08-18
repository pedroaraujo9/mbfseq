#' Plot Clustered Sequences
#'
#' Visualize the sequence clustering results as a heatmap, with cluster boundaries highlighted. Accepts either a fitted model object or raw clustering assignments.
#'
#' @param fit Optional. Output from fit_mbfseq. If provided, extracts clustering assignments and metadata for plotting.
#' @param G Optional. Number of clusters (G). Used to select the model from fit if fit is provided.
#' @param M Optional. Number of mixture components (M). Used to select the model from fit if fit is provided.
#' @param z Optional. Sequence cluster assignments. Used if fit is not provided.
#' @param w Optional. Cluster labels for each sequence. Used if fit is not provided.
#' @param id Vector of subject or sequence identifiers. Used if fit is not provided.
#' @param time Vector of time points for each observation. Used if fit is not provided.
#'
#' @return A ggplot2 object showing a heatmap of sequence clusters over time, with cluster boundaries marked in red.
#'
#' @details
#' If a fitted model object is provided, the function extracts cluster assignments and metadata automatically. Otherwise, provide z, w, id, and time directly.
#'
#' @examples
#' \dontrun{
#' # Using a fitted model:
#' plot_seq_cluster(fit = result, G = 2, M = 4)
#'
#' # Using raw assignments:
#' plot_seq_cluster(z = z, w = w, id = id, time = time)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile geom_hline labs theme_minimal scale_y_discrete scale_x_continuous geom_text coord_cartesian theme margin geom_segment expansion
#' @importFrom dplyr arrange group_by slice_tail mutate filter summarise left_join
#' @importFrom viridis scale_fill_viridis
#' @importFrom stats median
#' @export
plot_seq_cluster = function(fit = NULL,
                            G = NULL,
                            M = NULL,
                            z = NULL,
                            w = NULL,
                            id = NULL,
                            time = NULL) {

  if(!is.null(fit)) {
    model_name = paste0("G=", G, ", M=", M)
    model = fit$models[[model_name]]
    post_summ = fit %>% posterior_summary(M = M, G = G)
    z_levels = fit$model_data$z_levels
    w = model$w_class
    z = model$z_class %>% factor(labels = z_levels)
    id = fit$args$id
    time = fit$args$time

    z_mode = post_summ[[model_name]]$prob_stage_summ %>%
      dplyr::group_by(time, w) %>%
      dplyr::summarise(z_mode = cat[which.max(mean)], .groups = "drop") %>%
      dplyr::arrange(w, time)

    cluster_init_end = data.frame(w = w, id = names(w)) %>%
      dplyr::arrange(w) %>%
      dplyr::mutate(id = factor(id, levels = id)) %>%
      dplyr::arrange(w, id) %>%
      dplyr::group_by(w) %>%
      dplyr::summarise(min_id = min(as.integer(id)), max_id = max(as.integer(id)))

    # Create a variable for mode difference
    mode_change = z_mode %>%
      dplyr::group_by(w) %>%
      dplyr::mutate(mode_diff = c(0, z_mode %>% as.integer() %>% diff())) %>%
      as.data.frame() %>%
      tibble::as_tibble() %>%
      dplyr::filter(mode_diff != 0) %>%
      dplyr::select(time, w) %>%
      dplyr::left_join(cluster_init_end, by = c("w"))

    model_change_gg = ggplot2::geom_segment(
      data = mode_change,
      ggplot2::aes(x = time - 0.5, xend = time - 0.5, y = min_id - 0.5, yend = max_id + 0.5),
      inherit.aes = F,
      color = "red", linetype = 1
    )

  }else{
    model_change_gg = ggplot2::theme_minimal()
  }

  id_label = id %>% unique()

  cut_point = data.frame(w = w, id = names(w)) %>%
    dplyr::mutate(id = id %>% factor(levels = id_label)) %>%
    dplyr::arrange(w, id) %>%
    dplyr::group_by(w) %>%
    dplyr::slice_tail(n = 1) %>%
    .$id

  mid_points = data.frame(w = w, id = names(w)) %>%
    dplyr::arrange(w) %>%
    dplyr::mutate(id = factor(id, levels = id)) %>%
    dplyr::group_by(w) %>%
    dplyr::summarise(mid_point = round(stats::median(as.integer(id)), 0)) %>%
    dplyr::mutate(w = paste0("Cluster ", w))

  data.frame(id = id, time = time, z = z) %>%
    dplyr::mutate(id = factor(id, levels = names(sort(w)))) %>%
    ggplot2::ggplot(ggplot2::aes(x=time, y=id, fill=factor(z))) +
    ggplot2::geom_tile(color="gray20") +
    viridis::scale_fill_viridis(discrete = T) +
    ggplot2::labs(x="Time", y="Id", fill=expression(Z[it])) +
    ggplot2::geom_hline(yintercept = which(unique(names(sort(w))) %in% cut_point) + 0.5,
                        color = "red", linewidth = 1) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(0)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0)) +
    model_change_gg +
    ggplot2::theme_minimal() +
    ggplot2::geom_text(
      data = mid_points,
      ggplot2::aes(x = max(time), y = mid_point, label = w),
      inherit.aes = FALSE,
      hjust = -0.3,
      angle = 0
    ) +
    ggplot2::coord_cartesian(clip = 'off') + # allow text outside plot
    ggplot2::theme(
      legend.position = "top",
      plot.margin = ggplot2::margin(5, 60, 5, 5) # extra space on right
    )

}
