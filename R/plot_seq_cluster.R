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
 #' @importFrom ggplot2 ggplot aes geom_tile geom_hline labs theme_minimal scale_y_discrete scale_x_continuous
 #' @importFrom dplyr arrange group_by slice_tail mutate filter
 #' @importFrom viridis scale_fill_viridis
 #' @export
plot_seq_cluster = function(fit = NULL,
                            G = NULL,
                            M = NULL,
                            z = NULL,
                            w = NULL,
                            id = NULL,
                            time = NULL) {

  if(!is.null(fit)) {
    model = fit$models[[paste0("G=", G, ", M=", M)]]
    w = model$w_class
    z = model$z_class
    id = fit$args$id
    time = fit$args$time
  }

  id_label = id %>% unique()

  cut_point = data.frame(w = w, id = names(w)) %>%
    mutate(id = id %>% factor(id, levels = id_label)) %>%
    dplyr::arrange(w, id) %>%
    dplyr::group_by(w) %>%
    dplyr::slice_tail(n = 1) %>%
    .$id

  data.frame(id = id, time = time, z = z) %>%
    mutate(id = factor(id, levels = names(sort(w)))) %>%
    ggplot2::ggplot(aes(x=time, y=id, fill=factor(z))) +
    ggplot2::geom_tile(color="gray20") +
    viridis::scale_fill_viridis(discrete = T) +
    ggplot2::labs(x="Time", y="Id", fill=expression(Z[it])) +
    ggplot2::geom_hline(yintercept = which(unique(names(sort(w))) %in% cut_point) + 0.5, color = "red") +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(0)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(0))
}
