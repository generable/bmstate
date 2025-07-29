# Visualize a graph
transition_matrix_plot <- function(f, terminal_states, null_state,
                                   edge_labs, col_term = "firebrick1",
                                   col_null = "steelblue2", ...) {
  cn <- colnames(f)
  idx_term <- which(colnames(f) %in% terminal_states)
  idx_null <- which(colnames(f) %in% null_state)
  N <- length(cn)
  color <- rep("black", N)
  color[idx_term] <- col_term
  color[idx_null] <- col_null
  lcol <- rep("gray74", nrow(f))
  lcol[idx_term] <- col_term
  lcol[idx_null] <- col_null
  net <- network::network(f, directed = TRUE, loops = TRUE)

  # Create plot
  qgraph::qgraph(net,
    color = lcol,
    edge.labels = edge_labs,
    edge.color = "gray10",
    labels = cn,
    layout = "circle",
    ...
  )
}
