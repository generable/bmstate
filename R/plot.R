# Visualize a graph
transition_matrix_plot <- function(f, terminal_states, null_state,
                                   edge_labs, col_term = "firebrick",
                                   col_null = "steelblue2", ...) {
  cn <- colnames(f)
  idx_term <- which(colnames(f) %in% terminal_states)
  idx_null <- which(colnames(f) %in% null_state)
  N <- length(cn)
  color <- rep("black", N)
  color[idx_term] <- col_term
  color[idx_null] <- col_null
  acol <- matrix("black", nrow(f), ncol(f))
  lcol <- acol
  lcol[, idx_term] <- col_term
  lcol[, idx_null] <- col_null

  # Create plot
  qgraph::qgraph(f,
    edge.labels = edge_labs, label.color = color,
    edge.color = acol,
    fade = FALSE,
    layout = "circle",
    ...
  )
}
