# ==============================================================================
# FINAL HSMM BRAIN-STATE TOPOLOGY ANALYSIS — V5 VISUAL THRESHOLD SENSITIVITY
# Reviewer-shareable reproducible script
#
# Purpose
#   1. Load the four state-specific correlation matrices from the final K = 4
#      HSMM.
#   2. Construct non-negative weighted undirected networks for topology metrics.
#   3. Calculate global efficiency, local efficiency, and consensus modularity.
#   4. Calculate additional descriptive graph metrics.
#   5. Visualize the strongest 15%, 20%, and 25% of positive edges within each state.
#
# Network definitions
#   - Nodes: 20 regions of interest (ROIs).
#   - Edge weights: state-specific Pearson correlations obtained from the fitted
#     Gaussian emission covariance matrices.
#   - Negative correlations: set to zero only for graph-theoretical analysis.
#   - Self-connections: excluded.
#
# Metric networks versus visualization networks
#   - Topological metrics are calculated from the complete positive weighted
#     network after negative edges are set to zero.
#   - Network figures retain the strongest 15%, 20%, and 25% of positive edges
#     separately within each state. These proportional thresholds give equal
#     graph density across states and are used only for visualization.
#   - The primary topology metrics remain calculated from the complete positive
#     weighted networks and are not recalculated after thresholding.
#
# Required package
#   igraph
#
# Required input
#   Preferred:
#     Updated_State_Results/correlation_matrix_state_1.csv
#     ...
#     Updated_State_Results/correlation_matrix_state_4.csv
#
#   Fallback:
#     FINAL_SELECTED_HSMM_K04.rds
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. User settings
# ------------------------------------------------------------------------------

# Run this script from the repository root. The root can alternatively be set
# through the HSMM_PROJECT_ROOT environment variable.
PROJECT_ROOT <- normalizePath(
  Sys.getenv("HSMM_PROJECT_ROOT", unset = "."),
  winslash = "/",
  mustWork = TRUE
)
PROJECT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "HSMM_reviewer_reproducible_analysis"
)

STATE_RESULTS_DIR <- file.path(
  PROJECT_DIR,
  "Updated_State_Results"
)

MODEL_FILE <- file.path(
  PROJECT_DIR,
  "FINAL_SELECTED_HSMM_K04.rds"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "Updated_Topology_Results"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

K <- 4L

# Repeated Louvain runs used to obtain a consensus partition.
N_LOUVAIN_REPETITIONS <- 150L
CONSENSUS_THRESHOLD <- 0.50
COMMUNITY_SEED <- 20260715L

# Proportional thresholds used only for network visualization and reviewer
# sensitivity analysis. Primary graph metrics remain based on the complete
# positive weighted networks.
VISUALIZATION_EDGE_PROPORTIONS <- c(
  0.15,
  0.20,
  0.25
)

# "independent_fr" reproduces a separate Fruchterman-Reingold layout for each
# state. "common_fr" uses one layout for all states and facilitates direct
# visual comparison.
FIGURE_LAYOUT_MODE <- "independent_fr"

LAYOUT_SEED <- 20260716L

CREATE_FIGURES <- TRUE

# For reproducibility and user control, packages are not installed automatically.
# Install igraph before running this script.
AUTO_INSTALL_IGRAPH <- FALSE


# ------------------------------------------------------------------------------
# 2. Install/load igraph
# ------------------------------------------------------------------------------

if (!requireNamespace("igraph", quietly = TRUE)) {
  if (!isTRUE(AUTO_INSTALL_IGRAPH)) {
    stop(
      "Package 'igraph' is required but is not installed.\n",
      "Install it using install.packages('igraph'), then rerun this script."
    )
  }

  install.packages(
    "igraph",
    dependencies = TRUE
  )
}

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop(
    "Package 'igraph' could not be loaded after installation."
  )
}

suppressPackageStartupMessages(
  library(igraph)
)


# ------------------------------------------------------------------------------
# 3. ROI names
# ------------------------------------------------------------------------------

ROI_NAMES <- c(
  "FC2",
  "FC3",
  "FC4",
  "FC6",
  "FC7",
  "FC8",
  "IFG1",
  "IFG4",
  "MC-brodman1-3 L",
  "MC-brodman1-3 R",
  "MC-brodman 4 R",
  "MC-brodman 5 L",
  "MC-brodman 5 R",
  "MC-brodman 6 L",
  "MC-brodman 6 R",
  "PC1",
  "PC2",
  "STG L",
  "THALAMUS L",
  "THALAMUS R"
)


# ------------------------------------------------------------------------------
# 4. Helper functions for loading and validating matrices
# ------------------------------------------------------------------------------

as_plain_symmetric_matrix <- function(
    x,
    matrix_name = "matrix"
) {
  x <- as.matrix(
    x
  )

  if (
    length(dim(x)) != 2L ||
    nrow(x) != ncol(x)
  ) {
    stop(
      matrix_name,
      " is not a square matrix."
    )
  }

  nr <- nrow(
    x
  )

  nc <- ncol(
    x
  )

  row_names <- rownames(
    x
  )

  column_names <- colnames(
    x
  )

  x <- matrix(
    as.numeric(x),
    nrow = nr,
    ncol = nc
  )

  if (any(!is.finite(x))) {
    stop(
      matrix_name,
      " contains non-finite values."
    )
  }

  x <- (
    x + t(x)
  ) / 2

  rownames(x) <- row_names
  colnames(x) <- column_names

  storage.mode(x) <- "double"

  x
}


covariance_to_correlation <- function(
    covariance_matrix
) {
  covariance_matrix <- as_plain_symmetric_matrix(
    covariance_matrix,
    matrix_name = "Covariance matrix"
  )

  variances <- diag(
    covariance_matrix
  )

  if (
    any(!is.finite(variances)) ||
    any(variances <= 0)
  ) {
    stop(
      "A fitted covariance matrix contains a non-positive variance."
    )
  }

  correlation_matrix <- stats::cov2cor(
    covariance_matrix
  )

  correlation_matrix <- matrix(
    as.numeric(correlation_matrix),
    nrow = nrow(covariance_matrix),
    ncol = ncol(covariance_matrix)
  )

  correlation_matrix[
    correlation_matrix > 1
  ] <- 1

  correlation_matrix[
    correlation_matrix < -1
  ] <- -1

  correlation_matrix <- (
    correlation_matrix +
      t(correlation_matrix)
  ) / 2

  diag(correlation_matrix) <- 1

  correlation_matrix
}


read_state_correlation_csv <- function(
    state
) {
  correlation_file <- file.path(
    STATE_RESULTS_DIR,
    sprintf(
      "correlation_matrix_state_%d.csv",
      state
    )
  )

  if (!file.exists(correlation_file)) {
    return(NULL)
  }

  matrix_data <- utils::read.csv(
    correlation_file,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  as_plain_symmetric_matrix(
    matrix_data,
    matrix_name = paste(
      "State",
      state,
      "correlation matrix"
    )
  )
}


# ------------------------------------------------------------------------------
# 5. Load final state-specific correlation matrices
# ------------------------------------------------------------------------------

correlation_matrices <- lapply(
  seq_len(K),
  read_state_correlation_csv
)

all_CSV_files_available <- all(
  vapply(
    correlation_matrices,
    function(x) !is.null(x),
    FUN.VALUE = logical(1)
  )
)

if (!all_CSV_files_available) {
  if (!file.exists(MODEL_FILE)) {
    stop(
      "One or more state correlation CSV files are missing and the final ",
      "fitted model was not found:\n",
      MODEL_FILE
    )
  }

  final_fit <- readRDS(
    MODEL_FILE
  )

  sigma_list <-
    final_fit$model$parms.emission$sigma

  if (
    is.null(sigma_list) ||
    length(sigma_list) != K
  ) {
    stop(
      "The fitted model does not contain four state covariance matrices."
    )
  }

  correlation_matrices <- lapply(
    sigma_list,
    covariance_to_correlation
  )

  matrix_source <- "FINAL_SELECTED_HSMM_K04.rds"
} else {
  matrix_source <- paste0(
    "correlation_matrix_state_1.csv to ",
    "correlation_matrix_state_4.csv"
  )
}

p <- nrow(
  correlation_matrices[[1]]
)

if (p != length(ROI_NAMES)) {
  stop(
    "The loaded matrices contain ",
    p,
    " ROIs, but ROI_NAMES contains ",
    length(ROI_NAMES),
    " labels."
  )
}

for (state in seq_len(K)) {
  matrix_state <- as_plain_symmetric_matrix(
    correlation_matrices[[state]],
    matrix_name = paste(
      "State",
      state,
      "correlation matrix"
    )
  )

  if (
    nrow(matrix_state) != p ||
    ncol(matrix_state) != p
  ) {
    stop(
      "Correlation matrices do not have consistent dimensions."
    )
  }

  diag(matrix_state) <- 1

  dimnames(matrix_state) <- list(
    ROI_NAMES,
    ROI_NAMES
  )

  correlation_matrices[[state]] <- matrix_state
}

names(correlation_matrices) <- paste0(
  "State_",
  seq_len(K)
)

# Familiar object names retained for compatibility with previous scripts.
Corr_list <- correlation_matrices

Cor1 <- Corr_list[[1]]
Cor2 <- Corr_list[[2]]
Cor3 <- Corr_list[[3]]
Cor4 <- Corr_list[[4]]


# ------------------------------------------------------------------------------
# 6. Construct complete positive weighted adjacency matrices
# ------------------------------------------------------------------------------

prepare_positive_adjacency <- function(
    correlation_matrix
) {
  adjacency <- as_plain_symmetric_matrix(
    correlation_matrix,
    matrix_name = "Correlation matrix"
  )

  adjacency[
    adjacency < 0
  ] <- 0

  diag(adjacency) <- 0

  adjacency <- (
    adjacency + t(adjacency)
  ) / 2

  dimnames(adjacency) <- list(
    ROI_NAMES,
    ROI_NAMES
  )

  adjacency
}


positive_adjacency_matrices <- lapply(
  Corr_list,
  prepare_positive_adjacency
)

names(positive_adjacency_matrices) <- names(
  Corr_list
)


create_weighted_graph <- function(
    adjacency
) {
  graph <- igraph::graph_from_adjacency_matrix(
    adjacency,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  igraph::V(graph)$name <- ROI_NAMES

  graph
}


full_positive_graphs <- lapply(
  positive_adjacency_matrices,
  create_weighted_graph
)

names(full_positive_graphs) <- names(
  positive_adjacency_matrices
)


# ------------------------------------------------------------------------------
# 7. Weighted shortest-path helper
# ------------------------------------------------------------------------------

weighted_distance_matrix <- function(
    graph
) {
  number_of_nodes <- igraph::vcount(
    graph
  )

  if (number_of_nodes == 0L) {
    return(
      matrix(
        numeric(0),
        0,
        0
      )
    )
  }

  if (igraph::ecount(graph) == 0L) {
    distance_matrix <- matrix(
      Inf,
      nrow = number_of_nodes,
      ncol = number_of_nodes
    )

    diag(distance_matrix) <- 0

    return(
      distance_matrix
    )
  }

  edge_weights <- igraph::E(graph)$weight

  if (
    any(!is.finite(edge_weights)) ||
    any(edge_weights <= 0)
  ) {
    stop(
      "All graph edge weights must be finite and strictly positive."
    )
  }

  igraph::distances(
    graph,
    mode = "all",
    weights = 1 / edge_weights,
    algorithm = "automatic"
  )
}


# ------------------------------------------------------------------------------
# 8. Global efficiency
# ------------------------------------------------------------------------------

global_efficiency_weighted <- function(
    graph
) {
  number_of_nodes <- igraph::vcount(
    graph
  )

  if (number_of_nodes < 2L) {
    return(NA_real_)
  }

  distance_matrix <- weighted_distance_matrix(
    graph
  )

  inverse_distances <- 1 / distance_matrix

  inverse_distances[
    !is.finite(inverse_distances)
  ] <- 0

  diag(inverse_distances) <- 0

  sum(inverse_distances) / (
    number_of_nodes *
      (number_of_nodes - 1)
  )
}


# ------------------------------------------------------------------------------
# 9. Local efficiency
#
# This implements the article-aligned weighted neighborhood formulation used
# in the previous analysis. For each node, the induced graph among its
# neighbors is extracted; subgraph weights are cube-root transformed before
# weighted shortest-path efficiency is calculated.
# ------------------------------------------------------------------------------

local_efficiency_weighted <- function(
    graph
) {
  number_of_nodes <- igraph::vcount(
    graph
  )

  if (number_of_nodes < 3L) {
    return(NA_real_)
  }

  node_local_efficiency <- rep(
    NA_real_,
    number_of_nodes
  )

  for (node in seq_len(number_of_nodes)) {
    neighbor_vertices <- igraph::neighbors(
      graph,
      node,
      mode = "all"
    )

    number_of_neighbors <- length(
      neighbor_vertices
    )

    if (number_of_neighbors < 2L) {
      next
    }

    neighbor_subgraph <- igraph::induced_subgraph(
      graph,
      vids = neighbor_vertices
    )

    if (igraph::ecount(neighbor_subgraph) == 0L) {
      node_local_efficiency[node] <- 0
      next
    }

    igraph::E(neighbor_subgraph)$weight <-
      igraph::E(neighbor_subgraph)$weight ^ (
        1 / 3
      )

    distance_matrix <- weighted_distance_matrix(
      neighbor_subgraph
    )

    inverse_distances <- 1 / distance_matrix

    inverse_distances[
      !is.finite(inverse_distances)
    ] <- 0

    diag(inverse_distances) <- 0

    node_local_efficiency[node] <-
      sum(inverse_distances) / (
        number_of_neighbors *
          (number_of_neighbors - 1)
      )
  }

  valid_values <- node_local_efficiency[
    is.finite(node_local_efficiency)
  ]

  if (length(valid_values) == 0L) {
    return(NA_real_)
  }

  mean(
    valid_values
  )
}


# ------------------------------------------------------------------------------
# 10. Consensus modularity
# ------------------------------------------------------------------------------

consensus_modularity <- function(
    graph,
    repetitions = 150L,
    agreement_threshold = 0.50,
    seed = 1L
) {
  number_of_nodes <- igraph::vcount(
    graph
  )

  if (
    number_of_nodes < 2L ||
    igraph::ecount(graph) == 0L
  ) {
    return(
      list(
        modularity = NA_real_,
        membership = rep(
          NA_integer_,
          number_of_nodes
        ),
        agreement_matrix = matrix(
          NA_real_,
          number_of_nodes,
          number_of_nodes
        ),
        number_of_communities = NA_integer_
      )
    )
  }

  set.seed(
    seed
  )

  membership_matrix <- matrix(
    NA_integer_,
    nrow = number_of_nodes,
    ncol = repetitions
  )

  for (repetition in seq_len(repetitions)) {
    community_result <- igraph::cluster_louvain(
      graph,
      weights = igraph::E(graph)$weight
    )

    membership_matrix[
      ,
      repetition
    ] <- igraph::membership(
      community_result
    )
  }

  agreement_matrix <- matrix(
    0,
    nrow = number_of_nodes,
    ncol = number_of_nodes
  )

  for (repetition in seq_len(repetitions)) {
    one_membership <- membership_matrix[
      ,
      repetition
    ]

    agreement_matrix <- agreement_matrix +
      outer(
        one_membership,
        one_membership,
        FUN = "=="
      )
  }

  agreement_matrix <- agreement_matrix /
    repetitions

  diag(agreement_matrix) <- 0

  agreement_matrix[
    agreement_matrix < agreement_threshold
  ] <- 0

  consensus_graph <-
    igraph::graph_from_adjacency_matrix(
      agreement_matrix,
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )

  if (igraph::ecount(consensus_graph) == 0L) {
    final_membership <- seq_len(
      number_of_nodes
    )

    final_modularity <- 0
  } else {
    final_community <- igraph::cluster_louvain(
      consensus_graph,
      weights = igraph::E(
        consensus_graph
      )$weight
    )

    final_membership <- igraph::membership(
      final_community
    )

    # Article-aligned value: modularity of the final consensus graph.
    final_modularity <- igraph::modularity(
      final_community
    )
  }

  list(
    modularity = as.numeric(
      final_modularity
    ),
    membership = as.integer(
      final_membership
    ),
    agreement_matrix = agreement_matrix,
    number_of_communities = length(
      unique(final_membership)
    )
  )
}


# ------------------------------------------------------------------------------
# 11. Additional graph metrics
# ------------------------------------------------------------------------------

average_path_distance_weighted <- function(
    graph
) {
  number_of_nodes <- igraph::vcount(
    graph
  )

  if (number_of_nodes < 2L) {
    return(NA_real_)
  }

  distance_matrix <- weighted_distance_matrix(
    graph
  )

  diag(distance_matrix) <- NA_real_

  finite_distances <- distance_matrix[
    is.finite(distance_matrix)
  ]

  if (length(finite_distances) == 0L) {
    return(NA_real_)
  }

  mean(
    finite_distances
  )
}


weighted_clustering_coefficient <- function(
    graph
) {
  if (igraph::vcount(graph) < 3L) {
    return(NA_real_)
  }

  local_clustering <- igraph::transitivity(
    graph,
    type = "weighted",
    isolates = "zero"
  )

  local_clustering <- local_clustering[
    is.finite(local_clustering)
  ]

  if (length(local_clustering) == 0L) {
    return(NA_real_)
  }

  mean(
    local_clustering
  )
}


mean_weighted_betweenness <- function(
    graph
) {
  if (
    igraph::vcount(graph) < 3L ||
    igraph::ecount(graph) == 0L
  ) {
    return(NA_real_)
  }

  betweenness_values <- igraph::betweenness(
    graph,
    directed = FALSE,
    weights = 1 / igraph::E(graph)$weight,
    normalized = TRUE
  )

  mean(
    betweenness_values
  )
}


# ------------------------------------------------------------------------------
# 12. Calculate state-level topology
# ------------------------------------------------------------------------------

GlobalEff <- rep(
  NA_real_,
  K
)

LocalEff <- rep(
  NA_real_,
  K
)

Modularity <- rep(
  NA_real_,
  K
)

NumberCommunities <- rep(
  NA_integer_,
  K
)

MeanDegree <- rep(
  NA_real_,
  K
)

MeanStrength <- rep(
  NA_real_,
  K
)

Cost <- rep(
  NA_real_,
  K
)

AveragePathLength <- rep(
  NA_real_,
  K
)

ClusteringCoefficient <- rep(
  NA_real_,
  K
)

MeanBetweenness <- rep(
  NA_real_,
  K
)

NumberPositiveEdges <- integer(
  K
)

NetworkDensity <- rep(
  NA_real_,
  K
)

consensus_results <- vector(
  mode = "list",
  length = K
)

node_metric_rows <- list()
node_metric_counter <- 1L

for (state in seq_len(K)) {
  graph <- full_positive_graphs[[state]]

  GlobalEff[state] <-
    global_efficiency_weighted(
      graph
    )

  LocalEff[state] <-
    local_efficiency_weighted(
      graph
    )

  consensus_results[[state]] <-
    consensus_modularity(
      graph = graph,
      repetitions =
        N_LOUVAIN_REPETITIONS,
      agreement_threshold =
        CONSENSUS_THRESHOLD,
      seed =
        COMMUNITY_SEED + state
    )

  Modularity[state] <-
    consensus_results[[state]]$modularity

  NumberCommunities[state] <-
    consensus_results[[state]]$number_of_communities

  node_degree <- igraph::degree(
    graph,
    mode = "all",
    loops = FALSE
  )

  node_strength <- igraph::strength(
    graph,
    mode = "all",
    loops = FALSE,
    weights = igraph::E(graph)$weight
  )

  node_betweenness <- if (
    igraph::ecount(graph) > 0L
  ) {
    igraph::betweenness(
      graph,
      directed = FALSE,
      weights =
        1 / igraph::E(graph)$weight,
      normalized = TRUE
    )
  } else {
    rep(
      NA_real_,
      p
    )
  }

  node_clustering <- igraph::transitivity(
    graph,
    type = "weighted",
    isolates = "zero"
  )

  MeanDegree[state] <- mean(
    node_degree
  )

  MeanStrength[state] <- mean(
    node_strength
  )

  Cost[state] <- MeanDegree[state] / (
    p - 1
  )

  AveragePathLength[state] <-
    average_path_distance_weighted(
      graph
    )

  ClusteringCoefficient[state] <-
    weighted_clustering_coefficient(
      graph
    )

  MeanBetweenness[state] <-
    mean_weighted_betweenness(
      graph
    )

  NumberPositiveEdges[state] <-
    igraph::ecount(
      graph
    )

  NetworkDensity[state] <-
    igraph::edge_density(
      graph,
      loops = FALSE
    )

  for (node in seq_len(p)) {
    node_metric_rows[[node_metric_counter]] <- data.frame(
      state = state,
      ROI_number = node,
      ROI = ROI_NAMES[node],
      degree = as.numeric(
        node_degree[node]
      ),
      strength = as.numeric(
        node_strength[node]
      ),
      weighted_clustering = as.numeric(
        node_clustering[node]
      ),
      normalized_betweenness = as.numeric(
        node_betweenness[node]
      ),
      consensus_community =
        consensus_results[[state]]$membership[node],
      stringsAsFactors = FALSE
    )

    node_metric_counter <-
      node_metric_counter + 1L
  }
}


Topology_Table <- data.frame(
  State = paste0(
    "State ",
    seq_len(K)
  ),
  Global_Efficiency = GlobalEff,
  Local_Efficiency = LocalEff,
  Modularity = Modularity,
  Number_of_Communities =
    NumberCommunities,
  stringsAsFactors = FALSE
)

Additional_Graph_Metrics_Table <- data.frame(
  State = paste0(
    "State ",
    seq_len(K)
  ),
  Positive_Edges = NumberPositiveEdges,
  Density = NetworkDensity,
  Mean_Degree = MeanDegree,
  Cost = Cost,
  Mean_Strength = MeanStrength,
  Average_Path_Length =
    AveragePathLength,
  Weighted_Clustering_Coefficient =
    ClusteringCoefficient,
  Mean_Normalized_Betweenness =
    MeanBetweenness,
  stringsAsFactors = FALSE
)

Node_Metrics_Table <- do.call(
  rbind,
  node_metric_rows
)

utils::write.csv(
  Topology_Table,
  file.path(
    OUTPUT_DIR,
    "state_topology_primary_metrics.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  Additional_Graph_Metrics_Table,
  file.path(
    OUTPUT_DIR,
    "state_topology_additional_metrics.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  Node_Metrics_Table,
  file.path(
    OUTPUT_DIR,
    "state_topology_node_metrics.csv"
  ),
  row.names = FALSE
)


# Save consensus partitions and agreement matrices.
for (state in seq_len(K)) {
  membership_table <- data.frame(
    ROI_number = seq_len(p),
    ROI = ROI_NAMES,
    consensus_community =
      consensus_results[[state]]$membership,
    stringsAsFactors = FALSE
  )

  utils::write.csv(
    membership_table,
    file.path(
      OUTPUT_DIR,
      sprintf(
        "State_%d_consensus_communities.csv",
        state
      )
    ),
    row.names = FALSE
  )

  agreement_matrix <-
    consensus_results[[state]]$agreement_matrix

  dimnames(agreement_matrix) <- list(
    ROI_NAMES,
    ROI_NAMES
  )

  utils::write.csv(
    agreement_matrix,
    file.path(
      OUTPUT_DIR,
      sprintf(
        "State_%d_consensus_agreement_matrix.csv",
        state
      )
    ),
    row.names = TRUE
  )
}


# ------------------------------------------------------------------------------
# 13. Visualization-only threshold sensitivity: 15%, 20%, and 25%
#
# This section repeats the same graph-visualization procedure at three
# proportional densities. It does NOT recalculate global efficiency, local
# efficiency, modularity, or the additional topology metrics.
# ------------------------------------------------------------------------------

threshold_to_fixed_density <- function(
    positive_adjacency,
    proportion
) {
  if (
    length(proportion) != 1L ||
    !is.finite(proportion) ||
    proportion <= 0 ||
    proportion > 1
  ) {
    stop(
      "Each retained-edge proportion must be a single value in (0, 1]."
    )
  }

  number_of_nodes <- nrow(
    positive_adjacency
  )

  total_possible_edges <- (
    number_of_nodes *
      (number_of_nodes - 1)
  ) / 2

  # An integer number of edges is required.
  # With 20 ROIs:
  #   15% of 190 = 28.5 -> 29 edges
  #   20% of 190 = 38 edges
  #   25% of 190 = 47.5 -> 48 edges
  target_number_of_edges <- max(
    1L,
    as.integer(
      ceiling(
        proportion *
          total_possible_edges
      )
    )
  )

  upper_indices <- which(
    upper.tri(
      positive_adjacency,
      diag = FALSE
    ),
    arr.ind = TRUE
  )

  upper_weights <- positive_adjacency[
    upper.tri(
      positive_adjacency,
      diag = FALSE
    )
  ]

  positive_positions <- which(
    is.finite(upper_weights) &
      upper_weights > 0
  )

  thresholded <- matrix(
    0,
    nrow = number_of_nodes,
    ncol = number_of_nodes,
    dimnames = dimnames(
      positive_adjacency
    )
  )

  if (length(positive_positions) == 0L) {
    return(
      list(
        adjacency = thresholded,
        threshold = NA_real_,
        retained_edges = 0L,
        target_edges =
          target_number_of_edges,
        available_positive_edges = 0L,
        total_possible_edges =
          as.integer(total_possible_edges),
        requested_proportion =
          proportion,
        achieved_density = 0
      )
    )
  }

  if (
    length(positive_positions) <
      target_number_of_edges
  ) {
    warning(
      "A state contains fewer positive edges than the requested fixed-density ",
      "target. All available positive edges will be retained."
    )
  }

  number_to_retain <- min(
    target_number_of_edges,
    length(positive_positions)
  )

  ranked_positions <- positive_positions[
    order(
      upper_weights[
        positive_positions
      ],
      decreasing = TRUE,
      na.last = NA
    )
  ]

  retained_positions <- ranked_positions[
    seq_len(
      number_to_retain
    )
  ]

  for (position in retained_positions) {
    row_index <- upper_indices[
      position,
      1
    ]

    column_index <- upper_indices[
      position,
      2
    ]

    weight <- upper_weights[
      position
    ]

    thresholded[
      row_index,
      column_index
    ] <- weight

    thresholded[
      column_index,
      row_index
    ] <- weight
  }

  retained_threshold <- min(
    upper_weights[
      retained_positions
    ]
  )

  list(
    adjacency = thresholded,
    threshold = retained_threshold,
    retained_edges =
      as.integer(number_to_retain),
    target_edges =
      as.integer(target_number_of_edges),
    available_positive_edges =
      as.integer(length(positive_positions)),
    total_possible_edges =
      as.integer(total_possible_edges),
    requested_proportion =
      proportion,
    achieved_density =
      number_to_retain /
        total_possible_edges
  )
}


threshold_percentages <- as.integer(
  round(
    100 *
      VISUALIZATION_EDGE_PROPORTIONS
  )
)

threshold_labels <- as.character(
  threshold_percentages
)

names(VISUALIZATION_EDGE_PROPORTIONS) <-
  threshold_labels


threshold_results_by_level <- vector(
  mode = "list",
  length = length(
    VISUALIZATION_EDGE_PROPORTIONS
  )
)

names(threshold_results_by_level) <-
  threshold_labels

thresholded_adjacency_by_level <-
  threshold_results_by_level

visualization_graphs_by_level <-
  threshold_results_by_level


for (
  threshold_index in
  seq_along(
    VISUALIZATION_EDGE_PROPORTIONS
  )
) {
  threshold_label <-
    threshold_labels[
      threshold_index
    ]

  threshold_proportion <-
    VISUALIZATION_EDGE_PROPORTIONS[
      threshold_index
    ]

  threshold_results_by_level[[threshold_label]] <-
    lapply(
      positive_adjacency_matrices,
      threshold_to_fixed_density,
      proportion =
        threshold_proportion
    )

  thresholded_adjacency_by_level[[threshold_label]] <-
    lapply(
      threshold_results_by_level[[threshold_label]],
      function(result) {
        result$adjacency
      }
    )

  visualization_graphs_by_level[[threshold_label]] <-
    lapply(
      thresholded_adjacency_by_level[[threshold_label]],
      create_weighted_graph
    )
}


# ------------------------------------------------------------------------------
# 14. Save threshold summaries and adjacency matrices
# ------------------------------------------------------------------------------

threshold_summary_rows <- list()
threshold_summary_counter <- 1L

for (
  threshold_index in
  seq_along(
    VISUALIZATION_EDGE_PROPORTIONS
  )
) {
  threshold_label <-
    threshold_labels[
      threshold_index
    ]

  for (state in seq_len(K)) {
    one_result <-
      threshold_results_by_level[[threshold_label]][[state]]

    threshold_summary_rows[[threshold_summary_counter]] <-
      data.frame(
        Threshold_Percent =
          threshold_percentages[
            threshold_index
          ],
        Threshold_Proportion =
          VISUALIZATION_EDGE_PROPORTIONS[
            threshold_index
          ],
        State = paste0(
          "State ",
          state
        ),
        Total_Possible_Edges =
          one_result$total_possible_edges,
        Positive_Edges_Available =
          one_result$available_positive_edges,
        Target_Edges =
          one_result$target_edges,
        Edges_Retained =
          one_result$retained_edges,
        Minimum_Retained_Correlation =
          one_result$threshold,
        Achieved_Density =
          one_result$achieved_density,
        stringsAsFactors = FALSE
      )

    threshold_summary_counter <-
      threshold_summary_counter + 1L

    utils::write.csv(
      thresholded_adjacency_by_level[[threshold_label]][[state]],
      file.path(
        OUTPUT_DIR,
        sprintf(
          "State_%d_top%dpercent_adjacency.csv",
          state,
          threshold_percentages[
            threshold_index
          ]
        )
      ),
      row.names = TRUE
    )
  }
}

threshold_summary <- do.call(
  rbind,
  threshold_summary_rows
)

utils::write.csv(
  threshold_summary,
  file.path(
    OUTPUT_DIR,
    "network_visualization_threshold_sensitivity_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 15. Figure helpers
# ------------------------------------------------------------------------------

rescale_numeric <- function(
    x,
    to = c(
      0,
      1
    )
) {
  x <- as.numeric(
    x
  )

  if (length(x) == 0L) {
    return(
      numeric(0)
    )
  }

  finite_x <- x[
    is.finite(x)
  ]

  if (length(finite_x) == 0L) {
    return(
      rep(
        mean(to),
        length(x)
      )
    )
  }

  x_min <- min(
    finite_x
  )

  x_max <- max(
    finite_x
  )

  if (x_max == x_min) {
    return(
      rep(
        mean(to),
        length(x)
      )
    )
  }

  to[1] + (
    x - x_min
  ) * (
    to[2] - to[1]
  ) / (
    x_max - x_min
  )
}


create_common_layout <- function(
    adjacency_list,
    seed
) {
  mean_adjacency <- Reduce(
    "+",
    adjacency_list
  ) / length(
    adjacency_list
  )

  union_graph <- create_weighted_graph(
    mean_adjacency
  )

  set.seed(
    seed
  )

  igraph::layout_with_fr(
    union_graph,
    weights = igraph::E(
      union_graph
    )$weight,
    niter = 1000
  )
}


if (
  !FIGURE_LAYOUT_MODE %in%
    c(
      "independent_fr",
      "common_fr"
    )
) {
  stop(
    "FIGURE_LAYOUT_MODE must be 'independent_fr' or 'common_fr'."
  )
}


plot_state_network <- function(
    graph,
    state,
    layout_coordinates = NULL,
    layout_seed
) {
  if (is.null(layout_coordinates)) {
    set.seed(
      layout_seed + state
    )

    layout_coordinates <- igraph::layout_with_fr(
      graph,
      weights = if (
        igraph::ecount(graph) > 0L
      ) {
        igraph::E(graph)$weight
      } else {
        NULL
      },
      niter = 1000
    )
  }

  if (igraph::ecount(graph) > 0L) {
    edge_weights <- igraph::E(
      graph
    )$weight

    edge_widths <- rescale_numeric(
      edge_weights,
      to = c(
        1.5,
        6
      )
    )

    edge_alpha <- rescale_numeric(
      edge_weights,
      to = c(
        0.35,
        0.95
      )
    )

    base_edge_RGB <- grDevices::col2rgb(
      "#B2182B"
    ) / 255

    edge_colours <- grDevices::rgb(
      red = rep(
        base_edge_RGB[1, 1],
        length(edge_alpha)
      ),
      green = rep(
        base_edge_RGB[2, 1],
        length(edge_alpha)
      ),
      blue = rep(
        base_edge_RGB[3, 1],
        length(edge_alpha)
      ),
      alpha = edge_alpha
    )
  } else {
    edge_widths <- NULL
    edge_colours <- NULL
  }

  plot(
    graph,
    layout = layout_coordinates,
    vertex.size = 18,
    vertex.color = "#FFF2CC",
    vertex.frame.color = "#4D4D4D",
    vertex.label = seq_len(p),
    vertex.label.cex = 0.80,
    vertex.label.font = 2,
    vertex.label.color = "black",
    edge.width = edge_widths,
    edge.color = edge_colours,
    main = paste(
      "State",
      state
    ),
    main.font = 2,
    asp = 1,
    margin = 0.10
  )
}


draw_network_figure <- function(
    graph_list,
    adjacency_list,
    threshold_percent,
    seed
) {
  graphics::layout(
    matrix(
      c(
        1,
        2,
        5,
        3,
        4,
        5
      ),
      nrow = 2,
      byrow = TRUE
    ),
    widths = c(
      1,
      1,
      1.10
    )
  )

  graphics::par(
    mar = c(
      1,
      1,
      3,
      1
    ),
    oma = c(
      0,
      0,
      2,
      0
    )
  )

  common_layout <- if (
    FIGURE_LAYOUT_MODE == "common_fr"
  ) {
    create_common_layout(
      adjacency_list,
      seed = seed
    )
  } else {
    NULL
  }

  for (state in seq_len(K)) {
    plot_state_network(
      graph =
        graph_list[[state]],
      state = state,
      layout_coordinates =
        common_layout,
      layout_seed =
        seed
    )
  }

  graphics::par(
    mar = c(
      1,
      1,
      1,
      1
    )
  )

  graphics::plot.new()

  graphics::legend(
    "center",
    legend = paste0(
      seq_len(p),
      "  ",
      ROI_NAMES
    ),
    title =
      "Regions of Interest (ROIs)",
    title.font = 2,
    cex = 0.82,
    text.font = 1,
    bty = "n",
    x.intersp = 0.5,
    y.intersp = 1.0
  )

  graphics::mtext(
    paste0(
      "Strongest ",
      threshold_percent,
      "% of positive connections"
    ),
    side = 3,
    outer = TRUE,
    line = 0.4,
    font = 2,
    cex = 1.1
  )
}


# ------------------------------------------------------------------------------
# 16. Save network figures at 15%, 20%, and 25%
# ------------------------------------------------------------------------------

if (isTRUE(CREATE_FIGURES)) {
  for (
    threshold_index in
    seq_along(
      VISUALIZATION_EDGE_PROPORTIONS
    )
  ) {
    threshold_label <-
      threshold_labels[
        threshold_index
      ]

    threshold_percent <-
      threshold_percentages[
        threshold_index
      ]

    graph_list <-
      visualization_graphs_by_level[[threshold_label]]

    adjacency_list <-
      thresholded_adjacency_by_level[[threshold_label]]

    figure_seed <-
      LAYOUT_SEED +
        1000L *
          threshold_index

    grDevices::png(
      filename = file.path(
        OUTPUT_DIR,
        sprintf(
          "State_Networks_Top%dPercent.png",
          threshold_percent
        )
      ),
      width = 4200,
      height = 2800,
      res = 300
    )

    draw_network_figure(
      graph_list =
        graph_list,
      adjacency_list =
        adjacency_list,
      threshold_percent =
        threshold_percent,
      seed =
        figure_seed
    )

    grDevices::dev.off()


    grDevices::pdf(
      file = file.path(
        OUTPUT_DIR,
        sprintf(
          "State_Networks_Top%dPercent.pdf",
          threshold_percent
        )
      ),
      width = 14,
      height = 9,
      onefile = TRUE
    )

    draw_network_figure(
      graph_list =
        graph_list,
      adjacency_list =
        adjacency_list,
      threshold_percent =
        threshold_percent,
      seed =
        figure_seed
    )

    grDevices::dev.off()
  }
}


# ------------------------------------------------------------------------------
# 17. Rounded console tables
# ------------------------------------------------------------------------------

Topology_Table_Rounded <- Topology_Table

Topology_Table_Rounded[
  c(
    "Global_Efficiency",
    "Local_Efficiency",
    "Modularity"
  )
] <- lapply(
  Topology_Table_Rounded[
    c(
      "Global_Efficiency",
      "Local_Efficiency",
      "Modularity"
    )
  ],
  function(x) {
    round(
      x,
      6
    )
  }
)


Additional_Graph_Metrics_Table_Rounded <-
  Additional_Graph_Metrics_Table

additional_numeric_columns <- setdiff(
  names(
    Additional_Graph_Metrics_Table_Rounded
  ),
  c(
    "State",
    "Positive_Edges"
  )
)

Additional_Graph_Metrics_Table_Rounded[
  additional_numeric_columns
] <- lapply(
  Additional_Graph_Metrics_Table_Rounded[
    additional_numeric_columns
  ],
  function(x) {
    round(
      x,
      6
    )
  }
)


Threshold_Summary_Rounded <-
  threshold_summary

threshold_summary_numeric_columns <- c(
  "Threshold_Proportion",
  "Minimum_Retained_Correlation",
  "Achieved_Density"
)

Threshold_Summary_Rounded[
  threshold_summary_numeric_columns
] <- lapply(
  Threshold_Summary_Rounded[
    threshold_summary_numeric_columns
  ],
  function(x) {
    round(
      x,
      6
    )
  }
)


# ------------------------------------------------------------------------------
# 18. Reproducibility information
# ------------------------------------------------------------------------------

analysis_settings <- data.frame(
  setting = c(
    "matrix_source",
    "selected_K",
    "number_of_ROIs",
    "metric_network",
    "visualization_thresholds",
    "threshold_edge_count_denominator",
    "threshold_use",
    "negative_edge_handling",
    "self_connections",
    "global_efficiency_definition",
    "local_efficiency_definition",
    "community_algorithm",
    "Louvain_repetitions",
    "consensus_agreement_threshold",
    "community_seed",
    "figure_layout_mode",
    "layout_seed"
  ),
  value = c(
    matrix_source,
    K,
    p,
    "Complete positive weighted undirected network",
    paste(
      paste0(
        threshold_percentages,
        "%"
      ),
      collapse = ", "
    ),
    paste0(
      "All 190 possible undirected ROI pairs; identical edge count ",
      "across states at each threshold"
    ),
    paste0(
      "Visualization only; graph metrics are not recomputed after ",
      "thresholding"
    ),
    "Negative correlations set to zero for graph analysis only",
    "Excluded",
    paste0(
      "Mean inverse weighted shortest-path length; ",
      "edge length = 1 / correlation"
    ),
    paste0(
      "Mean neighborhood efficiency after cube-root transformation ",
      "of subgraph weights"
    ),
    "Repeated Louvain followed by thresholded agreement consensus",
    N_LOUVAIN_REPETITIONS,
    CONSENSUS_THRESHOLD,
    COMMUNITY_SEED,
    FIGURE_LAYOUT_MODE,
    LAYOUT_SEED
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  analysis_settings,
  file.path(
    OUTPUT_DIR,
    "topology_analysis_settings.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "topology_analysis_sessionInfo.txt"
  )
)


# ------------------------------------------------------------------------------
# 19. Final console output
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL BRAIN-STATE TOPOLOGY ANALYSIS COMPLETED\n")
cat("VISUAL SENSITIVITY AT 15%, 20%, AND 25%\n")
cat("============================================================\n")
cat("Matrix source:", matrix_source, "\n")
cat("States:", K, "\n")
cat("ROIs:", p, "\n")
cat(
  "Visualization thresholds:",
  paste(
    paste0(
      threshold_percentages,
      "%"
    ),
    collapse = ", "
  ),
  "\n"
)
cat(
  "Results directory:\n",
  normalizePath(
    OUTPUT_DIR,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
cat("============================================================\n\n")

cat("Primary topology metrics from the complete positive networks:\n")
print(
  Topology_Table_Rounded,
  row.names = FALSE
)

cat("\nVisualization threshold and density summary:\n")
print(
  Threshold_Summary_Rounded,
  row.names = FALSE
)

cat(
  "\nNote: 15%, 20%, and 25% thresholds were used only to create ",
  "the network figures. Primary topology metrics were not recalculated ",
  "on thresholded networks.\n",
  sep = ""
)
