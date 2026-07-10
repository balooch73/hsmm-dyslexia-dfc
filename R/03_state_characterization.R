# ==============================================================================
# HSMM STATE-SPECIFIC FUNCTIONAL CHARACTERIZATION
# Reviewer-shareable reproducible script
#
# Purpose
#   1. Load the final selected K = 4 HSMM.
#   2. Extract state-specific Gaussian covariance and mean parameters.
#   3. Convert covariance matrices to Pearson correlation matrices.
#   4. Plot all four state-specific correlation matrices using one common scale.
#   5. Summarize connectivity and mean-signal characteristics.
#
# Required input
#   FINAL_SELECTED_HSMM_K04.rds
#
# Required packages
#   None beyond base R.
#
# Important
#   Run this script from the repository root, or set HSMM_PROJECT_ROOT.
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

MODEL_FILE <- file.path(
  PROJECT_DIR,
  "FINAL_SELECTED_HSMM_K04.rds"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "Updated_State_Results"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ROI order used during model estimation.
# The order must match the columns of the input time-series data.
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
# 2. Load and validate the final fitted model
# ------------------------------------------------------------------------------

if (!file.exists(MODEL_FILE)) {
  stop(
    "The fitted HSMM file was not found:\n",
    MODEL_FILE
  )
}

final_fit <- readRDS(MODEL_FILE)

emission_parameters <- final_fit$model$parms.emission

if (
  is.null(emission_parameters) ||
  is.null(emission_parameters$sigma) ||
  is.null(emission_parameters$mu)
) {
  stop(
    "The RDS file does not contain the expected HSMM emission parameters."
  )
}

sigma_list <- emission_parameters$sigma
mu_list <- emission_parameters$mu

K <- length(sigma_list)

if (K != 4L) {
  stop(
    "The loaded model contains ",
    K,
    " states; this script expects the selected K = 4 model."
  )
}

first_sigma <- as.matrix(sigma_list[[1]])

if (length(dim(first_sigma)) != 2L || nrow(first_sigma) != ncol(first_sigma)) {
  stop("The state covariance matrices are not square matrices.")
}

p <- nrow(first_sigma)

if (length(mu_list) != K) {
  stop("The numbers of state covariance matrices and state mean vectors differ.")
}

if (length(ROI_NAMES) != p) {
  stored_names <- colnames(first_sigma)

  if (!is.null(stored_names) && length(stored_names) == p) {
    roi_names <- stored_names
    warning(
      "ROI_NAMES did not match the model dimension; column names stored in ",
      "the fitted model were used."
    )
  } else {
    stop(
      "ROI_NAMES contains ",
      length(ROI_NAMES),
      " labels, whereas the model contains ",
      p,
      " variables. Check the ROI order before continuing."
    )
  }
} else {
  roi_names <- ROI_NAMES
}

cat("Final HSMM loaded successfully.\n")
cat("Number of states:", K, "\n")
cat("Number of ROIs:", p, "\n")


# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------

as_plain_symmetric_matrix <- function(x, matrix_name = "matrix") {
  x <- as.matrix(x)

  if (length(dim(x)) != 2L || nrow(x) != ncol(x)) {
    stop(matrix_name, " is not a square matrix.")
  }

  nr <- nrow(x)
  nc <- ncol(x)

  # Remove any additional class attributes, including class "shrinkage".
  x <- matrix(
    as.numeric(x),
    nrow = nr,
    ncol = nc
  )

  if (any(!is.finite(x))) {
    stop(matrix_name, " contains non-finite values.")
  }

  x <- (x + t(x)) / 2
  storage.mode(x) <- "double"

  x
}


covariance_to_correlation <- function(S) {
  S <- as_plain_symmetric_matrix(
    S,
    matrix_name = "Covariance matrix"
  )

  variances <- diag(S)

  if (any(!is.finite(variances)) || any(variances <= 0)) {
    stop(
      "A covariance matrix contains a non-positive or non-finite variance."
    )
  }

  R <- stats::cov2cor(S)

  # Explicit conversion preserves the matrix dimensions and removes
  # non-standard classes.
  R <- matrix(
    as.numeric(R),
    nrow = nrow(S),
    ncol = ncol(S)
  )

  # Correct only negligible floating-point excursions outside [-1, 1].
  R[R > 1] <- 1
  R[R < -1] <- -1

  R <- (R + t(R)) / 2
  diag(R) <- 1
  storage.mode(R) <- "double"

  R
}


get_unique_correlations <- function(R) {
  R <- as.matrix(R)
  R[upper.tri(R, diag = FALSE)]
}


safe_mean <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(NA_real_)
  }

  mean(x)
}


mean_positive <- function(x) {
  safe_mean(x[x > 0])
}


mean_negative_signed <- function(x) {
  safe_mean(x[x < 0])
}


mean_negative_magnitude <- function(x) {
  values <- x[x < 0]

  if (length(values) == 0L) {
    return(NA_real_)
  }

  safe_mean(abs(values))
}


# ------------------------------------------------------------------------------
# 4. Extract state-specific covariance and correlation matrices
# ------------------------------------------------------------------------------

covariance_matrices <- lapply(
  seq_len(K),
  function(state) {
    S <- as_plain_symmetric_matrix(
      sigma_list[[state]],
      matrix_name = paste("State", state, "covariance matrix")
    )

    if (nrow(S) != p) {
      stop("Unexpected covariance-matrix dimension in state ", state, ".")
    }

    dimnames(S) <- list(
      roi_names,
      roi_names
    )

    S
  }
)

correlation_matrices <- lapply(
  covariance_matrices,
  function(S) {
    R <- covariance_to_correlation(S)

    dimnames(R) <- list(
      roi_names,
      roi_names
    )

    R
  }
)

names(covariance_matrices) <- paste0(
  "State_",
  seq_len(K)
)

names(correlation_matrices) <- paste0(
  "State_",
  seq_len(K)
)

stopifnot(
  all(
    vapply(
      correlation_matrices,
      is.matrix,
      FUN.VALUE = logical(1)
    )
  )
)

# Familiar object names retained for transparent inspection.
cov1 <- covariance_matrices[[1]]
cov2 <- covariance_matrices[[2]]
cov3 <- covariance_matrices[[3]]
cov4 <- covariance_matrices[[4]]

Cor1 <- correlation_matrices[[1]]
Cor2 <- correlation_matrices[[2]]
Cor3 <- correlation_matrices[[3]]
Cor4 <- correlation_matrices[[4]]


# ------------------------------------------------------------------------------
# 5. Save state-specific matrices
# ------------------------------------------------------------------------------

for (state in seq_len(K)) {
  utils::write.csv(
    covariance_matrices[[state]],
    file.path(
      OUTPUT_DIR,
      sprintf("covariance_matrix_state_%d.csv", state)
    ),
    row.names = TRUE
  )

  utils::write.csv(
    correlation_matrices[[state]],
    file.path(
      OUTPUT_DIR,
      sprintf("correlation_matrix_state_%d.csv", state)
    ),
    row.names = TRUE
  )
}


# ------------------------------------------------------------------------------
# 6. Plot the four state-specific correlation matrices
# ------------------------------------------------------------------------------

heatmap_colours <- grDevices::colorRampPalette(
  c(
    "#2166AC",
    "#67A9CF",
    "white",
    "#EF8A62",
    "#B2182B"
  )
)(201)


plot_state_matrix <- function(
    R,
    state_number,
    use_raster = TRUE
) {
  number_of_ROIs <- ncol(R)

  # Reverse the row order so that the first ROI appears at the top.
  plot_matrix <- R[
    number_of_ROIs:1,
    ,
    drop = FALSE
  ]

  graphics::image(
    x = seq_len(number_of_ROIs),
    y = seq_len(number_of_ROIs),
    z = t(plot_matrix),
    zlim = c(-1, 1),
    col = heatmap_colours,
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = paste("State", state_number),
    useRaster = use_raster
  )

  graphics::axis(
    side = 1,
    at = seq_len(number_of_ROIs),
    labels = roi_names,
    las = 2,
    cex.axis = 0.55
  )

  graphics::axis(
    side = 2,
    at = seq_len(number_of_ROIs),
    labels = rev(roi_names),
    las = 2,
    cex.axis = 0.55
  )

  graphics::box()
}


draw_all_state_matrices <- function(use_raster = TRUE) {
  graphics::layout(
    matrix(
      c(
        1, 2, 5,
        3, 4, 5
      ),
      nrow = 2,
      byrow = TRUE
    ),
    widths = c(1, 1, 0.12)
  )

  for (state in seq_len(K)) {
    graphics::par(
      mar = c(9, 9, 3, 1)
    )

    plot_state_matrix(
      correlation_matrices[[state]],
      state_number = state,
      use_raster = use_raster
    )
  }

  # Shared colour legend.
  graphics::par(
    mar = c(5, 1, 3, 4)
  )

  legend_values <- seq(
    -1,
    1,
    length.out = length(heatmap_colours)
  )

  graphics::image(
    x = 1,
    y = legend_values,
    z = matrix(
      legend_values,
      nrow = 1,
      ncol = length(legend_values)
    ),
    col = heatmap_colours,
    zlim = c(-1, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    useRaster = use_raster
  )

  graphics::axis(
    side = 4,
    at = seq(-1, 1, by = 0.5),
    las = 1
  )

  graphics::mtext(
    "Pearson correlation",
    side = 4,
    line = 2.5
  )
}


# High-resolution PNG.
grDevices::png(
  filename = file.path(
    OUTPUT_DIR,
    "All_States_Correlation_Matrices.png"
  ),
  width = 4800,
  height = 3600,
  res = 300
)

draw_all_state_matrices(
  use_raster = TRUE
)

grDevices::dev.off()


# Vector PDF suitable for manuscript preparation.
grDevices::pdf(
  file = file.path(
    OUTPUT_DIR,
    "All_States_Correlation_Matrices.pdf"
  ),
  width = 16,
  height = 12,
  onefile = TRUE
)

draw_all_state_matrices(
  use_raster = FALSE
)

grDevices::dev.off()


# ------------------------------------------------------------------------------
# 7. State-specific connectivity and mean-signal summaries
# ------------------------------------------------------------------------------

state_summary_rows <- lapply(
  seq_len(K),
  function(state) {
    R <- correlation_matrices[[state]]
    correlations <- get_unique_correlations(R)

    mu <- as.numeric(mu_list[[state]])

    if (length(mu) != p) {
      stop(
        "The mean-vector length is incorrect for state ",
        state,
        "."
      )
    }

    if (any(!is.finite(mu))) {
      stop(
        "The mean vector contains non-finite values in state ",
        state,
        "."
      )
    }

    data.frame(
      state = state,

      mean_all_off_diagonal_correlations =
        safe_mean(correlations),

      mean_positive_correlation =
        mean_positive(correlations),

      mean_negative_correlation_signed =
        mean_negative_signed(correlations),

      mean_negative_correlation_magnitude =
        mean_negative_magnitude(correlations),

      number_positive_correlations =
        sum(correlations > 0),

      number_negative_correlations =
        sum(correlations < 0),

      number_zero_correlations =
        sum(correlations == 0),

      mean_all_ROI_signals =
        safe_mean(mu),

      mean_positive_ROI_signal =
        mean_positive(mu),

      mean_negative_ROI_signal_signed =
        mean_negative_signed(mu),

      mean_negative_ROI_signal_magnitude =
        mean_negative_magnitude(mu),

      number_positive_ROIs =
        sum(mu > 0),

      number_negative_ROIs =
        sum(mu < 0),

      number_zero_ROIs =
        sum(mu == 0),

      stringsAsFactors = FALSE
    )
  }
)

state_summary <- do.call(
  rbind,
  state_summary_rows
)

utils::write.csv(
  state_summary,
  file.path(
    OUTPUT_DIR,
    "state_correlation_and_mean_signal_summary.csv"
  ),
  row.names = FALSE
)


# Reviewer-friendly rounded table.
state_summary_rounded <- state_summary

continuous_columns <- setdiff(
  names(state_summary_rounded),
  c(
    "state",
    "number_positive_correlations",
    "number_negative_correlations",
    "number_zero_correlations",
    "number_positive_ROIs",
    "number_negative_ROIs",
    "number_zero_ROIs"
  )
)

state_summary_rounded[
  continuous_columns
] <- lapply(
  state_summary_rounded[
    continuous_columns
  ],
  function(x) round(x, 6)
)

utils::write.csv(
  state_summary_rounded,
  file.path(
    OUTPUT_DIR,
    "state_correlation_and_mean_signal_summary_rounded.csv"
  ),
  row.names = FALSE
)

print(
  state_summary_rounded,
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8. Save state-specific mean signal vectors
# ------------------------------------------------------------------------------

mean_signal_table <- data.frame(
  ROI = roi_names,
  State_1 = as.numeric(mu_list[[1]]),
  State_2 = as.numeric(mu_list[[2]]),
  State_3 = as.numeric(mu_list[[3]]),
  State_4 = as.numeric(mu_list[[4]]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

utils::write.csv(
  mean_signal_table,
  file.path(
    OUTPUT_DIR,
    "state_specific_mean_signal_vectors.csv"
  ),
  row.names = FALSE
)


# Long-format table can be convenient for plotting or supplementary material.
mean_signal_long <- do.call(
  rbind,
  lapply(
    seq_len(K),
    function(state) {
      data.frame(
        state = state,
        ROI = roi_names,
        mean_signal = as.numeric(mu_list[[state]]),
        stringsAsFactors = FALSE
      )
    }
  )
)

utils::write.csv(
  mean_signal_long,
  file.path(
    OUTPUT_DIR,
    "state_specific_mean_signal_vectors_long.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 9. Reproducibility information
# ------------------------------------------------------------------------------

analysis_settings <- data.frame(
  setting = c(
    "model_file",
    "number_of_states",
    "number_of_ROIs",
    "number_of_unique_ROI_pairs",
    "correlation_type",
    "post_estimation_smoothing"
  ),
  value = c(
    normalizePath(
      MODEL_FILE,
      winslash = "/",
      mustWork = TRUE
    ),
    K,
    p,
    p * (p - 1) / 2,
    "Pearson correlations obtained using cov2cor on each fitted covariance matrix",
    "None"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  analysis_settings,
  file.path(
    OUTPUT_DIR,
    "state_characterization_analysis_settings.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "state_characterization_sessionInfo.txt"
  )
)


# ------------------------------------------------------------------------------
# 10. Completion message
# ------------------------------------------------------------------------------

cat("\nState-characterization analysis completed successfully.\n")
cat(
  "Results directory:\n",
  normalizePath(
    OUTPUT_DIR,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
