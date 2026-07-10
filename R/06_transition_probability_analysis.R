# ==============================================================================
# FINAL HSMM RUN-LEVEL STATE-TRANSITION ANALYSIS
# Reviewer-shareable reproducible script
#
# Input:
#   decoded_states_final_K04.csv
#
# Methodological notes:
#   - Repeated Viterbi labels are collapsed into explicit state runs.
#   - Only transitions between consecutive different runs are counted.
#   - Self-transitions are structurally excluded, as appropriate for an HSMM.
#   - First and last runs are retained because boundary censoring affects run
#     duration, not the identity of observed transitions between adjacent runs.
#   - Group labels are permuted at the participant level.
#   - Bootstrap resampling is stratified within TD and Dyslexia groups.
#   - No package beyond base R is required.
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

DECODED_FILE <- file.path(
  PROJECT_DIR,
  "decoded_states_final_K04.csv"
)

OUTPUT_DIR <- file.path(
  PROJECT_DIR,
  "Updated_Transition_Results"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

K <- 4L
N_PERMUTATIONS <- 5000L
PERMUTATION_SEED <- 20260713L
N_BOOTSTRAPS <- 2000L
BOOTSTRAP_SEED <- 20260714L
DESCRIPTIVE_DIFFERENCE_THRESHOLD <- 0.075
CREATE_FIGURES <- TRUE


# ------------------------------------------------------------------------------
# 2. Read and validate decoded-state data
# ------------------------------------------------------------------------------

if (!file.exists(DECODED_FILE)) {
  stop("The decoded-state file was not found:\n", DECODED_FILE)
}

decoded_states <- utils::read.csv(
  DECODED_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "subject_id",
  "group",
  "file",
  "timepoint",
  "state"
)

missing_columns <- setdiff(
  required_columns,
  names(decoded_states)
)

if (length(missing_columns) > 0L) {
  stop(
    "Missing column(s) in decoded_states_final_K04.csv: ",
    paste(missing_columns, collapse = ", ")
  )
}

decoded_states$subject_id <- as.integer(decoded_states$subject_id)
decoded_states$timepoint <- as.integer(decoded_states$timepoint)
decoded_states$state <- as.integer(decoded_states$state)

if (
  any(!is.finite(decoded_states$subject_id)) ||
  any(!is.finite(decoded_states$timepoint)) ||
  any(!is.finite(decoded_states$state))
) {
  stop("The decoded-state file contains invalid numeric values.")
}

if (any(!decoded_states$state %in% seq_len(K))) {
  stop("State labels must be integers from 1 to ", K, ".")
}


# ------------------------------------------------------------------------------
# 3. Standardize group labels
# ------------------------------------------------------------------------------

standardize_group <- function(x) {
  original <- trimws(as.character(x))
  lower <- tolower(original)
  output <- original

  output[
    lower %in% c(
      "control", "controls", "ctrl", "td",
      "typically developing", "typically_developing"
    )
  ] <- "TD"

  output[
    lower %in% c(
      "patient", "patients", "dyslexia", "dyslexic"
    )
  ] <- "Dyslexia"

  output
}

decoded_states$group <- standardize_group(decoded_states$group)
GROUP_LEVELS <- c("TD", "Dyslexia")

unexpected_groups <- setdiff(
  unique(decoded_states$group),
  GROUP_LEVELS
)

if (length(unexpected_groups) > 0L) {
  stop(
    "Unexpected group label(s): ",
    paste(unexpected_groups, collapse = ", "),
    ". Expected groups are TD and Dyslexia."
  )
}


# ------------------------------------------------------------------------------
# 4. Reconstruct participant sequences
# ------------------------------------------------------------------------------

decoded_states <- decoded_states[
  order(decoded_states$subject_id, decoded_states$timepoint),
  ,
  drop = FALSE
]
rownames(decoded_states) <- NULL

subject_ids <- sort(unique(decoded_states$subject_id))
n_subjects <- length(subject_ids)

subject_sequences <- vector("list", n_subjects)
names(subject_sequences) <- as.character(subject_ids)
subject_information_rows <- vector("list", n_subjects)

for (i in seq_along(subject_ids)) {
  subject_id <- subject_ids[i]
  one_subject <- decoded_states[
    decoded_states$subject_id == subject_id,
    ,
    drop = FALSE
  ]

  if (!identical(one_subject$timepoint, seq_len(nrow(one_subject)))) {
    stop(
      "Time points are not consecutive from 1 for subject ",
      subject_id,
      "."
    )
  }

  if (length(unique(one_subject$group)) != 1L) {
    stop("More than one group label was found for subject ", subject_id, ".")
  }

  if (length(unique(one_subject$file)) != 1L) {
    stop("More than one filename was found for subject ", subject_id, ".")
  }

  subject_sequences[[i]] <- as.integer(one_subject$state)

  subject_information_rows[[i]] <- data.frame(
    subject_index = i,
    subject_id = subject_id,
    group = one_subject$group[1],
    file = one_subject$file[1],
    number_of_timepoints = nrow(one_subject),
    stringsAsFactors = FALSE
  )
}

subject_information <- do.call(rbind, subject_information_rows)
rownames(subject_information) <- NULL
subject_information$group <- factor(
  subject_information$group,
  levels = GROUP_LEVELS
)

if (any(is.na(subject_information$group))) {
  stop("At least one participant has an invalid group label.")
}

TD_indices <- which(subject_information$group == "TD")
Dyslexia_indices <- which(subject_information$group == "Dyslexia")

if (length(TD_indices) == 0L || length(Dyslexia_indices) == 0L) {
  stop("Both TD and Dyslexia groups must contain participants.")
}


# ------------------------------------------------------------------------------
# 5. Participant-level run-transition counts
# ------------------------------------------------------------------------------

state_names <- paste0("State_", seq_len(K))

empty_transition_matrix <- function() {
  matrix(
    0,
    nrow = K,
    ncol = K,
    dimnames = list(From = state_names, To = state_names)
  )
}

count_run_level_transitions <- function(state_sequence) {
  state_sequence <- as.integer(state_sequence)
  state_sequence <- state_sequence[!is.na(state_sequence)]
  count_matrix <- empty_transition_matrix()

  if (length(state_sequence) < 2L) {
    return(count_matrix)
  }

  run_states <- as.integer(rle(state_sequence)$values)

  if (length(run_states) < 2L) {
    return(count_matrix)
  }

  from_states <- run_states[-length(run_states)]
  to_states <- run_states[-1L]

  for (transition_index in seq_along(from_states)) {
    from_state <- from_states[transition_index]
    to_state <- to_states[transition_index]

    if (
      from_state >= 1L && from_state <= K &&
      to_state >= 1L && to_state <= K &&
      from_state != to_state
    ) {
      count_matrix[from_state, to_state] <-
        count_matrix[from_state, to_state] + 1L
    }
  }

  diag(count_matrix) <- 0
  count_matrix
}

subject_transition_counts <- lapply(
  subject_sequences,
  count_run_level_transitions
)

names(subject_transition_counts) <- as.character(subject_ids)


# ------------------------------------------------------------------------------
# 6. Row-normalization: P(To | From)
# ------------------------------------------------------------------------------

row_normalize_transition_counts <- function(count_matrix) {
  count_matrix <- as.matrix(count_matrix)

  probability_matrix <- matrix(
    NA_real_,
    nrow = K,
    ncol = K,
    dimnames = dimnames(count_matrix)
  )

  row_totals <- rowSums(count_matrix, na.rm = TRUE)
  estimable_rows <- which(row_totals > 0)

  if (length(estimable_rows) > 0L) {
    probability_matrix[estimable_rows, ] <-
      count_matrix[estimable_rows, , drop = FALSE] /
      row_totals[estimable_rows]

    probability_matrix[cbind(estimable_rows, estimable_rows)] <- 0
  }

  probability_matrix
}

subject_transition_probabilities <- lapply(
  subject_transition_counts,
  row_normalize_transition_counts
)


# ------------------------------------------------------------------------------
# 7. Save participant-level transition results
# ------------------------------------------------------------------------------

subject_transition_rows <- list()
row_counter <- 1L

for (i in seq_len(n_subjects)) {
  count_matrix <- subject_transition_counts[[i]]
  probability_matrix <- subject_transition_probabilities[[i]]

  for (from_state in seq_len(K)) {
    for (to_state in seq_len(K)) {
      if (from_state == to_state) {
        next
      }

      subject_transition_rows[[row_counter]] <- data.frame(
        subject_id = subject_information$subject_id[i],
        group = as.character(subject_information$group[i]),
        file = subject_information$file[i],
        from_state = from_state,
        to_state = to_state,
        transition_count = count_matrix[from_state, to_state],
        total_transitions_from_state = sum(count_matrix[from_state, ]),
        conditional_probability = probability_matrix[from_state, to_state],
        stringsAsFactors = FALSE
      )

      row_counter <- row_counter + 1L
    }
  }
}

subject_transition_table <- do.call(rbind, subject_transition_rows)

utils::write.csv(
  subject_transition_table,
  file.path(
    OUTPUT_DIR,
    "subject_level_transition_counts_and_probabilities.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 8. Group-level pooled transition probabilities
# ------------------------------------------------------------------------------

sum_transition_matrices <- function(participant_indices) {
  total_matrix <- empty_transition_matrix()

  for (participant_index in participant_indices) {
    total_matrix <- total_matrix +
      subject_transition_counts[[participant_index]]
  }

  total_matrix
}

TD_counts <- sum_transition_matrices(TD_indices)
Dyslexia_counts <- sum_transition_matrices(Dyslexia_indices)

TD_probabilities <- row_normalize_transition_counts(TD_counts)
Dyslexia_probabilities <- row_normalize_transition_counts(Dyslexia_counts)

observed_difference <- Dyslexia_probabilities - TD_probabilities
diag(observed_difference) <- NA_real_

utils::write.csv(
  TD_counts,
  file.path(OUTPUT_DIR, "TD_transition_count_matrix.csv"),
  row.names = TRUE
)

utils::write.csv(
  Dyslexia_counts,
  file.path(OUTPUT_DIR, "Dyslexia_transition_count_matrix.csv"),
  row.names = TRUE
)

utils::write.csv(
  TD_probabilities,
  file.path(OUTPUT_DIR, "TD_transition_probability_matrix.csv"),
  row.names = TRUE
)

utils::write.csv(
  Dyslexia_probabilities,
  file.path(OUTPUT_DIR, "Dyslexia_transition_probability_matrix.csv"),
  row.names = TRUE
)

utils::write.csv(
  observed_difference,
  file.path(
    OUTPUT_DIR,
    "transition_probability_difference_Dyslexia_minus_TD.csv"
  ),
  row.names = TRUE
)


# ------------------------------------------------------------------------------
# 9. Participant-level permutation test
# ------------------------------------------------------------------------------

aggregate_by_labels <- function(labels) {
  labels <- as.character(labels)
  TD_total <- empty_transition_matrix()
  Dyslexia_total <- empty_transition_matrix()

  for (i in seq_len(n_subjects)) {
    if (labels[i] == "TD") {
      TD_total <- TD_total + subject_transition_counts[[i]]
    } else if (labels[i] == "Dyslexia") {
      Dyslexia_total <- Dyslexia_total + subject_transition_counts[[i]]
    } else {
      stop("Unexpected label encountered during aggregation.")
    }
  }

  list(TD = TD_total, Dyslexia = Dyslexia_total)
}

set.seed(PERMUTATION_SEED)

permuted_differences <- array(
  NA_real_,
  dim = c(K, K, N_PERMUTATIONS),
  dimnames = list(
    From = state_names,
    To = state_names,
    Permutation = NULL
  )
)

original_group_labels <- as.character(subject_information$group)

for (permutation_index in seq_len(N_PERMUTATIONS)) {
  permuted_labels <- sample(
    original_group_labels,
    size = n_subjects,
    replace = FALSE
  )

  permuted_counts <- aggregate_by_labels(permuted_labels)

  permuted_TD_probabilities <- row_normalize_transition_counts(
    permuted_counts$TD
  )

  permuted_Dyslexia_probabilities <- row_normalize_transition_counts(
    permuted_counts$Dyslexia
  )

  one_difference <-
    permuted_Dyslexia_probabilities - permuted_TD_probabilities

  diag(one_difference) <- NA_real_
  permuted_differences[, , permutation_index] <- one_difference

  if (
    permutation_index %% 500L == 0L ||
    permutation_index == N_PERMUTATIONS
  ) {
    cat(
      "Completed permutation",
      permutation_index,
      "of",
      N_PERMUTATIONS,
      "\n"
    )
  }
}

permutation_p_values <- matrix(
  NA_real_,
  nrow = K,
  ncol = K,
  dimnames = list(From = state_names, To = state_names)
)

valid_permutation_counts <- matrix(
  0L,
  nrow = K,
  ncol = K,
  dimnames = dimnames(permutation_p_values)
)

for (from_state in seq_len(K)) {
  for (to_state in seq_len(K)) {
    if (
      from_state == to_state ||
      !is.finite(observed_difference[from_state, to_state])
    ) {
      next
    }

    null_values <- permuted_differences[
      from_state,
      to_state,
      
    ]

    null_values <- null_values[is.finite(null_values)]
    valid_permutation_counts[from_state, to_state] <- length(null_values)

    if (length(null_values) == 0L) {
      next
    }

    permutation_p_values[from_state, to_state] <-
      (
        1 + sum(
          abs(null_values) >=
            abs(observed_difference[from_state, to_state])
        )
      ) /
      (1 + length(null_values))
  }
}

FDR_adjusted_p_values <- matrix(
  NA_real_,
  nrow = K,
  ncol = K,
  dimnames = dimnames(permutation_p_values)
)

estimable_p_indices <- which(is.finite(permutation_p_values))

if (length(estimable_p_indices) > 0L) {
  FDR_adjusted_p_values[estimable_p_indices] <- stats::p.adjust(
    permutation_p_values[estimable_p_indices],
    method = "BH"
  )
}


# ------------------------------------------------------------------------------
# 10. Stratified participant bootstrap
# ------------------------------------------------------------------------------

set.seed(BOOTSTRAP_SEED)

bootstrap_TD_probabilities <- array(
  NA_real_,
  dim = c(K, K, N_BOOTSTRAPS)
)

bootstrap_Dyslexia_probabilities <- array(
  NA_real_,
  dim = c(K, K, N_BOOTSTRAPS)
)

bootstrap_differences <- array(
  NA_real_,
  dim = c(K, K, N_BOOTSTRAPS)
)

for (bootstrap_index in seq_len(N_BOOTSTRAPS)) {
  sampled_TD_indices <- sample(
    TD_indices,
    size = length(TD_indices),
    replace = TRUE
  )

  sampled_Dyslexia_indices <- sample(
    Dyslexia_indices,
    size = length(Dyslexia_indices),
    replace = TRUE
  )

  bootstrap_TD_counts <- sum_transition_matrices(sampled_TD_indices)
  bootstrap_Dyslexia_counts <- sum_transition_matrices(
    sampled_Dyslexia_indices
  )

  bootstrap_TD_matrix <- row_normalize_transition_counts(
    bootstrap_TD_counts
  )

  bootstrap_Dyslexia_matrix <- row_normalize_transition_counts(
    bootstrap_Dyslexia_counts
  )

  bootstrap_difference_matrix <-
    bootstrap_Dyslexia_matrix - bootstrap_TD_matrix

  diag(bootstrap_difference_matrix) <- NA_real_

  bootstrap_TD_probabilities[, , bootstrap_index] <-
    bootstrap_TD_matrix

  bootstrap_Dyslexia_probabilities[, , bootstrap_index] <-
    bootstrap_Dyslexia_matrix

  bootstrap_differences[, , bootstrap_index] <-
    bootstrap_difference_matrix

  if (
    bootstrap_index %% 500L == 0L ||
    bootstrap_index == N_BOOTSTRAPS
  ) {
    cat(
      "Completed bootstrap",
      bootstrap_index,
      "of",
      N_BOOTSTRAPS,
      "\n"
    )
  }
}

safe_quantile <- function(values, probability) {
  values <- values[is.finite(values)]

  if (length(values) == 0L) {
    return(NA_real_)
  }

  unname(
    stats::quantile(
      values,
      probs = probability,
      names = FALSE,
      type = 7
    )
  )
}


# ------------------------------------------------------------------------------
# 11. Final transition-comparison table
# ------------------------------------------------------------------------------

transition_result_rows <- list()
result_counter <- 1L

for (from_state in seq_len(K)) {
  for (to_state in seq_len(K)) {
    if (from_state == to_state) {
      next
    }

    TD_bootstrap_values <- bootstrap_TD_probabilities[
      from_state,
      to_state,
      
    ]

    Dyslexia_bootstrap_values <- bootstrap_Dyslexia_probabilities[
      from_state,
      to_state,
      
    ]

    difference_bootstrap_values <- bootstrap_differences[
      from_state,
      to_state,
      
    ]

    transition_result_rows[[result_counter]] <- data.frame(
      from_state = from_state,
      to_state = to_state,
      TD_count = TD_counts[from_state, to_state],
      TD_probability = TD_probabilities[from_state, to_state],
      TD_CI_lower = safe_quantile(TD_bootstrap_values, 0.025),
      TD_CI_upper = safe_quantile(TD_bootstrap_values, 0.975),
      TD_valid_bootstraps = sum(is.finite(TD_bootstrap_values)),
      Dyslexia_count = Dyslexia_counts[from_state, to_state],
      Dyslexia_probability =
        Dyslexia_probabilities[from_state, to_state],
      Dyslexia_CI_lower =
        safe_quantile(Dyslexia_bootstrap_values, 0.025),
      Dyslexia_CI_upper =
        safe_quantile(Dyslexia_bootstrap_values, 0.975),
      Dyslexia_valid_bootstraps =
        sum(is.finite(Dyslexia_bootstrap_values)),
      difference_Dyslexia_minus_TD =
        observed_difference[from_state, to_state],
      difference_CI_lower =
        safe_quantile(difference_bootstrap_values, 0.025),
      difference_CI_upper =
        safe_quantile(difference_bootstrap_values, 0.975),
      difference_valid_bootstraps =
        sum(is.finite(difference_bootstrap_values)),
      permutation_p_value =
        permutation_p_values[from_state, to_state],
      FDR_adjusted_p_value =
        FDR_adjusted_p_values[from_state, to_state],
      valid_permutations =
        valid_permutation_counts[from_state, to_state],
      stringsAsFactors = FALSE
    )

    result_counter <- result_counter + 1L
  }
}

TransitionResults <- do.call(rbind, transition_result_rows)

TransitionResults$significant_FDR_0_05 <- ifelse(
  is.finite(TransitionResults$FDR_adjusted_p_value),
  TransitionResults$FDR_adjusted_p_value < 0.05,
  NA
)

TransitionResults$large_descriptive_difference <- ifelse(
  is.finite(TransitionResults$difference_Dyslexia_minus_TD),
  abs(TransitionResults$difference_Dyslexia_minus_TD) >=
    DESCRIPTIVE_DIFFERENCE_THRESHOLD,
  NA
)

TransitionResults$direction <- ifelse(
  is.na(TransitionResults$difference_Dyslexia_minus_TD),
  "Not estimable",
  ifelse(
    TransitionResults$difference_Dyslexia_minus_TD > 0,
    "Higher in Dyslexia",
    ifelse(
      TransitionResults$difference_Dyslexia_minus_TD < 0,
      "Higher in TD",
      "No difference"
    )
  )
)

TransitionResults <- TransitionResults[
  order(
    -abs(TransitionResults$difference_Dyslexia_minus_TD),
    na.last = TRUE
  ),
  ,
  drop = FALSE
]
rownames(TransitionResults) <- NULL

utils::write.csv(
  TransitionResults,
  file.path(
    OUTPUT_DIR,
    "transition_probability_group_comparison.csv"
  ),
  row.names = FALSE
)

TopTransitionDifferences <- TransitionResults[
  !is.na(TransitionResults$large_descriptive_difference) &
    TransitionResults$large_descriptive_difference,
  ,
  drop = FALSE
]

utils::write.csv(
  TopTransitionDifferences,
  file.path(
    OUTPUT_DIR,
    "large_transition_probability_differences.csv"
  ),
  row.names = FALSE
)

TransitionResults_rounded <- TransitionResults

columns_to_round <- c(
  "TD_probability",
  "TD_CI_lower",
  "TD_CI_upper",
  "Dyslexia_probability",
  "Dyslexia_CI_lower",
  "Dyslexia_CI_upper",
  "difference_Dyslexia_minus_TD",
  "difference_CI_lower",
  "difference_CI_upper",
  "permutation_p_value",
  "FDR_adjusted_p_value"
)

TransitionResults_rounded[columns_to_round] <- lapply(
  TransitionResults_rounded[columns_to_round],
  function(x) round(x, 6)
)


# ------------------------------------------------------------------------------
# 12. Base-R heatmaps
# ------------------------------------------------------------------------------

plot_probability_heatmap <- function(
    probability_matrix,
    main_title,
    maximum_probability
) {
  plot_matrix <- probability_matrix[K:1, , drop = FALSE]
  colours <- grDevices::colorRampPalette(
    c("white", "orange", "red")
  )(100)

  graphics::image(
    x = seq_len(K),
    y = seq_len(K),
    z = t(plot_matrix),
    zlim = c(0, maximum_probability),
    col = colours,
    axes = FALSE,
    xlab = "To state",
    ylab = "From state",
    main = main_title
  )

  graphics::axis(
    side = 1,
    at = seq_len(K),
    labels = paste0("State ", seq_len(K))
  )

  graphics::axis(
    side = 2,
    at = seq_len(K),
    labels = rev(paste0("State ", seq_len(K)))
  )

  graphics::box()

  for (from_state in seq_len(K)) {
    for (to_state in seq_len(K)) {
      value <- probability_matrix[from_state, to_state]

      if (is.finite(value)) {
        graphics::text(
          x = to_state,
          y = K - from_state + 1L,
          labels = sprintf("%.3f", value),
          cex = 0.9
        )
      }
    }
  }
}

plot_difference_heatmap <- function(difference_matrix, main_title) {
  finite_values <- difference_matrix[is.finite(difference_matrix)]

  if (length(finite_values) == 0L) {
    stop("No finite transition differences are available for plotting.")
  }

  maximum_absolute_difference <- max(abs(finite_values))
  plot_matrix <- difference_matrix[K:1, , drop = FALSE]
  colours <- grDevices::colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(101)

  graphics::image(
    x = seq_len(K),
    y = seq_len(K),
    z = t(plot_matrix),
    zlim = c(
      -maximum_absolute_difference,
      maximum_absolute_difference
    ),
    col = colours,
    axes = FALSE,
    xlab = "To state",
    ylab = "From state",
    main = main_title
  )

  graphics::axis(
    side = 1,
    at = seq_len(K),
    labels = paste0("State ", seq_len(K))
  )

  graphics::axis(
    side = 2,
    at = seq_len(K),
    labels = rev(paste0("State ", seq_len(K)))
  )

  graphics::box()

  for (from_state in seq_len(K)) {
    for (to_state in seq_len(K)) {
      value <- difference_matrix[from_state, to_state]

      if (is.finite(value)) {
        graphics::text(
          x = to_state,
          y = K - from_state + 1L,
          labels = sprintf("%+.3f", value),
          cex = 0.9
        )
      }
    }
  }
}

if (isTRUE(CREATE_FIGURES)) {
  maximum_probability <- max(
    c(TD_probabilities, Dyslexia_probabilities),
    na.rm = TRUE
  )

  grDevices::png(
    filename = file.path(
      OUTPUT_DIR,
      "Group_Transition_Probability_Heatmaps.png"
    ),
    width = 3200,
    height = 1600,
    res = 300
  )

  graphics::par(
    mfrow = c(1, 2),
    mar = c(5, 5, 4, 2)
  )

  plot_probability_heatmap(
    TD_probabilities,
    "TD Transition Probabilities",
    maximum_probability
  )

  plot_probability_heatmap(
    Dyslexia_probabilities,
    "Dyslexia Transition Probabilities",
    maximum_probability
  )

  grDevices::dev.off()

  grDevices::pdf(
    file = file.path(
      OUTPUT_DIR,
      "Group_Transition_Probability_Heatmaps.pdf"
    ),
    width = 12,
    height = 6,
    onefile = TRUE
  )

  graphics::par(
    mfrow = c(1, 2),
    mar = c(5, 5, 4, 2)
  )

  plot_probability_heatmap(
    TD_probabilities,
    "TD Transition Probabilities",
    maximum_probability
  )

  plot_probability_heatmap(
    Dyslexia_probabilities,
    "Dyslexia Transition Probabilities",
    maximum_probability
  )

  grDevices::dev.off()

  grDevices::png(
    filename = file.path(
      OUTPUT_DIR,
      "Transition_Difference_Dyslexia_minus_TD.png"
    ),
    width = 1800,
    height = 1800,
    res = 300
  )

  graphics::par(mar = c(5, 5, 4, 2))

  plot_difference_heatmap(
    observed_difference,
    "Transition Difference: Dyslexia minus TD"
  )

  grDevices::dev.off()

  grDevices::pdf(
    file = file.path(
      OUTPUT_DIR,
      "Transition_Difference_Dyslexia_minus_TD.pdf"
    ),
    width = 7,
    height = 7,
    onefile = TRUE
  )

  graphics::par(mar = c(5, 5, 4, 2))

  plot_difference_heatmap(
    observed_difference,
    "Transition Difference: Dyslexia minus TD"
  )

  grDevices::dev.off()
}


# ------------------------------------------------------------------------------
# 13. Reproducibility information
# ------------------------------------------------------------------------------

analysis_settings <- data.frame(
  setting = c(
    "decoded_state_file",
    "selected_K",
    "transition_definition",
    "self_transitions",
    "boundary_runs_removed",
    "group_probability_estimator",
    "number_of_permutations",
    "permutation_seed",
    "permutation_unit",
    "permutation_test",
    "multiple_testing_adjustment",
    "number_of_bootstraps",
    "bootstrap_seed",
    "bootstrap_method",
    "descriptive_difference_threshold"
  ),
  value = c(
    normalizePath(
      DECODED_FILE,
      winslash = "/",
      mustWork = TRUE
    ),
    K,
    paste0(
      "Transitions between consecutive state runs after collapsing ",
      "repeated Viterbi labels with rle"
    ),
    "Excluded structurally",
    paste0(
      "No; boundary censoring concerns duration, whereas observed ",
      "between-run transitions are fully observed"
    ),
    paste0(
      "Pooled run-level transition counts within group, followed by ",
      "row normalization: P(To | From)"
    ),
    N_PERMUTATIONS,
    PERMUTATION_SEED,
    "Participant-level group labels",
    "Two-sided cell-wise test with plus-one correction",
    "Benjamini-Hochberg FDR across estimable off-diagonal cells",
    N_BOOTSTRAPS,
    BOOTSTRAP_SEED,
    paste0(
      "Stratified participant bootstrap within TD and Dyslexia; ",
      "percentile 95% confidence intervals"
    ),
    DESCRIPTIVE_DIFFERENCE_THRESHOLD
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  analysis_settings,
  file.path(OUTPUT_DIR, "transition_analysis_settings.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "transition_analysis_sessionInfo.txt"
  )
)


# ------------------------------------------------------------------------------
# 14. Final console output
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL RUN-LEVEL TRANSITION ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Participants:", n_subjects, "\n")
cat("TD participants:", length(TD_indices), "\n")
cat("Dyslexia participants:", length(Dyslexia_indices), "\n")
cat("Permutation replicates:", N_PERMUTATIONS, "\n")
cat("Bootstrap replicates:", N_BOOTSTRAPS, "\n")
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

cat("TD transition probabilities P(To | From):\n")
print(round(TD_probabilities, 4))

cat("\nDyslexia transition probabilities P(To | From):\n")
print(round(Dyslexia_probabilities, 4))

cat("\nTransition-comparison results, ordered by absolute difference:\n")
print(TransitionResults_rounded, row.names = FALSE)
