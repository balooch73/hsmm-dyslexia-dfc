# ==============================================================================
# HSMM STATE OCCUPANCY AND VITERBI-SEQUENCE ANALYSIS
# Corrected reviewer-shareable reproducible script
#
# This version uses decoded_states_final_K04.csv as the single source of truth.
# Subject-specific sequence lengths are derived directly from the decoded data,
# so no "number_of_timepoints" column is required in subject_information.csv.
#
# Required package:
#   Base R only
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
  "Updated_Occupancy_Results_v2"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

K <- 4L
N_PERMUTATIONS <- 5000L
PERMUTATION_SEED <- 20260711L
MAX_TIME_TO_PLOT <- 200L

STATE_COLOURS <- c(
  "#0072B2",
  "#E69F00",
  "#009E73",
  "#CC79A7"
)

names(STATE_COLOURS) <- paste0(
  "State ",
  seq_len(K)
)


# ------------------------------------------------------------------------------
# 2. Read and validate decoded-state data
# ------------------------------------------------------------------------------

if (!file.exists(DECODED_FILE)) {
  stop(
    "The decoded-state file was not found:\n",
    DECODED_FILE
  )
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

decoded_states$subject_id <- as.integer(
  decoded_states$subject_id
)

decoded_states$timepoint <- as.integer(
  decoded_states$timepoint
)

decoded_states$state <- as.integer(
  decoded_states$state
)

if (
  any(!is.finite(decoded_states$subject_id)) ||
  any(!is.finite(decoded_states$timepoint)) ||
  any(!is.finite(decoded_states$state))
) {
  stop(
    "The decoded-state file contains invalid numeric values."
  )
}

if (any(!decoded_states$state %in% seq_len(K))) {
  stop(
    "State labels must be integers from 1 to ",
    K,
    "."
  )
}


# ------------------------------------------------------------------------------
# 3. Standardize group labels
# ------------------------------------------------------------------------------

standardize_group <- function(x) {
  original <- trimws(
    as.character(x)
  )

  lower <- tolower(
    original
  )

  out <- original

  out[
    lower %in% c(
      "control",
      "controls",
      "ctrl",
      "td",
      "typically developing",
      "typically_developing"
    )
  ] <- "TD"

  out[
    lower %in% c(
      "patient",
      "patients",
      "dyslexia",
      "dyslexic"
    )
  ] <- "Dyslexia"

  out
}

decoded_states$group <- standardize_group(
  decoded_states$group
)

expected_groups <- c(
  "TD",
  "Dyslexia"
)

unexpected_groups <- setdiff(
  unique(decoded_states$group),
  expected_groups
)

if (length(unexpected_groups) > 0L) {
  stop(
    "Unexpected group label(s): ",
    paste(unexpected_groups, collapse = ", "),
    ". Expected labels are TD and Dyslexia."
  )
}


# ------------------------------------------------------------------------------
# 4. Sort data and derive subject information directly from decoded data
# ------------------------------------------------------------------------------

decoded_states <- decoded_states[
  order(
    decoded_states$subject_id,
    decoded_states$timepoint
  ),
  ,
  drop = FALSE
]

rownames(decoded_states) <- NULL

subject_ids <- sort(
  unique(decoded_states$subject_id)
)

subject_sequences <- vector(
  mode = "list",
  length = length(subject_ids)
)

names(subject_sequences) <- as.character(
  subject_ids
)

subject_information_rows <- vector(
  mode = "list",
  length = length(subject_ids)
)

for (i in seq_along(subject_ids)) {
  subject_id <- subject_ids[i]

  one_subject <- decoded_states[
    decoded_states$subject_id == subject_id,
    ,
    drop = FALSE
  ]

  expected_timepoints <- seq_len(
    nrow(one_subject)
  )

  if (!identical(one_subject$timepoint, expected_timepoints)) {
    stop(
      "Time points are not consecutive from 1 for subject ",
      subject_id,
      "."
    )
  }

  if (length(unique(one_subject$group)) != 1L) {
    stop(
      "More than one group label was found for subject ",
      subject_id,
      "."
    )
  }

  if (length(unique(one_subject$file)) != 1L) {
    stop(
      "More than one filename was found for subject ",
      subject_id,
      "."
    )
  }

  subject_sequences[[i]] <- one_subject

  subject_information_rows[[i]] <- data.frame(
    subject_id = subject_id,
    group = one_subject$group[1],
    file = one_subject$file[1],
    number_of_timepoints = nrow(one_subject),
    stringsAsFactors = FALSE
  )
}

subject_information <- do.call(
  rbind,
  subject_information_rows
)

rownames(subject_information) <- NULL

Lengths <- as.integer(
  subject_information$number_of_timepoints
)

names(Lengths) <- as.character(
  subject_information$subject_id
)

n_subjects <- nrow(
  subject_information
)

if (sum(Lengths) != nrow(decoded_states)) {
  stop(
    "The sum of subject-specific sequence lengths does not match ",
    "the total number of decoded observations."
  )
}

utils::write.csv(
  subject_information,
  file.path(
    OUTPUT_DIR,
    "subject_information_derived_from_decoded_states.csv"
  ),
  row.names = FALSE
)

cat("Decoded data validated successfully.\n")
cat("Number of subjects:", n_subjects, "\n")
cat("Total decoded observations:", nrow(decoded_states), "\n")
cat(
  "TD participants:",
  sum(subject_information$group == "TD"),
  "\n"
)
cat(
  "Dyslexia participants:",
  sum(subject_information$group == "Dyslexia"),
  "\n"
)


# ------------------------------------------------------------------------------
# 5. Save subject-level Viterbi sequences
# ------------------------------------------------------------------------------

StateSeqDF <- decoded_states[
  ,
  c(
    "subject_id",
    "group",
    "file",
    "timepoint",
    "state"
  ),
  drop = FALSE
]

utils::write.csv(
  StateSeqDF,
  file.path(
    OUTPUT_DIR,
    "subject_level_Viterbi_state_sequences.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 6. Calculate subject-level state occupancy
# ------------------------------------------------------------------------------

occupancy_matrix <- matrix(
  0,
  nrow = n_subjects,
  ncol = K,
  dimnames = list(
    paste0(
      "Subject_",
      subject_information$subject_id
    ),
    paste0(
      "State_",
      seq_len(K)
    )
  )
)

for (i in seq_len(n_subjects)) {
  states_i <- subject_sequences[[i]]$state

  occupancy_matrix[i, ] <- (
    tabulate(
      states_i,
      nbins = K
    ) /
      length(states_i)
  ) * 100
}

occupancy_sums <- rowSums(
  occupancy_matrix
)

if (any(abs(occupancy_sums - 100) > 1e-8)) {
  stop(
    "At least one participant's occupancy percentages do not sum to 100."
  )
}

subject_occupancy <- data.frame(
  subject_id = subject_information$subject_id,
  group = subject_information$group,
  file = subject_information$file,
  number_of_timepoints =
    subject_information$number_of_timepoints,
  occupancy_matrix,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

utils::write.csv(
  subject_occupancy,
  file.path(
    OUTPUT_DIR,
    "subject_level_state_occupancy_percent.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7. Group-level occupancy summaries
# ------------------------------------------------------------------------------

group_order <- c(
  "TD",
  "Dyslexia"
)

group_summary_rows <- list()
counter <- 1L

for (group_name in group_order) {
  group_indices <- which(
    subject_occupancy$group == group_name
  )

  if (length(group_indices) == 0L) {
    stop(
      "No participants were found in group: ",
      group_name
    )
  }

  for (state in seq_len(K)) {
    values <- occupancy_matrix[
      group_indices,
      state
    ]

    group_summary_rows[[counter]] <- data.frame(
      group = group_name,
      state = state,
      n = length(values),
      mean_percent = mean(values),
      standard_deviation = stats::sd(values),
      standard_error = stats::sd(values) /
        sqrt(length(values)),
      median_percent = stats::median(values),
      IQR_low = unname(
        stats::quantile(
          values,
          probs = 0.25,
          names = FALSE
        )
      ),
      IQR_high = unname(
        stats::quantile(
          values,
          probs = 0.75,
          names = FALSE
        )
      ),
      minimum = min(values),
      maximum = max(values),
      stringsAsFactors = FALSE
    )

    counter <- counter + 1L
  }
}

group_occupancy_summary <- do.call(
  rbind,
  group_summary_rows
)

utils::write.csv(
  group_occupancy_summary,
  file.path(
    OUTPUT_DIR,
    "group_state_occupancy_summary.csv"
  ),
  row.names = FALSE
)

group_occupancy_summary_rounded <-
  group_occupancy_summary

columns_to_round_group <- setdiff(
  names(group_occupancy_summary_rounded),
  c(
    "group",
    "state",
    "n"
  )
)

group_occupancy_summary_rounded[
  columns_to_round_group
] <- lapply(
  group_occupancy_summary_rounded[
    columns_to_round_group
  ],
  function(x) round(x, 4)
)


# ------------------------------------------------------------------------------
# 8. Permutation tests for group differences in occupancy
# ------------------------------------------------------------------------------

permutation_test_mean_difference <- function(
    values,
    group,
    B,
    seed
) {
  values <- as.numeric(
    values
  )

  group <- factor(
    group,
    levels = c(
      "TD",
      "Dyslexia"
    )
  )

  if (
    any(!is.finite(values)) ||
    any(is.na(group))
  ) {
    stop(
      "Invalid values were passed to the permutation test."
    )
  }

  observed_difference <- (
    mean(
      values[
        group == "Dyslexia"
      ]
    ) -
      mean(
        values[
          group == "TD"
        ]
      )
  )

  set.seed(seed)

  permuted_differences <- numeric(
    B
  )

  for (b in seq_len(B)) {
    permuted_group <- sample(
      group,
      size = length(group),
      replace = FALSE
    )

    permuted_differences[b] <- (
      mean(
        values[
          permuted_group == "Dyslexia"
        ]
      ) -
        mean(
          values[
            permuted_group == "TD"
          ]
        )
    )
  }

  p_value <- (
    sum(
      abs(permuted_differences) >=
        abs(observed_difference)
    ) + 1
  ) / (
    B + 1
  )

  list(
    observed_difference = observed_difference,
    p_value = p_value
  )
}

permutation_results <- vector(
  mode = "list",
  length = K
)

for (state in seq_len(K)) {
  permutation_results[[state]] <-
    permutation_test_mean_difference(
      values = occupancy_matrix[, state],
      group = subject_occupancy$group,
      B = N_PERMUTATIONS,
      seed = PERMUTATION_SEED + state
    )
}

occupancy_permutation_table <- data.frame(
  state = seq_len(K),

  mean_TD_percent = vapply(
    seq_len(K),
    function(state) {
      mean(
        occupancy_matrix[
          subject_occupancy$group == "TD",
          state
        ]
      )
    },
    FUN.VALUE = numeric(1)
  ),

  SD_TD = vapply(
    seq_len(K),
    function(state) {
      stats::sd(
        occupancy_matrix[
          subject_occupancy$group == "TD",
          state
        ]
      )
    },
    FUN.VALUE = numeric(1)
  ),

  mean_Dyslexia_percent = vapply(
    seq_len(K),
    function(state) {
      mean(
        occupancy_matrix[
          subject_occupancy$group == "Dyslexia",
          state
        ]
      )
    },
    FUN.VALUE = numeric(1)
  ),

  SD_Dyslexia = vapply(
    seq_len(K),
    function(state) {
      stats::sd(
        occupancy_matrix[
          subject_occupancy$group == "Dyslexia",
          state
        ]
      )
    },
    FUN.VALUE = numeric(1)
  ),

  mean_difference_Dyslexia_minus_TD =
    vapply(
      permutation_results,
      function(result) {
        result$observed_difference
      },
      FUN.VALUE = numeric(1)
    ),

  permutation_p_value =
    vapply(
      permutation_results,
      function(result) {
        result$p_value
      },
      FUN.VALUE = numeric(1)
    ),

  stringsAsFactors = FALSE
)

occupancy_permutation_table$FDR_adjusted_p_value <-
  stats::p.adjust(
    occupancy_permutation_table$permutation_p_value,
    method = "BH"
  )

occupancy_permutation_table$significant_FDR_0_05 <-
  occupancy_permutation_table$FDR_adjusted_p_value < 0.05

utils::write.csv(
  occupancy_permutation_table,
  file.path(
    OUTPUT_DIR,
    "state_occupancy_permutation_tests.csv"
  ),
  row.names = FALSE
)

occupancy_permutation_table_rounded <-
  occupancy_permutation_table

columns_to_round_test <- c(
  "mean_TD_percent",
  "SD_TD",
  "mean_Dyslexia_percent",
  "SD_Dyslexia",
  "mean_difference_Dyslexia_minus_TD",
  "permutation_p_value",
  "FDR_adjusted_p_value"
)

occupancy_permutation_table_rounded[
  columns_to_round_test
] <- lapply(
  occupancy_permutation_table_rounded[
    columns_to_round_test
  ],
  function(x) round(x, 6)
)


# ------------------------------------------------------------------------------
# 9. Construct subject-by-time state matrices
# ------------------------------------------------------------------------------

create_state_matrix <- function(
    selected_subject_ids,
    maximum_time
) {
  selected_subject_ids <- as.integer(
    selected_subject_ids
  )

  state_matrix <- matrix(
    NA_integer_,
    nrow = length(selected_subject_ids),
    ncol = maximum_time,
    dimnames = list(
      as.character(selected_subject_ids),
      as.character(seq_len(maximum_time))
    )
  )

  for (i in seq_along(selected_subject_ids)) {
    subject_id <- selected_subject_ids[i]

    states_i <- subject_sequences[[as.character(subject_id)]]$state

    usable_length <- min(
      length(states_i),
      maximum_time
    )

    state_matrix[
      i,
      seq_len(usable_length)
    ] <- states_i[
      seq_len(usable_length)
    ]
  }

  state_matrix
}

plot_time_limit <- min(
  MAX_TIME_TO_PLOT,
  max(Lengths)
)

all_subject_ids <- subject_information$subject_id

dyslexia_subject_ids <-
  subject_information$subject_id[
    subject_information$group == "Dyslexia"
  ]

TD_subject_ids <-
  subject_information$subject_id[
    subject_information$group == "TD"
  ]

StateMatrix_All <- create_state_matrix(
  all_subject_ids,
  plot_time_limit
)

StateMatrix_Dyslexia <- create_state_matrix(
  dyslexia_subject_ids,
  plot_time_limit
)

StateMatrix_TD <- create_state_matrix(
  TD_subject_ids,
  plot_time_limit
)


# ------------------------------------------------------------------------------
# 10. Plot discrete Viterbi state-sequence heatmaps
# ------------------------------------------------------------------------------

plot_state_sequence_heatmap <- function(
    state_matrix,
    main_title
) {
  state_matrix <- as.matrix(
    state_matrix
  )

  n_rows <- nrow(
    state_matrix
  )

  n_columns <- ncol(
    state_matrix
  )

  if (
    n_rows == 0L ||
    n_columns == 0L
  ) {
    stop(
      "The state matrix supplied to the heatmap is empty."
    )
  }

  image_matrix <- t(
    state_matrix[
      n_rows:1,
      ,
      drop = FALSE
    ]
  )

  graphics::par(
    mar = c(
      5,
      6,
      4,
      8
    ),
    xpd = NA
  )

  graphics::image(
    x = seq_len(n_columns),
    y = seq_len(n_rows),
    z = image_matrix,
    zlim = c(
      0.5,
      K + 0.5
    ),
    col = STATE_COLOURS,
    axes = FALSE,
    xlab = "Time point",
    ylab = "Participant",
    main = main_title,
    useRaster = TRUE
  )

  x_ticks <- unique(
    c(
      1L,
      seq(
        10L,
        n_columns,
        by = 10L
      ),
      n_columns
    )
  )

  graphics::axis(
    side = 1,
    at = x_ticks,
    labels = x_ticks,
    cex.axis = 0.75
  )

  graphics::axis(
    side = 2,
    at = seq_len(n_rows),
    labels = rev(
      rownames(state_matrix)
    ),
    las = 1,
    cex.axis = 0.65
  )

  graphics::box()

  graphics::legend(
    "right",
    inset = c(
      -0.18,
      0
    ),
    legend = paste(
      "State",
      seq_len(K)
    ),
    fill = STATE_COLOURS,
    border = NA,
    bty = "n",
    title = "State"
  )
}

save_sequence_heatmap <- function(
    state_matrix,
    filename_stem,
    main_title
) {
  image_height <- max(
    1600,
    45 * nrow(state_matrix) + 600
  )

  grDevices::png(
    filename = file.path(
      OUTPUT_DIR,
      paste0(
        filename_stem,
        ".png"
      )
    ),
    width = 3000,
    height = image_height,
    res = 300
  )

  plot_state_sequence_heatmap(
    state_matrix,
    main_title
  )

  grDevices::dev.off()

  grDevices::pdf(
    file = file.path(
      OUTPUT_DIR,
      paste0(
        filename_stem,
        ".pdf"
      )
    ),
    width = 12,
    height = max(
      7,
      0.22 * nrow(state_matrix) + 3
    ),
    onefile = TRUE
  )

  plot_state_sequence_heatmap(
    state_matrix,
    main_title
  )

  grDevices::dev.off()
}

save_sequence_heatmap(
  StateMatrix_All,
  "Viterbi_State_Sequences_All_Participants",
  "Viterbi State Sequences for All Participants"
)

save_sequence_heatmap(
  StateMatrix_Dyslexia,
  "Viterbi_State_Sequences_Dyslexia",
  "Viterbi State Sequences for Dyslexia Participants"
)

save_sequence_heatmap(
  StateMatrix_TD,
  "Viterbi_State_Sequences_TD",
  "Viterbi State Sequences for TD Participants"
)


# ------------------------------------------------------------------------------
# 11. Grouped occupancy bar plot
# ------------------------------------------------------------------------------

mean_matrix <- matrix(
  NA_real_,
  nrow = length(group_order),
  ncol = K,
  dimnames = list(
    group_order,
    paste0(
      "State ",
      seq_len(K)
    )
  )
)

SE_matrix <- mean_matrix

for (group_name in group_order) {
  for (state in seq_len(K)) {
    one_row <- group_occupancy_summary[
      group_occupancy_summary$group == group_name &
        group_occupancy_summary$state == state,
      ,
      drop = FALSE
    ]

    if (nrow(one_row) != 1L) {
      stop(
        "Could not identify one unique summary row for ",
        group_name,
        ", state ",
        state,
        "."
      )
    }

    mean_matrix[
      group_name,
      state
    ] <- one_row$mean_percent

    SE_matrix[
      group_name,
      state
    ] <- one_row$standard_error
  }
}

draw_occupancy_barplot <- function() {
  graphics::par(
    mar = c(
      5,
      5,
      4,
      1
    )
  )

  upper_limit <- max(
    mean_matrix + SE_matrix,
    na.rm = TRUE
  ) * 1.15

  bar_positions <- graphics::barplot(
    height = mean_matrix,
    beside = TRUE,
    col = c(
      "grey75",
      "grey35"
    ),
    border = "black",
    ylim = c(
      0,
      upper_limit
    ),
    ylab = "Mean occupancy (%)",
    xlab = "Brain state",
    main = "Mean State Occupancy by Group",
    names.arg = paste(
      "State",
      seq_len(K)
    ),
    legend.text = group_order,
    args.legend = list(
      x = "topright",
      bty = "n",
      title = "Group"
    )
  )

  graphics::arrows(
    x0 = as.vector(bar_positions),
    y0 = as.vector(
      mean_matrix - SE_matrix
    ),
    x1 = as.vector(bar_positions),
    y1 = as.vector(
      mean_matrix + SE_matrix
    ),
    angle = 90,
    code = 3,
    length = 0.04
  )
}

grDevices::png(
  filename = file.path(
    OUTPUT_DIR,
    "Mean_State_Occupancy_by_Group.png"
  ),
  width = 2400,
  height = 1800,
  res = 300
)

draw_occupancy_barplot()

grDevices::dev.off()

grDevices::pdf(
  file = file.path(
    OUTPUT_DIR,
    "Mean_State_Occupancy_by_Group.pdf"
  ),
  width = 8,
  height = 6,
  onefile = TRUE
)

draw_occupancy_barplot()

grDevices::dev.off()


# ------------------------------------------------------------------------------
# 12. Reproducibility information
# ------------------------------------------------------------------------------

analysis_settings <- data.frame(
  setting = c(
    "decoded_state_file",
    "selected_K",
    "number_of_permutations",
    "permutation_seed",
    "permutation_test",
    "multiple_testing_adjustment",
    "maximum_timepoints_displayed",
    "occupancy_definition"
  ),
  value = c(
    normalizePath(
      DECODED_FILE,
      winslash = "/",
      mustWork = TRUE
    ),
    K,
    N_PERMUTATIONS,
    PERMUTATION_SEED,
    paste0(
      "Two-sided permutation test of mean difference; ",
      "Dyslexia minus TD; plus-one p-value correction"
    ),
    "Benjamini-Hochberg FDR across four states",
    plot_time_limit,
    paste0(
      "100 multiplied by the number of time points decoded to a state, ",
      "divided by the participant's total number of decoded time points"
    )
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  analysis_settings,
  file.path(
    OUTPUT_DIR,
    "occupancy_analysis_settings.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "occupancy_analysis_sessionInfo.txt"
  )
)


# ------------------------------------------------------------------------------
# 13. Final output
# ------------------------------------------------------------------------------

cat("\nOccupancy and Viterbi-sequence analysis completed successfully.\n")
cat(
  "Results directory:\n",
  normalizePath(
    OUTPUT_DIR,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)

cat("\nGroup occupancy summary:\n")
print(
  group_occupancy_summary_rounded,
  row.names = FALSE
)

cat("\nPermutation-test results:\n")
print(
  occupancy_permutation_table_rounded,
  row.names = FALSE
)
