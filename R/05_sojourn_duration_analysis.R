# ==============================================================================
# FINAL HSMM SOJOURN-DURATION ANALYSIS
# Reviewer-shareable reproducible script
#
# Input:
#   decoded_states_final_K04.csv
#
# This script:
#   1. reconstructs each participant's final Viterbi sequence;
#   2. extracts complete (uncensored) state runs;
#   3. removes the first and last run of every participant;
#   4. summarizes sojourn durations by group and state;
#   5. compares group-specific duration distributions using
#      D_KL(Dyslexia || TD);
#   6. obtains p-values by permuting participant-level group labels;
#   7. applies Benjamini-Hochberg FDR correction only to estimable states.
#
# Important:
#   If a state has no complete runs in either group, its between-group sojourn
#   distribution is not estimable. Such a state is reported as NA and is not
#   included in the FDR correction.
#
# Required packages:
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
  "Updated_Sojourn_Results"
)

dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

K <- 4L

N_PERMUTATIONS <- 5000L
PERMUTATION_SEED <- 20260712L

# Maximum represented duration. Durations above this value are placed in the
# final bin. With 200 time points per participant, 200 is the natural upper cap.
MAX_DURATION_CAP <- 200L

# Small pseudocount used only after confirming that both groups contain runs.
ALPHA <- 1e-3

# Minimum evidence required in each group before testing a state.
MIN_RUNS_PER_GROUP <- 1L
MIN_PARTICIPANTS_WITH_RUNS <- 1L

# Optional figures.
CREATE_FIGURES <- TRUE


# ------------------------------------------------------------------------------
# 2. Read and validate final decoded-state data
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

  output <- original

  output[
    lower %in% c(
      "control",
      "controls",
      "ctrl",
      "td",
      "typically developing",
      "typically_developing"
    )
  ] <- "TD"

  output[
    lower %in% c(
      "patient",
      "patients",
      "dyslexia",
      "dyslexic"
    )
  ] <- "Dyslexia"

  output
}

decoded_states$group <- standardize_group(
  decoded_states$group
)

GROUP_LEVELS <- c(
  "TD",
  "Dyslexia"
)

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
# 4. Sort data and reconstruct subject sequences
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

  if (
    !identical(
      one_subject$timepoint,
      seq_len(nrow(one_subject))
    )
  ) {
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

  subject_sequences[[i]] <- as.integer(
    one_subject$state
  )

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

subject_information$group <- factor(
  subject_information$group,
  levels = GROUP_LEVELS
)

n_subjects <- nrow(
  subject_information
)


# ------------------------------------------------------------------------------
# 5. Extract complete, uncensored sojourn runs
# ------------------------------------------------------------------------------

extract_complete_sojourns <- function(state_sequence) {
  state_sequence <- as.integer(
    state_sequence
  )

  state_sequence <- state_sequence[
    !is.na(state_sequence)
  ]

  if (length(state_sequence) == 0L) {
    return(
      data.frame(
        state = integer(0),
        duration = integer(0),
        original_run_index = integer(0)
      )
    )
  }

  runs <- rle(
    state_sequence
  )

  run_table <- data.frame(
    state = as.integer(runs$values),
    duration = as.integer(runs$lengths),
    original_run_index = seq_along(runs$values),
    stringsAsFactors = FALSE
  )

  # The first and last runs are boundary-censored.
  # If there are only one or two runs, no complete internal run remains.
  if (nrow(run_table) <= 2L) {
    return(
      run_table[
        0,
        ,
        drop = FALSE
      ]
    )
  }

  run_table[
    2:(nrow(run_table) - 1L),
    ,
    drop = FALSE
  ]
}


sojourn_rows <- list()
row_counter <- 1L

for (i in seq_len(n_subjects)) {
  one_subject_runs <- extract_complete_sojourns(
    subject_sequences[[i]]
  )

  if (nrow(one_subject_runs) == 0L) {
    next
  }

  one_subject_runs$subject_id <-
    subject_information$subject_id[i]

  one_subject_runs$group <-
    as.character(subject_information$group[i])

  one_subject_runs$file <-
    subject_information$file[i]

  one_subject_runs$complete_run_index <-
    seq_len(nrow(one_subject_runs))

  sojourn_rows[[row_counter]] <- one_subject_runs[
    ,
    c(
      "subject_id",
      "group",
      "file",
      "state",
      "duration",
      "original_run_index",
      "complete_run_index"
    ),
    drop = FALSE
  ]

  row_counter <- row_counter + 1L
}

if (length(sojourn_rows) == 0L) {
  stop(
    "No complete internal sojourn runs remained after removing boundary runs."
  )
}

sojourn_long <- do.call(
  rbind,
  sojourn_rows
)

rownames(sojourn_long) <- NULL

utils::write.csv(
  sojourn_long,
  file.path(
    OUTPUT_DIR,
    "complete_subject_level_sojourn_runs.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 6. Summarize sojourn durations by group and state
# ------------------------------------------------------------------------------

safe_summary <- function(values) {
  values <- as.numeric(
    values
  )

  values <- values[
    is.finite(values) &
      values > 0
  ]

  if (length(values) == 0L) {
    return(
      c(
        n_runs = 0,
        mean = NA_real_,
        SD = NA_real_,
        median = NA_real_,
        IQR_low = NA_real_,
        IQR_high = NA_real_,
        minimum = NA_real_,
        maximum = NA_real_
      )
    )
  }

  c(
    n_runs = length(values),
    mean = mean(values),
    SD = if (length(values) > 1L) stats::sd(values) else NA_real_,
    median = stats::median(values),
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
    maximum = max(values)
  )
}


summary_rows <- list()
summary_counter <- 1L

for (group_name in GROUP_LEVELS) {
  for (state in seq_len(K)) {
    state_group_runs <- sojourn_long[
      sojourn_long$group == group_name &
        sojourn_long$state == state,
      ,
      drop = FALSE
    ]

    duration_summary <- safe_summary(
      state_group_runs$duration
    )

    summary_rows[[summary_counter]] <- data.frame(
      group = group_name,
      state = state,
      n_total_participants =
        sum(subject_information$group == group_name),
      n_participants_with_complete_runs =
        length(unique(state_group_runs$subject_id)),
      n_complete_runs =
        as.integer(duration_summary["n_runs"]),
      mean_duration =
        as.numeric(duration_summary["mean"]),
      standard_deviation =
        as.numeric(duration_summary["SD"]),
      median_duration =
        as.numeric(duration_summary["median"]),
      IQR_low =
        as.numeric(duration_summary["IQR_low"]),
      IQR_high =
        as.numeric(duration_summary["IQR_high"]),
      minimum =
        as.numeric(duration_summary["minimum"]),
      maximum =
        as.numeric(duration_summary["maximum"]),
      stringsAsFactors = FALSE
    )

    summary_counter <- summary_counter + 1L
  }
}

sojourn_group_summary <- do.call(
  rbind,
  summary_rows
)

utils::write.csv(
  sojourn_group_summary,
  file.path(
    OUTPUT_DIR,
    "sojourn_duration_group_summary.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 7. Functions for discrete duration distributions and KL divergence
# ------------------------------------------------------------------------------

get_duration_probability <- function(
    durations,
    bins,
    alpha
) {
  durations <- as.integer(
    durations
  )

  durations <- durations[
    is.finite(durations) &
      durations > 0
  ]

  if (length(durations) == 0L) {
    return(NULL)
  }

  capped_durations <- pmin(
    durations,
    max(bins)
  )

  counts <- tabulate(
    capped_durations,
    nbins = length(bins)
  )

  probabilities <- (
    counts + alpha
  ) / (
    sum(counts) +
      alpha * length(bins)
  )

  as.numeric(
    probabilities
  )
}


kl_divergence <- function(
    p,
    q
) {
  if (
    is.null(p) ||
    is.null(q) ||
    length(p) != length(q)
  ) {
    return(NA_real_)
  }

  if (
    any(!is.finite(p)) ||
    any(!is.finite(q)) ||
    any(p <= 0) ||
    any(q <= 0)
  ) {
    return(NA_real_)
  }

  p <- p / sum(p)
  q <- q / sum(q)

  sum(
    p * log(p / q)
  )
}


maximum_observed_duration <- max(
  sojourn_long$duration,
  na.rm = TRUE
)

maximum_duration <- min(
  MAX_DURATION_CAP,
  maximum_observed_duration
)

if (
  !is.finite(maximum_duration) ||
  maximum_duration < 1L
) {
  stop(
    "No valid positive sojourn durations were found."
  )
}

duration_bins <- seq_len(
  maximum_duration
)


# ------------------------------------------------------------------------------
# 8. Observed state-specific KL divergences
# ------------------------------------------------------------------------------

state_is_estimable <- logical(
  K
)

observed_KL <- rep(
  NA_real_,
  K
)

TD_run_counts <- integer(
  K
)

Dyslexia_run_counts <- integer(
  K
)

TD_participant_counts <- integer(
  K
)

Dyslexia_participant_counts <- integer(
  K
)

for (state in seq_len(K)) {
  TD_rows <- sojourn_long[
    sojourn_long$group == "TD" &
      sojourn_long$state == state,
    ,
    drop = FALSE
  ]

  dyslexia_rows <- sojourn_long[
    sojourn_long$group == "Dyslexia" &
      sojourn_long$state == state,
    ,
    drop = FALSE
  ]

  TD_run_counts[state] <- nrow(
    TD_rows
  )

  Dyslexia_run_counts[state] <- nrow(
    dyslexia_rows
  )

  TD_participant_counts[state] <- length(
    unique(TD_rows$subject_id)
  )

  Dyslexia_participant_counts[state] <- length(
    unique(dyslexia_rows$subject_id)
  )

  state_is_estimable[state] <- (
    TD_run_counts[state] >= MIN_RUNS_PER_GROUP &&
      Dyslexia_run_counts[state] >= MIN_RUNS_PER_GROUP &&
      TD_participant_counts[state] >= MIN_PARTICIPANTS_WITH_RUNS &&
      Dyslexia_participant_counts[state] >= MIN_PARTICIPANTS_WITH_RUNS
  )

  if (!state_is_estimable[state]) {
    next
  }

  p_dyslexia <- get_duration_probability(
    dyslexia_rows$duration,
    bins = duration_bins,
    alpha = ALPHA
  )

  q_TD <- get_duration_probability(
    TD_rows$duration,
    bins = duration_bins,
    alpha = ALPHA
  )

  observed_KL[state] <- kl_divergence(
    p = p_dyslexia,
    q = q_TD
  )
}


# ------------------------------------------------------------------------------
# 9. Participant-label permutation test
# ------------------------------------------------------------------------------

set.seed(
  PERMUTATION_SEED
)

permuted_KL <- matrix(
  NA_real_,
  nrow = N_PERMUTATIONS,
  ncol = K,
  dimnames = list(
    NULL,
    paste0(
      "State_",
      seq_len(K)
    )
  )
)

original_groups <- as.character(
  subject_information$group
)

names(original_groups) <- as.character(
  subject_information$subject_id
)

for (b in seq_len(N_PERMUTATIONS)) {
  permuted_groups <- sample(
    original_groups,
    size = length(original_groups),
    replace = FALSE
  )

  names(permuted_groups) <- names(
    original_groups
  )

  run_permuted_group <- permuted_groups[
    as.character(sojourn_long$subject_id)
  ]

  for (state in seq_len(K)) {
    if (!state_is_estimable[state]) {
      next
    }

    durations_TD <- sojourn_long$duration[
      sojourn_long$state == state &
        run_permuted_group == "TD"
    ]

    durations_Dyslexia <- sojourn_long$duration[
      sojourn_long$state == state &
        run_permuted_group == "Dyslexia"
    ]

    if (
      length(durations_TD) < MIN_RUNS_PER_GROUP ||
      length(durations_Dyslexia) < MIN_RUNS_PER_GROUP
    ) {
      next
    }

    participants_TD <- unique(
      sojourn_long$subject_id[
        sojourn_long$state == state &
          run_permuted_group == "TD"
      ]
    )

    participants_Dyslexia <- unique(
      sojourn_long$subject_id[
        sojourn_long$state == state &
          run_permuted_group == "Dyslexia"
      ]
    )

    if (
      length(participants_TD) < MIN_PARTICIPANTS_WITH_RUNS ||
      length(participants_Dyslexia) < MIN_PARTICIPANTS_WITH_RUNS
    ) {
      next
    }

    p_permuted_Dyslexia <- get_duration_probability(
      durations_Dyslexia,
      bins = duration_bins,
      alpha = ALPHA
    )

    q_permuted_TD <- get_duration_probability(
      durations_TD,
      bins = duration_bins,
      alpha = ALPHA
    )

    permuted_KL[b, state] <- kl_divergence(
      p = p_permuted_Dyslexia,
      q = q_permuted_TD
    )
  }
}


# ------------------------------------------------------------------------------
# 10. P-values and FDR correction
# ------------------------------------------------------------------------------

raw_p_values <- rep(
  NA_real_,
  K
)

valid_permutation_counts <- integer(
  K
)

for (state in seq_len(K)) {
  if (!state_is_estimable[state]) {
    next
  }

  valid_null_values <- permuted_KL[
    is.finite(permuted_KL[, state]),
    state
  ]

  valid_permutation_counts[state] <- length(
    valid_null_values
  )

  if (length(valid_null_values) == 0L) {
    next
  }

  # Right-tail test because larger KL indicates greater distributional
  # separation. The plus-one correction prevents p = 0.
  raw_p_values[state] <- (
    1 +
      sum(
        valid_null_values >= observed_KL[state]
      )
  ) / (
    1 +
      length(valid_null_values)
  )
}

FDR_p_values <- rep(
  NA_real_,
  K
)

estimable_p_indices <- which(
  is.finite(raw_p_values)
)

if (length(estimable_p_indices) > 0L) {
  FDR_p_values[
    estimable_p_indices
  ] <- stats::p.adjust(
    raw_p_values[
      estimable_p_indices
    ],
    method = "BH"
  )
}


# ------------------------------------------------------------------------------
# 11. Final inferential results table
# ------------------------------------------------------------------------------

TD_means <- vapply(
  seq_len(K),
  function(state) {
    values <- sojourn_long$duration[
      sojourn_long$group == "TD" &
        sojourn_long$state == state
    ]

    if (length(values) == 0L) {
      return(NA_real_)
    }

    mean(values)
  },
  FUN.VALUE = numeric(1)
)

Dyslexia_means <- vapply(
  seq_len(K),
  function(state) {
    values <- sojourn_long$duration[
      sojourn_long$group == "Dyslexia" &
        sojourn_long$state == state
    ]

    if (length(values) == 0L) {
      return(NA_real_)
    }

    mean(values)
  },
  FUN.VALUE = numeric(1)
)

SojournResults <- data.frame(
  state = seq_len(K),
  n_TD_complete_runs = TD_run_counts,
  n_Dyslexia_complete_runs = Dyslexia_run_counts,
  n_TD_participants_with_runs = TD_participant_counts,
  n_Dyslexia_participants_with_runs =
    Dyslexia_participant_counts,
  mean_TD_duration = TD_means,
  mean_Dyslexia_duration = Dyslexia_means,
  mean_difference_Dyslexia_minus_TD =
    Dyslexia_means - TD_means,
  estimable_between_group_KL = state_is_estimable,
  KL_Dyslexia_to_TD = observed_KL,
  valid_permutations = valid_permutation_counts,
  permutation_p_value = raw_p_values,
  FDR_adjusted_p_value = FDR_p_values,
  significant_FDR_0_05 =
    ifelse(
      is.finite(FDR_p_values),
      FDR_p_values < 0.05,
      NA
    ),
  stringsAsFactors = FALSE
)

SojournResults$interpretation_note <- ifelse(
  SojournResults$estimable_between_group_KL,
  "Between-group duration distribution was estimable.",
  paste0(
    "Not estimable: at least one group lacked sufficient complete ",
    "internal runs for this state."
  )
)

utils::write.csv(
  SojournResults,
  file.path(
    OUTPUT_DIR,
    "sojourn_KL_permutation_results.csv"
  ),
  row.names = FALSE
)

SojournResults_rounded <- SojournResults

numeric_columns_to_round <- c(
  "mean_TD_duration",
  "mean_Dyslexia_duration",
  "mean_difference_Dyslexia_minus_TD",
  "KL_Dyslexia_to_TD",
  "permutation_p_value",
  "FDR_adjusted_p_value"
)

SojournResults_rounded[
  numeric_columns_to_round
] <- lapply(
  SojournResults_rounded[
    numeric_columns_to_round
  ],
  function(x) round(x, 6)
)


# ------------------------------------------------------------------------------
# 12. Optional figures
# ------------------------------------------------------------------------------

if (isTRUE(CREATE_FIGURES)) {
  for (state in seq_len(K)) {
    TD_durations <- sojourn_long$duration[
      sojourn_long$group == "TD" &
        sojourn_long$state == state
    ]

    Dyslexia_durations <- sojourn_long$duration[
      sojourn_long$group == "Dyslexia" &
        sojourn_long$state == state
    ]

    all_durations <- c(
      TD_durations,
      Dyslexia_durations
    )

    if (length(all_durations) > 0L) {
      grDevices::png(
        filename = file.path(
          OUTPUT_DIR,
          sprintf(
            "State_%d_Sojourn_Duration_Histogram.png",
            state
          )
        ),
        width = 2400,
        height = 1800,
        res = 300
      )

      histogram_breaks <- seq(
        0.5,
        max(all_durations) + 0.5,
        by = 1
      )

      if (length(TD_durations) > 0L) {
        graphics::hist(
          TD_durations,
          breaks = histogram_breaks,
          col = grDevices::adjustcolor(
            "grey40",
            alpha.f = 0.45
          ),
          border = "white",
          xlim = c(
            0.5,
            max(all_durations) + 0.5
          ),
          main = paste(
            "State",
            state,
            "Complete Sojourn Durations"
          ),
          xlab = "Duration (time points)",
          ylab = "Number of complete runs"
        )

        if (length(Dyslexia_durations) > 0L) {
          graphics::hist(
            Dyslexia_durations,
            breaks = histogram_breaks,
            col = grDevices::adjustcolor(
              "grey75",
              alpha.f = 0.55
            ),
            border = "white",
            add = TRUE
          )
        }
      } else {
        graphics::hist(
          Dyslexia_durations,
          breaks = histogram_breaks,
          col = grDevices::adjustcolor(
            "grey75",
            alpha.f = 0.55
          ),
          border = "white",
          xlim = c(
            0.5,
            max(all_durations) + 0.5
          ),
          main = paste(
            "State",
            state,
            "Complete Sojourn Durations"
          ),
          xlab = "Duration (time points)",
          ylab = "Number of complete runs"
        )
      }

      graphics::legend(
        "topright",
        legend = c(
          "TD",
          "Dyslexia"
        ),
        fill = c(
          grDevices::adjustcolor(
            "grey40",
            alpha.f = 0.45
          ),
          grDevices::adjustcolor(
            "grey75",
            alpha.f = 0.55
          )
        ),
        border = NA,
        bty = "n"
      )

      grDevices::dev.off()
    }

    if (
      state_is_estimable[state] &&
      any(is.finite(permuted_KL[, state]))
    ) {
      null_values <- permuted_KL[
        is.finite(permuted_KL[, state]),
        state
      ]

      grDevices::png(
        filename = file.path(
          OUTPUT_DIR,
          sprintf(
            "State_%d_KL_Permutation_Null.png",
            state
          )
        ),
        width = 2400,
        height = 1800,
        res = 300
      )

      graphics::hist(
        null_values,
        breaks = 50,
        col = "grey80",
        border = "white",
        main = paste(
          "State",
          state,
          "Permutation Null Distribution"
        ),
        xlab = "KL divergence: Dyslexia || TD",
        ylab = "Frequency"
      )

      graphics::abline(
        v = observed_KL[state],
        lwd = 2,
        lty = 2
      )

      graphics::legend(
        "topright",
        legend = "Observed KL",
        lty = 2,
        lwd = 2,
        bty = "n"
      )

      grDevices::dev.off()
    }
  }
}


# ------------------------------------------------------------------------------
# 13. Reproducibility information
# ------------------------------------------------------------------------------

analysis_settings <- data.frame(
  setting = c(
    "decoded_state_file",
    "selected_K",
    "boundary_run_handling",
    "KL_direction",
    "duration_cap",
    "pseudocount_alpha",
    "number_of_permutations",
    "permutation_seed",
    "permutation_unit",
    "p_value_tail",
    "p_value_correction",
    "multiple_testing_adjustment"
  ),
  value = c(
    normalizePath(
      DECODED_FILE,
      winslash = "/",
      mustWork = TRUE
    ),
    K,
    paste0(
      "First and last run removed for every participant; ",
      "participants with <=2 total runs contribute no complete run."
    ),
    "D_KL(Dyslexia || TD)",
    maximum_duration,
    ALPHA,
    N_PERMUTATIONS,
    PERMUTATION_SEED,
    "Participant-level group labels",
    "Right tail with plus-one correction",
    "(1 + number of null statistics >= observed) / (1 + valid permutations)",
    "Benjamini-Hochberg FDR across estimable states only"
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  analysis_settings,
  file.path(
    OUTPUT_DIR,
    "sojourn_analysis_settings.csv"
  ),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUTPUT_DIR,
    "sojourn_analysis_sessionInfo.txt"
  )
)


# ------------------------------------------------------------------------------
# 14. Final output
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL SOJOURN-DURATION ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Participants:", n_subjects, "\n")
cat("Complete internal runs:", nrow(sojourn_long), "\n")
cat("Duration bins: 1 to", maximum_duration, "\n")
cat("Permutation replicates:", N_PERMUTATIONS, "\n")
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

cat("Group-level duration summary:\n")
print(
  sojourn_group_summary,
  row.names = FALSE
)

cat("\nKL permutation results:\n")
print(
  SojournResults_rounded,
  row.names = FALSE
)
