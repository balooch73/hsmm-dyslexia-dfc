# ==============================================================================
# LOAD THE FINAL BEST HSMM MODEL AND PREPARE OBJECTS FOR DOWNSTREAM RESULTS
#
# This script does NOT refit the HSMM.
# It only loads:
#   1) FINAL_SELECTED_HSMM_K04.rds
#   2) decoded_states_final_K04.csv
#
# It then creates the objects needed by the older result-analysis codes:
#   final_fit, testrun, yhat, Lengths, Group, patient_list, control_list,
#   StateSeqDF, subject_sequences, sigma_list, mu_list, K
# ==============================================================================

rm(list = ls())
invisible(gc())

suppressPackageStartupMessages(
  library(mhsmm)
)

# ------------------------------------------------------------------------------
# 1. Project paths
# ------------------------------------------------------------------------------

# Run this script from the repository root. The root can alternatively be set
# through the HSMM_PROJECT_ROOT environment variable.
PROJECT_ROOT <- normalizePath(
  Sys.getenv("HSMM_PROJECT_ROOT", unset = "."),
  winslash = "/",
  mustWork = TRUE
)
BASE_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "HSMM_reviewer_reproducible_analysis"
)

MODEL_FILE <- file.path(
  BASE_DIR,
  "FINAL_SELECTED_HSMM_K04.rds"
)

DECODED_FILE <- file.path(
  BASE_DIR,
  "decoded_states_final_K04.csv"
)

if (!file.exists(MODEL_FILE)) {
  stop(
    "The final fitted model was not found:\n",
    MODEL_FILE
  )
}

if (!file.exists(DECODED_FILE)) {
  stop(
    "The final decoded-state file was not found:\n",
    DECODED_FILE
  )
}


# ------------------------------------------------------------------------------
# 2. Load the final best model
# ------------------------------------------------------------------------------

final_fit <- readRDS(
  MODEL_FILE
)

if (
  is.null(final_fit$model) ||
  is.null(final_fit$model$parms.emission) ||
  is.null(final_fit$model$parms.emission$mu) ||
  is.null(final_fit$model$parms.emission$sigma)
) {
  stop(
    "The loaded RDS file does not contain the expected fitted HSMM structure."
  )
}

K <- length(
  final_fit$model$parms.emission$sigma
)

if (K != 4L) {
  stop(
    "The loaded model contains ",
    K,
    " states rather than the selected K = 4 model."
  )
}

mu_list <- final_fit$model$parms.emission$mu
sigma_list <- final_fit$model$parms.emission$sigma


# ------------------------------------------------------------------------------
# 3. Load the final Viterbi-decoded sequences
# ------------------------------------------------------------------------------

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
    paste(
      missing_columns,
      collapse = ", "
    )
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
    "The decoded-state file contains labels outside states 1 to ",
    K,
    "."
  )
}


# ------------------------------------------------------------------------------
# 4. Standardize group labels
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

unexpected_groups <- setdiff(
  unique(decoded_states$group),
  c(
    "TD",
    "Dyslexia"
  )
)

if (length(unexpected_groups) > 0L) {
  stop(
    "Unexpected group label(s): ",
    paste(
      unexpected_groups,
      collapse = ", "
    )
  )
}


# ------------------------------------------------------------------------------
# 5. Sort the decoded observations
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


# ------------------------------------------------------------------------------
# 6. Derive subject information and sequence lengths
# ------------------------------------------------------------------------------

subject_information_rows <- vector(
  mode = "list",
  length = length(subject_ids)
)

subject_sequences <- vector(
  mode = "list",
  length = length(subject_ids)
)

names(subject_sequences) <- as.character(
  subject_ids
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

  if (length(unique(one_subject$file)) != 1L) {
    stop(
      "More than one filename was found for subject ",
      subject_id,
      "."
    )
  }

  subject_sequences[[i]] <- one_subject

  subject_information_rows[[i]] <- data.frame(
    subject_index = i,
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

if (sum(Lengths) != nrow(decoded_states)) {
  stop(
    "The sequence lengths do not sum to the total number of decoded observations."
  )
}


# ------------------------------------------------------------------------------
# 7. Create group objects used in the previous result codes
# ------------------------------------------------------------------------------

Group <- factor(
  subject_information$group,
  levels = c(
    "TD",
    "Dyslexia"
  )
)

patient_list <- subject_information$file[
  Group == "Dyslexia"
]

control_list <- subject_information$file[
  Group == "TD"
]

n_patient <- length(
  patient_list
)

n_control <- length(
  control_list
)


# ------------------------------------------------------------------------------
# 8. Create the concatenated final state sequence
# ------------------------------------------------------------------------------

yhat <- as.integer(
  decoded_states$state
)

if (length(yhat) != sum(Lengths)) {
  stop(
    "The concatenated Viterbi sequence has an unexpected length."
  )
}


# ------------------------------------------------------------------------------
# 9. Compatibility object for the older codes
# ------------------------------------------------------------------------------

# The old scripts used an object called "testrun".
# This alias allows those scripts to access the final selected model.
testrun <- final_fit

# The old scripts also used testrun$yhat.
testrun$yhat <- yhat


# ------------------------------------------------------------------------------
# 10. Create the subject-level state-sequence table
# ------------------------------------------------------------------------------

subject_index_for_each_row <- match(
  decoded_states$subject_id,
  subject_ids
)

StateSeqDF <- data.frame(
  subject = subject_index_for_each_row,
  subject_id = decoded_states$subject_id,
  group = decoded_states$group,
  file = decoded_states$file,
  time = decoded_states$timepoint,
  state = decoded_states$state,
  stringsAsFactors = FALSE
)

if (nrow(StateSeqDF) != sum(Lengths)) {
  stop(
    "StateSeqDF has an unexpected number of rows."
  )
}


# ------------------------------------------------------------------------------
# 11. Optional direct access to fitted model components
# ------------------------------------------------------------------------------

transition_matrix <- final_fit$model$transition
sojourn_distribution <- final_fit$model$sojourn$d
initial_probabilities <- final_fit$model$J


# ------------------------------------------------------------------------------
# 12. Save the reconstructed subject-information table
# ------------------------------------------------------------------------------

utils::write.csv(
  subject_information,
  file.path(
    BASE_DIR,
    "subject_information_reconstructed.csv"
  ),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 13. Final checks and summary
# ------------------------------------------------------------------------------

cat("\n============================================================\n")
cat("FINAL BEST HSMM IS READY FOR DOWNSTREAM ANALYSES\n")
cat("============================================================\n")
cat("Selected number of states:", K, "\n")
cat("Number of participants:", nrow(subject_information), "\n")
cat("Dyslexia participants:", n_patient, "\n")
cat("TD participants:", n_control, "\n")
cat("Total decoded observations:", length(yhat), "\n")
cat("Minimum sequence length:", min(Lengths), "\n")
cat("Maximum sequence length:", max(Lengths), "\n")
cat(
  "Final model file:\n",
  normalizePath(
    MODEL_FILE,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
cat(
  "Decoded-state file:\n",
  normalizePath(
    DECODED_FILE,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n"
)
cat("============================================================\n\n")

cat("Objects available for the next analyses:\n")
cat(
  paste(
    c(
      "final_fit",
      "testrun",
      "yhat",
      "Lengths",
      "Group",
      "patient_list",
      "control_list",
      "n_patient",
      "n_control",
      "subject_information",
      "subject_sequences",
      "StateSeqDF",
      "K",
      "mu_list",
      "sigma_list",
      "transition_matrix",
      "sojourn_distribution",
      "initial_probabilities"
    ),
    collapse = ", "
  ),
  "\n"
)
