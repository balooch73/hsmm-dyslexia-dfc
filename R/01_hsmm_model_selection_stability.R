## ============================================================================
## Reproducible HSMM model selection and random-start stability analysis
##
## Model:
##   - Multivariate Gaussian emissions
##   - Nonparametric sojourn distributions
##   - Multiple observation sequences (one sequence per participant)
##
## Workflow:
##   1) Fit K = 2,...,10 with equal numbers of genuinely different starts.
##   2) Select K using the minimum BIC from the best eligible fit for each K.
##   3) Run additional starts for the selected K.
##   4) Select the highest-likelihood eligible fit as the final model.
##   5) Quantify random-start agreement after state matching.
##
## IMPORTANT:
##   - cov.shrink is used ONLY to create numerically stable starting covariance
##     matrices. Final emission parameters are estimated by mstep.mvnorm.
##   - No post-estimation smoothing is applied to the sojourn distributions.
## ============================================================================

rm(list = ls())
gc()

## ----------------------------------------------------------------------------
## 1. Packages
## ----------------------------------------------------------------------------

required_packages <- c("readxl", "mhsmm", "mvtnorm", "corpcor", "clue")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

## ----------------------------------------------------------------------------
## 2. User settings
## ----------------------------------------------------------------------------

# Run this script from the repository root. The root can alternatively be set
# through the HSMM_PROJECT_ROOT environment variable.
PROJECT_ROOT <- normalizePath(
  Sys.getenv("HSMM_PROJECT_ROOT", unset = "."),
  winslash = "/",
  mustWork = TRUE
)
DATA_DIR <- file.path(PROJECT_ROOT, "data")
PATIENT_DIR <- file.path(DATA_DIR, "Patients")
CONTROL_DIR <- file.path(DATA_DIR, "Controls")
OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "HSMM_reviewer_reproducible_analysis"
)

## Candidate number of states
K_VALUES <- 2:10

## Equal number of starts used for BIC model selection
N_STARTS_MODEL_SELECTION <- 10L

## Total number of attempted starts for the selected K
N_STARTS_SELECTED_K <- 50L

## EM settings
MAXIT <- 800L
M_CAP <- 200L

## Reproducibility
MASTER_SEED <- 20260710L

## Do not mix new results with an older run
if (dir.exists(OUTPUT_DIR)) {
  stop(
    "The output directory already exists. Rename or delete it before rerunning: ",
    OUTPUT_DIR
  )
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "fits"), recursive = TRUE, showWarnings = FALSE)

set.seed(MASTER_SEED)

## ----------------------------------------------------------------------------
## 3. Read and validate data
## ----------------------------------------------------------------------------

patient_files <- sort(list.files(
  PATIENT_DIR,
  pattern = "\\.xlsx$",
  full.names = TRUE,
  ignore.case = TRUE
))

control_files <- sort(list.files(
  CONTROL_DIR,
  pattern = "\\.xlsx$",
  full.names = TRUE,
  ignore.case = TRUE
))

if (length(patient_files) == 0L) stop("No .xlsx files found in Patients.")
if (length(control_files) == 0L) stop("No .xlsx files found in Controls.")

all_files <- c(patient_files, control_files)
group_labels <- c(
  rep("Dyslexia", length(patient_files)),
  rep("TD", length(control_files))
)

read_subject <- function(file) {
  dat <- readxl::read_excel(file, .name_repair = "minimal")

  numeric_columns <- vapply(dat, is.numeric, FUN.VALUE = logical(1))
  if (!all(numeric_columns)) {
    stop(
      "Non-numeric column(s) in ", basename(file), ": ",
      paste(names(dat)[!numeric_columns], collapse = ", ")
    )
  }

  x <- as.matrix(dat)
  storage.mode(x) <- "double"

  if (nrow(x) < 2L) stop("Too few time points in ", basename(file), ".")
  if (ncol(x) < 2L) stop("Too few variables in ", basename(file), ".")
  if (any(!is.finite(x))) stop("Missing/non-finite values in ", basename(file), ".")

  x
}

all_list <- lapply(all_files, read_subject)

p_values <- vapply(all_list, ncol, FUN.VALUE = integer(1))
if (length(unique(p_values)) != 1L) {
  stop("The number of variables differs across participants.")
}

p <- p_values[1L]
reference_names <- colnames(all_list[[1L]])

same_names <- vapply(
  all_list,
  function(x) identical(colnames(x), reference_names),
  FUN.VALUE = logical(1)
)

if (!all(same_names)) stop("Variable names/order differ across participants.")

Lengths <- as.integer(vapply(all_list, nrow, FUN.VALUE = integer(1)))
N_total <- sum(Lengths)
n_subjects <- length(all_list)

## Common maximum duration support for every candidate K
M <- min(M_CAP, max(Lengths))
if (M < 2L) stop("M must be at least 2.")

## Subject-wise z-scoring
zscore_subject <- function(x, file) {
  s <- apply(x, 2L, stats::sd)

  if (any(!is.finite(s)) || any(s <= sqrt(.Machine$double.eps))) {
    stop("Zero/invalid within-subject variance in ", basename(file), ".")
  }

  z <- scale(x, center = TRUE, scale = TRUE)
  z <- as.matrix(z)
  storage.mode(z) <- "double"

  if (any(!is.finite(z))) stop("Invalid z-scores in ", basename(file), ".")
  z
}

scaled_list <- Map(zscore_subject, all_list, all_files)
DataFinal <- do.call(rbind, scaled_list)
storage.mode(DataFinal) <- "double"

Datanew <- list(x = DataFinal, N = Lengths)
class(Datanew) <- "hsmm.data"

subject_start_indices <- c(1L, head(cumsum(Lengths), -1L) + 1L)

subject_information <- data.frame(
  subject_id = seq_along(all_files),
  group = group_labels,
  file = basename(all_files),
  n_timepoints = Lengths,
  stringsAsFactors = FALSE
)

write.csv(
  subject_information,
  file.path(OUTPUT_DIR, "subject_information.csv"),
  row.names = FALSE
)

cat("Subjects:", n_subjects, "\n")
cat("Variables:", p, "\n")
cat("Total observations:", N_total, "\n")
cat("Nonparametric duration support M:", M, "\n")

## ----------------------------------------------------------------------------
## 4. Helper functions
## ----------------------------------------------------------------------------

normalize_probability <- function(x, floor_value = 1e-10) {
  x <- as.numeric(x)
  x[!is.finite(x) | x < floor_value] <- floor_value
  x / sum(x)
}

random_probability <- function(n, alpha = 1) {
  normalize_probability(stats::rgamma(n, shape = alpha, rate = 1))
}

plain_numeric_matrix <- function(x, nrow_value, ncol_value) {
  out <- matrix(as.numeric(x), nrow = nrow_value, ncol = ncol_value)
  storage.mode(out) <- "double"
  out
}

make_positive_definite <- function(S, eigen_floor = 1e-6) {
  S <- plain_numeric_matrix(S, p, p)
  S <- (S + t(S)) / 2

  eig_min <- min(eigen(S, symmetric = TRUE, only.values = TRUE)$values)
  if (!is.finite(eig_min)) stop("Invalid covariance eigenvalues.")

  if (eig_min < eigen_floor) {
    S <- S + diag(eigen_floor - eig_min + 1e-8, p)
  }

  S <- (S + t(S)) / 2
  storage.mode(S) <- "double"
  S
}

## Shrinkage covariance used only for initialization
safe_shrink_covariance <- function(X, fallback) {
  S <- tryCatch(
    corpcor::cov.shrink(X, verbose = FALSE),
    error = function(e) fallback
  )

  S <- plain_numeric_matrix(S, p, p)  # removes class "shrinkage"
  make_positive_definite(S)
}

pooled_covariance <- corpcor::cov.shrink(DataFinal, verbose = FALSE)
pooled_covariance <- make_positive_definite(pooled_covariance)

## ----------------------------------------------------------------------------
## 5. Create one genuinely randomized starting model
## ----------------------------------------------------------------------------

create_starting_model <- function(K, seed) {
  set.seed(seed)

  minimum_cluster_size <- max(5L, p + 2L)
  km <- NULL

  for (attempt in seq_len(100L)) {
    center_rows <- sample.int(nrow(DataFinal), K, replace = FALSE)

    candidate <- tryCatch(
      stats::kmeans(
        DataFinal,
        centers = DataFinal[center_rows, , drop = FALSE],
        iter.max = 100L,
        nstart = 1L,
        algorithm = "Lloyd"
      ),
      error = function(e) NULL
    )

    if (is.null(candidate)) next

    sizes <- tabulate(candidate$cluster, nbins = K)
    if (all(sizes >= minimum_cluster_size)) {
      km <- candidate
      break
    }
  }

  if (is.null(km)) stop("Could not obtain a valid randomized k-means partition.")

  clusters <- as.integer(km$cluster)
  mu_initial <- vector("list", K)
  sigma_initial <- vector("list", K)

  for (state in seq_len(K)) {
    X_state <- DataFinal[clusters == state, , drop = FALSE]

    ## Different state means across starts
    mu_initial[[state]] <- as.numeric(
      colMeans(X_state) + stats::rnorm(p, mean = 0, sd = 0.02)
    )

    ## Positive-definite shrinkage covariance, used only as a starting value
    S <- safe_shrink_covariance(X_state, pooled_covariance)

    ## Small positive-definite perturbation makes starts genuinely different
    A <- matrix(stats::rnorm(p * p, mean = 0, sd = 0.003), p, p)
    S <- S + crossprod(A)
    S <- S * exp(stats::rnorm(1L, mean = 0, sd = 0.02))

    sigma_initial[[state]] <- make_positive_definite(S)
  }

  ## Initial-state probabilities are based on participant starting observations
  initial_counts <- tabulate(
    clusters[subject_start_indices],
    nbins = K
  )
  init_initial <- normalize_probability(initial_counts + 0.5)

  ## Random embedded transition matrix; self-transitions are zero in an HSMM
  transition_initial <- matrix(0, K, K)
  for (from_state in seq_len(K)) {
    targets <- setdiff(seq_len(K), from_state)
    transition_initial[from_state, targets] <- random_probability(
      length(targets),
      alpha = 1
    )
  }
  diag(transition_initial) <- 0

  ## Positive nonparametric starting duration distributions
  duration_initial <- matrix(0, nrow = M, ncol = K)
  lambda_max <- max(3, min(40, M / 3))

  for (state in seq_len(K)) {
    lambda_state <- stats::runif(1L, min = 2, max = lambda_max)
    d <- stats::dpois(seq_len(M), lambda = lambda_state)
    d <- 0.95 * normalize_probability(d) + 0.05 * rep(1 / M, M)
    duration_initial[, state] <- normalize_probability(d)
  }

  colnames(duration_initial) <- paste0("State", seq_len(K))

  mhsmm::hsmmspec(
    init = init_initial,
    transition = transition_initial,
    parms.emission = list(mu = mu_initial, sigma = sigma_initial),
    sojourn = list(d = duration_initial, type = "nonparametric"),
    dens.emission = mhsmm::dmvnorm.hsmm
  )
}

## ----------------------------------------------------------------------------
## 6. Validate one fitted model
## ----------------------------------------------------------------------------

validate_fit <- function(fit, K, tolerance = 1e-6) {
  if (is.null(fit$model) || is.null(fit$loglik)) return(FALSE)

  ll <- as.numeric(fit$loglik)
  if (length(ll) < 2L || any(!is.finite(ll))) return(FALSE)

  model <- fit$model

  if (length(model$init) != K || any(!is.finite(model$init))) return(FALSE)
  if (abs(sum(model$init) - 1) > tolerance) return(FALSE)

  P <- as.matrix(model$transition)
  if (!all(dim(P) == c(K, K)) || any(!is.finite(P))) return(FALSE)
  if (max(abs(diag(P))) > tolerance) return(FALSE)
  if (max(abs(rowSums(P) - 1)) > tolerance) return(FALSE)

  d <- as.matrix(model$sojourn$d)
  if (!all(dim(d) == c(M, K)) || any(!is.finite(d)) || any(d < 0)) return(FALSE)
  if (max(abs(colSums(d) - 1)) > tolerance) return(FALSE)

  mu <- model$parms.emission$mu
  Sigma <- model$parms.emission$sigma

  if (length(mu) != K || length(Sigma) != K) return(FALSE)

  for (state in seq_len(K)) {
    if (length(mu[[state]]) != p || any(!is.finite(mu[[state]]))) return(FALSE)

    S <- as.matrix(Sigma[[state]])
    if (!all(dim(S) == c(p, p)) || any(!is.finite(S))) return(FALSE)

    eig_min <- min(eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values)
    if (!is.finite(eig_min) || eig_min <= 0) return(FALSE)
  }

  TRUE
}

## ----------------------------------------------------------------------------
## 7. Fit one random start robustly
## ----------------------------------------------------------------------------

result_row <- function(
    K, start, seed, completed, stopped_before_maxit, valid,
    iterations = NA_integer_, logLik = NA_real_, final_change = NA_real_,
    warnings = "", error = "", fit_file = ""
) {
  data.frame(
    K = as.integer(K),
    start = as.integer(start),
    seed = as.integer(seed),
    completed = as.logical(completed),
    stopped_before_maxit = as.logical(stopped_before_maxit),
    valid = as.logical(valid),
    eligible = as.logical(completed && stopped_before_maxit && valid),
    iterations = as.integer(iterations),
    logLik = as.numeric(logLik),
    final_change = as.numeric(final_change),
    warnings = as.character(warnings),
    error = as.character(error),
    fit_file = as.character(fit_file),
    stringsAsFactors = FALSE
  )
}

fit_one_start <- function(K, start, seed, fit_directory) {
  cat("K =", K, "| start =", start, "| seed =", seed, "\n")

  starting_model <- tryCatch(
    create_starting_model(K, seed),
    error = function(e) e
  )

  if (inherits(starting_model, "error")) {
    return(list(
      fit = NULL,
      result = result_row(
        K, start, seed,
        completed = FALSE,
        stopped_before_maxit = FALSE,
        valid = FALSE,
        error = conditionMessage(starting_model)
      )
    ))
  }

  warning_messages <- character(0)

  fitted <- tryCatch(
    withCallingHandlers(
      mhsmm::hsmmfit(
        x = Datanew,
        model = starting_model,
        mstep = mhsmm::mstep.mvnorm,
        M = M,
        maxit = MAXIT,
        lock.transition = FALSE,
        lock.d = FALSE,
        graphical = FALSE
      ),
      warning = function(w) {
        warning_messages <<- c(warning_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  warning_text <- paste(unique(warning_messages), collapse = " | ")

  if (inherits(fitted, "error")) {
    return(list(
      fit = NULL,
      result = result_row(
        K, start, seed,
        completed = FALSE,
        stopped_before_maxit = FALSE,
        valid = FALSE,
        warnings = warning_text,
        error = conditionMessage(fitted)
      )
    ))
  }

  gamma_failure <- any(grepl(
    "NaNs detected in gamma",
    warning_messages,
    fixed = TRUE
  ))

  if (gamma_failure || is.null(fitted$loglik)) {
    return(list(
      fit = NULL,
      result = result_row(
        K, start, seed,
        completed = FALSE,
        stopped_before_maxit = FALSE,
        valid = FALSE,
        warnings = warning_text,
        error = if (gamma_failure) {
          "Numerical failure: NaNs detected in gamma."
        } else {
          "No log-likelihood trace was returned."
        }
      )
    ))
  }

  ll <- as.numeric(fitted$loglik)
  if (length(ll) < 2L || any(!is.finite(ll))) {
    return(list(
      fit = NULL,
      result = result_row(
        K, start, seed,
        completed = FALSE,
        stopped_before_maxit = FALSE,
        valid = FALSE,
        warnings = warning_text,
        error = "Invalid log-likelihood trace."
      )
    ))
  }

  iterations <- length(ll)
  stopped_before_maxit <- iterations < MAXIT
  valid <- validate_fit(fitted, K)
  eligible <- stopped_before_maxit && valid

  fit_file <- ""
  if (eligible) {
    fit_file <- file.path(
      fit_directory,
      sprintf("K%02d_start%03d_seed%d.rds", K, start, seed)
    )
    saveRDS(fitted, fit_file)
  }

  list(
    fit = if (eligible) fitted else NULL,
    result = result_row(
      K, start, seed,
      completed = TRUE,
      stopped_before_maxit = stopped_before_maxit,
      valid = valid,
      iterations = iterations,
      logLik = tail(ll, 1L),
      final_change = tail(diff(ll), 1L),
      warnings = warning_text,
      error = "",
      fit_file = fit_file
    )
  )
}

## ----------------------------------------------------------------------------
## 8. BIC parameter count
## ----------------------------------------------------------------------------

information_criteria <- function(log_likelihood, K) {
  n_initial <- K - 1L
  n_transition <- K * (K - 2L)
  n_means <- K * p
  n_covariances <- K * p * (p + 1L) / 2
  n_sojourn <- K * (M - 1L)

  n_parameters <- n_initial + n_transition + n_means + n_covariances + n_sojourn

  data.frame(
    n_initial = n_initial,
    n_transition = n_transition,
    n_means = n_means,
    n_covariances = n_covariances,
    n_sojourn = n_sojourn,
    n_parameters = n_parameters,
    AIC = -2 * log_likelihood + 2 * n_parameters,
    BIC = -2 * log_likelihood + log(N_total) * n_parameters
  )
}

## ----------------------------------------------------------------------------
## 9. Model selection: K = 2,...,10, with 10 starts per K
## ----------------------------------------------------------------------------

model_selection_rows <- list()
run_tables <- list()

for (K in K_VALUES) {
  cat("\n============================================================\n")
  cat("MODEL SELECTION: K =", K, "\n")
  cat("============================================================\n")

  fit_directory <- file.path(OUTPUT_DIR, "fits", sprintf("K_%02d", K))
  dir.create(fit_directory, recursive = TRUE, showWarnings = FALSE)

  results_K <- vector("list", N_STARTS_MODEL_SELECTION)

  for (start in seq_len(N_STARTS_MODEL_SELECTION)) {
    seed <- MASTER_SEED + K * 10000L + start
    one <- fit_one_start(K, start, seed, fit_directory)
    results_K[[start]] <- one$result
  }

  table_K <- do.call(rbind, results_K)
  run_tables[[as.character(K)]] <- table_K

  write.csv(
    table_K,
    file.path(OUTPUT_DIR, sprintf("random_starts_model_selection_K%02d.csv", K)),
    row.names = FALSE
  )

  eligible_K <- table_K[table_K$eligible & is.finite(table_K$logLik), , drop = FALSE]

  if (nrow(eligible_K) == 0L) {
    stop("No eligible model was obtained for K = ", K, ".")
  }

  best_index <- which.max(eligible_K$logLik)
  best_row <- eligible_K[best_index, , drop = FALSE]
  ic <- information_criteria(best_row$logLik, K)

  best_fit_copy <- file.path(OUTPUT_DIR, sprintf("best_model_selection_fit_K%02d.rds", K))
  file.copy(best_row$fit_file, best_fit_copy, overwrite = TRUE)

  model_selection_rows[[as.character(K)]] <- data.frame(
    K = K,
    attempted_starts = N_STARTS_MODEL_SELECTION,
    eligible_starts = nrow(eligible_K),
    best_start = best_row$start,
    best_seed = best_row$seed,
    best_logLik = best_row$logLik,
    n_parameters = ic$n_parameters,
    AIC = ic$AIC,
    BIC = ic$BIC,
    best_fit_file = normalizePath(best_fit_copy, winslash = "/", mustWork = TRUE),
    stringsAsFactors = FALSE
  )
}

model_selection <- do.call(rbind, model_selection_rows)
rownames(model_selection) <- NULL
model_selection <- model_selection[order(model_selection$K), ]
model_selection$delta_BIC <- model_selection$BIC - min(model_selection$BIC)
model_selection$BIC_rank <- rank(model_selection$BIC, ties.method = "min")

write.csv(
  model_selection,
  file.path(OUTPUT_DIR, "model_selection_BIC.csv"),
  row.names = FALSE
)

print(model_selection)

selected_row <- which.min(model_selection$BIC)
selected_K <- model_selection$K[selected_row]

cat("\nSelected K by minimum BIC:", selected_K, "\n")

## ----------------------------------------------------------------------------
## 10. Additional random starts for the selected K
## ----------------------------------------------------------------------------

selected_fit_directory <- file.path(OUTPUT_DIR, "fits", sprintf("K_%02d", selected_K))
dir.create(selected_fit_directory, recursive = TRUE, showWarnings = FALSE)

selected_start_results <- run_tables[[as.character(selected_K)]]

if (N_STARTS_SELECTED_K > N_STARTS_MODEL_SELECTION) {
  extra_starts <- (N_STARTS_MODEL_SELECTION + 1L):N_STARTS_SELECTED_K

  extra_results <- vector("list", length(extra_starts))

  for (i in seq_along(extra_starts)) {
    start <- extra_starts[i]
    seed <- MASTER_SEED + selected_K * 10000L + start
    one <- fit_one_start(selected_K, start, seed, selected_fit_directory)
    extra_results[[i]] <- one$result
  }

  selected_start_results <- rbind(
    selected_start_results,
    do.call(rbind, extra_results)
  )
}

write.csv(
  selected_start_results,
  file.path(OUTPUT_DIR, sprintf("random_starts_selected_K%02d.csv", selected_K)),
  row.names = FALSE
)

selected_eligible <- selected_start_results[
  selected_start_results$eligible & is.finite(selected_start_results$logLik),
  ,
  drop = FALSE
]

if (nrow(selected_eligible) < 2L) {
  stop("Fewer than two eligible solutions were obtained for the selected K.")
}

best_final_row <- selected_eligible[which.max(selected_eligible$logLik), , drop = FALSE]
final_fit <- readRDS(best_final_row$fit_file)

final_model_file <- file.path(
  OUTPUT_DIR,
  sprintf("FINAL_SELECTED_HSMM_K%02d.rds", selected_K)
)
saveRDS(final_fit, final_model_file)

final_ic <- information_criteria(best_final_row$logLik, selected_K)
final_model_summary <- data.frame(
  selected_K = selected_K,
  attempted_starts = nrow(selected_start_results),
  eligible_starts = nrow(selected_eligible),
  final_start = best_final_row$start,
  final_seed = best_final_row$seed,
  final_logLik = best_final_row$logLik,
  n_parameters = final_ic$n_parameters,
  AIC = final_ic$AIC,
  BIC = final_ic$BIC,
  final_model_file = normalizePath(final_model_file, winslash = "/", mustWork = TRUE),
  stringsAsFactors = FALSE
)

write.csv(
  final_model_summary,
  file.path(OUTPUT_DIR, "final_model_summary.csv"),
  row.names = FALSE
)

## ----------------------------------------------------------------------------
## 11. Decode the final model
## ----------------------------------------------------------------------------

decoded_final <- predict(final_fit, Datanew, method = "viterbi")
if (is.null(decoded_final$s) || length(decoded_final$s) != N_total) {
  stop("Final Viterbi sequence was not returned correctly.")
}

final_decoded_table <- data.frame(
  subject_id = rep(seq_along(Lengths), times = Lengths),
  group = rep(group_labels, times = Lengths),
  file = rep(basename(all_files), times = Lengths),
  timepoint = sequence(Lengths),
  state = as.integer(decoded_final$s),
  stringsAsFactors = FALSE
)

write.csv(
  final_decoded_table,
  file.path(OUTPUT_DIR, sprintf("decoded_states_final_K%02d.csv", selected_K)),
  row.names = FALSE
)

## ----------------------------------------------------------------------------
## 12. Random-start stability after state matching
## ----------------------------------------------------------------------------

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  choose2 <- function(z) z * (z - 1) / 2

  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(n)

  if (total_pairs == 0) return(NA_real_)

  expected <- sum_rows * sum_cols / total_pairs
  maximum <- 0.5 * (sum_rows + sum_cols)

  if (abs(maximum - expected) < .Machine$double.eps) return(1)
  (sum_cells - expected) / (maximum - expected)
}

upper_vector <- function(S) {
  S <- as.matrix(S)
  S[upper.tri(S)]
}

safe_correlation <- function(x, y) {
  value <- suppressWarnings(stats::cor(x, y, use = "pairwise.complete.obs"))
  if (!is.finite(value)) 0 else value
}

subject_occupancy <- function(states, K) {
  out <- matrix(NA_real_, nrow = length(Lengths), ncol = K)
  start <- 1L

  for (i in seq_along(Lengths)) {
    end <- start + Lengths[i] - 1L
    out[i, ] <- tabulate(states[start:end], nbins = K) / Lengths[i]
    start <- end + 1L
  }

  colnames(out) <- paste0("State", seq_len(K))
  out
}

match_candidate_to_reference <- function(reference_fit, candidate_fit, K) {
  ref_mu <- reference_fit$model$parms.emission$mu
  ref_sigma <- reference_fit$model$parms.emission$sigma
  cand_mu <- candidate_fit$model$parms.emission$mu
  cand_sigma <- candidate_fit$model$parms.emission$sigma

  similarity <- matrix(NA_real_, nrow = K, ncol = K)

  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      covariance_similarity <- safe_correlation(
        upper_vector(ref_sigma[[i]]),
        upper_vector(cand_sigma[[j]])
      )

      mean_similarity <- safe_correlation(
        as.numeric(ref_mu[[i]]),
        as.numeric(cand_mu[[j]])
      )

      similarity[i, j] <-
        0.5 * covariance_similarity +
        0.5 * mean_similarity
    }
  }

  ## solve_LSAP requires a finite nonnegative matrix. Convert the
  ## combined similarity from [-1, 1] to a score in [0, 1] and maximize it.
  similarity[!is.finite(similarity)] <- -1
  similarity <- pmax(-1, pmin(1, similarity))

  score <- (similarity + 1) / 2
  score[!is.finite(score)] <- 0
  score <- matrix(
    pmax(0, pmin(1, score)),
    nrow = K,
    ncol = K
  )

  assignment <- as.integer(
    clue::solve_LSAP(score, maximum = TRUE)
  )

  ## assignment[i] is the candidate state matched to reference state i
  assignment
}

relabel_candidate_states <- function(states, assignment, K) {
  old_to_new <- integer(K)
  old_to_new[assignment] <- seq_len(K)
  old_to_new[states]
}

## Read and decode every eligible selected-K fit once
stability_objects <- lapply(seq_len(nrow(selected_eligible)), function(i) {
  fit <- readRDS(selected_eligible$fit_file[i])
  decoded <- predict(fit, Datanew, method = "viterbi")

  if (is.null(decoded$s) || length(decoded$s) != N_total) {
    stop("Invalid decoded sequence in stability fit: ", selected_eligible$fit_file[i])
  }

  list(
    fit = fit,
    states = as.integer(decoded$s),
    start = selected_eligible$start[i],
    seed = selected_eligible$seed[i],
    logLik = selected_eligible$logLik[i]
  )
})

pairwise_rows <- list()
counter <- 1L

for (a in seq_len(length(stability_objects) - 1L)) {
  for (b in (a + 1L):length(stability_objects)) {
    A <- stability_objects[[a]]
    B <- stability_objects[[b]]

    assignment <- match_candidate_to_reference(A$fit, B$fit, selected_K)
    B_states_aligned <- relabel_candidate_states(B$states, assignment, selected_K)

    covariance_correlations <- numeric(selected_K)
    mean_correlations <- numeric(selected_K)

    for (state in seq_len(selected_K)) {
      candidate_state <- assignment[state]

      covariance_correlations[state] <- safe_correlation(
        upper_vector(A$fit$model$parms.emission$sigma[[state]]),
        upper_vector(B$fit$model$parms.emission$sigma[[candidate_state]])
      )

      mean_correlations[state] <- safe_correlation(
        as.numeric(A$fit$model$parms.emission$mu[[state]]),
        as.numeric(B$fit$model$parms.emission$mu[[candidate_state]])
      )
    }

    occupancy_A <- subject_occupancy(A$states, selected_K)
    occupancy_B <- subject_occupancy(B_states_aligned, selected_K)

    pairwise_rows[[counter]] <- data.frame(
      run_A = a,
      run_B = b,
      start_A = A$start,
      start_B = B$start,
      seed_A = A$seed,
      seed_B = B$seed,
      ARI = adjusted_rand_index(A$states, B_states_aligned),
      covariance_correlation_median = median(covariance_correlations),
      covariance_correlation_minimum = min(covariance_correlations),
      mean_vector_correlation_median = median(mean_correlations),
      mean_vector_correlation_minimum = min(mean_correlations),
      occupancy_correlation = safe_correlation(
        as.vector(occupancy_A),
        as.vector(occupancy_B)
      ),
      stringsAsFactors = FALSE
    )

    counter <- counter + 1L
  }
}

pairwise_stability <- do.call(rbind, pairwise_rows)

write.csv(
  pairwise_stability,
  file.path(OUTPUT_DIR, sprintf("pairwise_stability_selected_K%02d.csv", selected_K)),
  row.names = FALSE
)

summarize_metric <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) == 0L) {
    return(
      c(
        median = NA_real_,
        IQR_low = NA_real_,
        IQR_high = NA_real_,
        minimum = NA_real_,
        maximum = NA_real_
      )
    )
  }

  c(
    median = stats::median(x),
    IQR_low = unname(stats::quantile(x, 0.25)),
    IQR_high = unname(stats::quantile(x, 0.75)),
    minimum = min(x),
    maximum = max(x)
  )
}

metric_names <- c(
  "ARI",
  "covariance_correlation_median",
  "covariance_correlation_minimum",
  "mean_vector_correlation_median",
  "mean_vector_correlation_minimum",
  "occupancy_correlation"
)

stability_summary <- do.call(rbind, lapply(metric_names, function(metric) {
  values <- summarize_metric(pairwise_stability[[metric]])
  data.frame(
    metric = metric,
    median = unname(values["median"]),
    IQR_low = unname(values["IQR_low"]),
    IQR_high = unname(values["IQR_high"]),
    minimum = unname(values["minimum"]),
    maximum = unname(values["maximum"]),
    stringsAsFactors = FALSE
  )
}))

write.csv(
  stability_summary,
  file.path(OUTPUT_DIR, sprintf("stability_summary_selected_K%02d.csv", selected_K)),
  row.names = FALSE
)

print(stability_summary)

## ----------------------------------------------------------------------------
## 13. Similarity among states in the final fitted model
## ----------------------------------------------------------------------------

final_state_similarity_rows <- list()
counter <- 1L

for (i in seq_len(selected_K - 1L)) {
  for (j in (i + 1L):selected_K) {
    final_state_similarity_rows[[counter]] <- data.frame(
      state_1 = i,
      state_2 = j,
      covariance_pattern_correlation = safe_correlation(
        upper_vector(final_fit$model$parms.emission$sigma[[i]]),
        upper_vector(final_fit$model$parms.emission$sigma[[j]])
      ),
      mean_vector_correlation = safe_correlation(
        as.numeric(final_fit$model$parms.emission$mu[[i]]),
        as.numeric(final_fit$model$parms.emission$mu[[j]])
      ),
      stringsAsFactors = FALSE
    )
    counter <- counter + 1L
  }
}

final_state_similarity <- do.call(rbind, final_state_similarity_rows)

write.csv(
  final_state_similarity,
  file.path(OUTPUT_DIR, sprintf("within_model_state_similarity_K%02d.csv", selected_K)),
  row.names = FALSE
)

## ----------------------------------------------------------------------------
## 14. Figures and reproducibility information
## ----------------------------------------------------------------------------

png(
  file.path(OUTPUT_DIR, "BIC_by_K.png"),
  width = 1600,
  height = 1200,
  res = 200
)
plot(
  model_selection$K,
  model_selection$BIC,
  type = "b",
  pch = 16,
  xlab = "Number of states (K)",
  ylab = "BIC",
  main = "HSMM model selection"
)
abline(v = selected_K, lty = 2)
dev.off()

png(
  file.path(OUTPUT_DIR, "best_loglikelihood_by_K.png"),
  width = 1600,
  height = 1200,
  res = 200
)
plot(
  model_selection$K,
  model_selection$best_logLik,
  type = "b",
  pch = 16,
  xlab = "Number of states (K)",
  ylab = "Best eligible log-likelihood",
  main = "HSMM log-likelihood by K"
)
dev.off()

analysis_settings <- data.frame(
  setting = c(
    "K_values",
    "model_selection_starts_per_K",
    "selected_K_total_attempted_starts",
    "maximum_EM_iterations",
    "master_seed",
    "number_of_subjects",
    "number_of_variables",
    "total_observations",
    "maximum_sojourn_duration_M",
    "selected_K"
  ),
  value = c(
    paste(K_VALUES, collapse = ", "),
    N_STARTS_MODEL_SELECTION,
    N_STARTS_SELECTED_K,
    MAXIT,
    MASTER_SEED,
    n_subjects,
    p,
    N_total,
    M,
    selected_K
  ),
  stringsAsFactors = FALSE
)

write.csv(
  analysis_settings,
  file.path(OUTPUT_DIR, "analysis_settings.csv"),
  row.names = FALSE
)

package_versions <- data.frame(
  package = required_packages,
  version = vapply(
    required_packages,
    function(pkg) as.character(utils::packageVersion(pkg)),
    FUN.VALUE = character(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  package_versions,
  file.path(OUTPUT_DIR, "package_versions.csv"),
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = file.path(OUTPUT_DIR, "sessionInfo.txt")
)

cat("\n============================================================\n")
cat("Analysis completed successfully.\n")
cat("Selected K:", selected_K, "\n")
cat("Eligible selected-K solutions:", nrow(selected_eligible), "\n")
cat("Final model:", normalizePath(final_model_file, winslash = "/"), "\n")
cat("Results directory:", normalizePath(OUTPUT_DIR, winslash = "/"), "\n")
cat("============================================================\n")
