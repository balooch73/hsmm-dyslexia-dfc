# HSMM-based dynamic functional connectivity in developmental dyslexia

This repository contains the analysis code for a Hidden Semi-Markov Model (HSMM) study of task-based dynamic functional connectivity during a reading-aloud task in children with dyslexia and typically developing (TD) children.

## Repository contents

- `R/01_hsmm_model_selection_stability.R`: HSMM model selection for `K = 2–10`, randomized initialisations, final-model estimation, Viterbi decoding, and stability assessment.
- `R/02_load_final_hsmm.R`: optional utility that loads the final model and reconstructs objects for downstream or interactive analyses.
- `R/03_state_characterization.R`: state-specific correlation matrices and mean-signal summaries.
- `R/04_occupancy_sequence_analysis.R`: occupancy estimates, Viterbi-sequence figures, and participant-label permutation tests.
- `R/05_sojourn_duration_analysis.R`: complete internal-run extraction, sojourn summaries, KL-divergence tests, and FDR correction.
- `R/06_transition_probability_analysis.R`: run-level transition probabilities, participant-label permutation tests, and participant-level bootstrap confidence intervals.
- `R/07_topology_threshold_sensitivity.R`: graph-theoretical metrics from complete non-negative weighted networks and visual sensitivity analyses at 15%, 20%, and 25% proportional thresholds.

## Data availability and expected structure

Participant-level imaging data are not included in this repository. To run the full pipeline, place one Excel workbook per participant in the following structure:

```text
data/
├── Patients/   # children with dyslexia
└── Controls/   # typically developing children
```

Each workbook should contain a numeric time-by-ROI matrix: rows are time points, columns are ROIs, and all participants must have the same ROI columns in the same order. The original analysis used 20 ROIs and 200 time points per participant.

Do not upload participant-level data, decoded sequences, fitted `.rds` objects, or files containing identifiable information to the public repository.

## Software requirements

The analysis was written in R. Required packages are:

```r
install.packages(c("readxl", "mhsmm", "mvtnorm", "corpcor", "clue", "igraph"))
```

The scripts save `sessionInfo()` and package-version files with their outputs.

## Running the analysis

1. Clone or download the repository.
2. Set the working directory to the repository root.
3. Place the participant workbooks under `data/Patients/` and `data/Controls/`.
4. Run the scripts in numerical order.

From an R session started in the repository root:

```r
source("R/01_hsmm_model_selection_stability.R")
source("R/03_state_characterization.R")
source("R/04_occupancy_sequence_analysis.R")
source("R/05_sojourn_duration_analysis.R")
source("R/06_transition_probability_analysis.R")
source("R/07_topology_threshold_sensitivity.R")
```

Script 02 is optional and is not required by Scripts 03–07. The first script is computationally intensive because it fits candidate models with multiple randomized initialisations and performs an additional 50-start stability analysis for the selected state number.

By default, outputs are written to:

```text
outputs/HSMM_reviewer_reproducible_analysis/
```

The repository root can also be supplied through the environment variable `HSMM_PROJECT_ROOT`.

## Methodological notes

- HSMM emissions are multivariate Gaussian with state-specific mean vectors and covariance matrices.
- Sojourn distributions are modelled nonparametrically.
- Model selection is based on the minimum BIC from the highest-likelihood eligible solution for each candidate `K`.
- Sojourn analyses exclude each participant's first and last run to avoid boundary-censored durations.
- Transition analyses retain boundary runs because observed transitions between adjacent runs are identifiable.
- Negative correlations are set to zero only for graph analysis.
- Global efficiency, local efficiency, and modularity are calculated from complete non-negative weighted networks. Proportional thresholds of 15%, 20%, and 25% are used only for visualisation.

## Citation

Please cite the associated article when using this code. Full citation details will be added after publication.
