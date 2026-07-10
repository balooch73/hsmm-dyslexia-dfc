# Data directory

Participant-level data are not distributed with this repository.

Expected local structure:

```text
data/
├── Patients/   # one .xlsx file per child with dyslexia
└── Controls/   # one .xlsx file per TD child
```

Each `.xlsx` file must contain only numeric ROI time series, with time points in rows and identically ordered ROI columns across participants.
