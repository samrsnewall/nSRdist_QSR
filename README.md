# nSRdist_QSR

MATLAB code for estimating the probability distribution of normalised sedimentation rates (NSR) from radiocarbon-dated marine sediment cores. This is the code behind

> Newall, S. R., Lisiecki, L. E., Lee, T., & Lawrence, C. E. *Millennial-scale sedimentation rate variability of the deep ocean estimated from radiocarbon-dated sediment cores.* Quaternary Science Reviews (in review).

Code archive: [Zenodo DOI to be added]. Data deposit (Bchron outputs and the results files behind the figures): [Zenodo DOI to be added]. See "Data deposit" below.

## Scientific Background

This code uses radiocarbon dates from sediment cores to estimate **normalised sedimentation rates (NSR)** — the sedimentation rate at each point in a core divided by the average sedimentation rate for that core — and fits probability distributions to them across a large ensemble of cores.

The approach builds on the method of [Lin et al. (2014)](https://doi.org/10.1002/2014PA002713), which used Bchron age-depth models to estimate NSR. The primary development in this codebase is a probabilistic extension that propagates the full calibrated age uncertainty of each radiocarbon date, rather than using only the mode or median of the age model. Several NSR estimation methods are implemented and compared.

The code uses data from Lin et al. (2014), stored in the directory Lin2014Cores/BchronInput, and from the World Atlas created by Mulitza et al. (2022). 

## Pipeline Overview

The code operates in three stages:

**Stage 1 — `calcData.m`**

Reads core metadata and radiocarbon data, filters dates, handles doubly-dated depths and age reversals by constructing scenarios, runs Bchron age-depth modelling (or reads the Bchron outputs provided in `BchronFolders/`), and computes NSR histories using several methods. Results are saved to a `.mat` file in `Results/`.

**Stage 2 — `fitData.m`**

Loads the saved results from Stage 1 and fits four candidate probability distributions to the NSR data: a 2-component Mixture Lognormal (MLN), a Lognormal (LN), a Gamma, and an Inverse Gamma. Fitting is performed under various subset choices and weighting schemes, and results are saved to a new `.mat` file.

**Stage 3 — `ResultsScripts/results_QSR.m`**

Runs Stages 1 and 2 with the manuscript settings (or, with `runS.plotsOnly = true`, loads the saved results), then plots the figures and constructs the tables found in the manuscript. The other scripts in `ResultsScripts/` make the sensitivity tests and goodness-of-fit diagnostics in the response to reviewers (see "Running the Code").

### NSR Estimation Methods

| Method | Description |
|---|---|
| **BMedian** | Uses the median of the Bchron age-depth model at each radiocarbon dated depth. Analogous to, but more sensible than, the Lin et al. (2014) approach.|
| **BSamp** | Samples directly from the Bchron posterior age-depth distribution across many runs. |
| **RSR0 / RSR500 / RSR1000** | Random sampling from calibrated radiocarbon age probability distributions, with an optional minimum age-difference restriction (0, 500, or 1000 years) between consecutive dates. The dates are sampled from the shallowest to the deepest, each constrained to be older than the dates already sampled above it. |
| **RSR<dt>d2s / RSR<dt>random (sensitivity tests)** | The same, sampling the dates from the deepest to the shallowest, or in a random order. Used in the response to reviewers to test the sensitivity of the results to the sampling order. |
| **rejection (sensitivity test)** | Joint rejection sampling: every date is drawn from its full calibrated PDF at once and the draw is kept only if the ages increase with depth, so no sampling order is involved. Not feasible for cores with many overlapping PDFs, which are skipped. |
| **BMode (not used in manuscript)** | Uses the mode of the Bchron age-depth model|
| **LinNSR (not used in manuscript)** | Uses median of the Bchron age-depth model, linear interpolation, and calibrated ages, a replication of the Lin et al. (2014) approach|

---

## Internal Data Formats

### The nSR matrix format

Throughout the codebase, NSR data for a single core (or a single Bchron run) is stored as a **3-row matrix**. This format is used by `nSRBchron.m`, `oneCoreRSR.m`, and the functions that aggregate and fit distributions.

Each matrix has 3 rows and **(K + 1)** columns, where **K** is the number of depth intervals (i.e. the number of accepted dated levels minus one).

```
Row 1: [ NaN,  nSR_1,   nSR_2,   ..., nSR_K   ]
Row 2: [ d_1,  ddep_1,  ddep_2,  ..., ddep_K  ]
Row 3: [ a_1,  dage_1,  dage_2,  ..., dage_K  ]
```

| Row | Column 1 | Remaining columns |
|-----|----------|-------------------|
| 1 | `NaN` (block delimiter) | Normalised sedimentation rates (dimensionless) |
| 2 | First depth of sequence (cm) | Depth differences between consecutive dated levels (cm) |
| 3 | First age of sequence (yr BP) | Age differences between consecutive dated levels (yr) |

The `NaN` in column 1 of row 1 serves as a **block delimiter**. When multiple runs (e.g. BSamp Bchron iterations, or RSR random samples) are concatenated, each run's block is appended horizontally and the NaN marks the start of each new block. Downstream functions use `isnan(matrix(1,:))` to locate block boundaries. Row 2 col 1 holds the absolute first depth and row 3 col 1 holds the absolute first age; these header values are excluded from computations by filtering on `~isnan(matrix(1,:))`.

---

## Requirements

### MATLAB
- MATLAB (tested on MATLAB R2025b; no specific minimum version identified)
- **Statistics and Machine Learning Toolbox** — required for `fitgmdist` (mixture model fitting)
- **Parallel Computing Toolbox** — required for `parfor` loops in `calcData.m` (optional; the code will run without it if `parfor` is replaced with `for`)

### Third-Party MATLAB Function: MatCal
Radiocarbon calibration is performed using **MatCal** (Lougheed & Obrochta, 2016). (A modified version of this is stored within this repository so no need to download separately).
[https://github.com/bryanlougheed/matcal](https://github.com/bryanlougheed/matcal)

### R and Bchron
Bchron age-depth modelling is run via an R script called from MATLAB. You will need:
- **R** installed and the path to 'Rscript' must be provided in calcData.m

The R script 'runBchron.R' will install ('require()') the necessary needed R packages:
- The **Bchron** R package
- The **IntCal** R package

---

## Setup

### 1. Clone or download the repository

```bash
git clone https://github.com/samrsnewall/nSRdist_QSR.git
```

### 2. Run from the repository root

The scripts find the repository root from their own location (`S.sandboxPath` in `calcData.m`) and all other paths (Bchron folders, Lin2014 data, etc.) are relative to it, so run them with the repository root as MATLAB's current folder. `results_QSR.m` and the other scripts in `ResultsScripts/` change to the root themselves.

### 3. Set the Rscript path

For Mac/Linux:
Open Terminal and run 'which Rscript'. Copy and paste this to the `S.RscriptPath` variable in `calcData.m`


For Windows:
Open Command Prompt and run 'where Rscript'. Copy and paste this to the `S.RscriptPath` variable in `calcData.m`

```matlab
S.RscriptPath = "/path/to/Rscript";  % e.g. "/usr/local/bin/Rscript" on Mac/Linux
```

### 4. Ensure MatCal useable

Two options:
1) (Recommended) In `calcData.m` set `S.matcalFast` to true. This uses a modified (faster) copy of matcal that is included within this repo, so no need to install the add on.
2) If `S.matcalFast` set to false, code will need MatCal Add on (https://github.com/bryanlougheed/MatCal). Check whether already on path by entering `which matcal` in Command Window.

### 5. Input data (included)

The code reads core metadata from an Excel spreadsheet (`DataSheets/datasheet.xlsx`) and radiocarbon data from the Mulitza et al. (2022) World Atlas dataset (`WA_Foraminiferal_Isotopes_2022/`) and the Lin et al. (2014) supplementary database (`Lin2014Cores/`). All three are in the repository, at the paths `calcData.m` expects (`S.sheet`, `S.WApath`), so nothing needs setting.

### 6. Bchron outputs (from the data deposit)

`calcData.m` reads pre-computed Bchron outputs from `BchronFolders/<S.BchronFolderName>/` (the manuscript run is `All1_RLGtrue_DS0p05_Jun2`). That folder is 4.5 GB and is not in the repository: download it from the data deposit and place it at `BchronFolders/All1_RLGtrue_DS0p05_Jun2/`. Without it, `calcData.m` runs Bchron itself for every core through R (several hours, and the MCMC draws will differ from the manuscript's).

---

## Running the Code

### The manuscript figures and tables

Open `ResultsScripts/results_QSR.m` and run it. With `runS.plotsOnly = false` it runs `calcData.m` on the production core set (`S.coreSelection = "good"`, 69 cores), then `fitData.m` once for each weighting (depth, age, none), saving `Results/dataT_QSR.mat` and the fit files in `Results/results_QSR/`, and then makes the figures and tables. With `runS.plotsOnly = true` it loads those saved files (for example the ones in the data deposit, see below) and only plots.

The full run takes several hours (the RSRx sampling of 69 cores at 10 000 runs dominates; the three RSRx loops in `calcData.m` are `parfor`).

### Stage 1 and 2 on their own

Open and run `calcData.m`. Key settings at the top of the script control which cores are analysed, which methods are run, and where outputs are saved. The main output is a `.mat` file saved to `Results/`. Then open `fitData.m`, set `filePath` to point to the Stage 1 output file, and run it to fit the distributions and save the results. A driver script can instead set these through the `calcOverride` and `fitOverride` structs, as `results_QSR.m` does.

### The response to reviewers

The remaining scripts in `ResultsScripts/` make the sensitivity tests and diagnostics in the response to reviewers, from the files `results_QSR.m` saves. Each has a `plotsOnly` setting, `true` as committed, in which it loads its saved fit files (in the data deposit) and plots; set it `false` to remake the fit files.

- `samplingDirectionSensitivity_review1_allcores.m` — RSR500 histograms and fits for the three sampling orders (shallow to deep, deep to shallow, random) and the BICeff ranking of each RSRx method under each order. Reads the manuscript fit file directly.
- `samplingDirectionSensitivity_review1_reseeds.m` — reruns `calcData.m` with three new random seeds (about 2.5 h each) and compares the BICeff of the four distributions across seeds, for every RSRx method and order. Fit files in `Results/review1_seeds/`.
- `samplingDirectionSensitivity_review1_jointSampleCores.m` — compares the three sampling orders with the joint (rejection) sampler on the cores where rejection sampling is feasible; needs a second, rejection-only `calcData.m` run with a new seed (the commented-out section at the top shows how). Fit files in `Results/review1_jointSampleCores/`.
- `qq_tailGOF.m`, `residualsKSAD.m`, `quantileSegmentFits.m` — Q-Q plots, residual histograms with Kolmogorov-Smirnov and Anderson-Darling statistics, and quantile-separated fit curves for each distribution and method. Run from the `ResultsScripts/` folder; they read the manuscript fit file.
- `BICeff_nSensitivity.m` — BICeff of the four distributions against the effective number of radiocarbon pairs.

### A Note on Reproducibility

Some processes in Stage 1 (the RSRx age sampling and running Bchron) are stochastic. The Bchron outputs used in the manuscript are provided in the data deposit (`BchronFolders/All1_RLGtrue_DS0p05_Jun2/`) and are read rather than rerun (`S.BchronReDo = false`), so the BMedian and BSamp results are reproduced exactly. The RSRx loops run in parallel, so their random draws are not reproduced by a seed, and rerunning the code gives RSRx results that differ from the manuscript's at the level of Monte Carlo noise (in our reruns the BICeff values moved by less than 3 units and no ranking changed). The conclusions are reproducible; the exact numbers are not. The `calcData.m` and `fitData.m` output files behind the manuscript figures are in the data deposit, so the figures can be remade exactly with `runS.plotsOnly = true`.

---

## Data deposit

The files that are too large for the repository are deposited at [Zenodo DOI to be added]. Unpack them into the repository so that the paths are:

```
BchronFolders/All1_RLGtrue_DS0p05_Jun2/           # Bchron inputs and outputs, 109 cores (4.5 GB)
Results/dataT_QSR.mat                             # calcData output, 69 cores, 10 000 runs (1.0 GB)
Results/dataT_QSR_coreInfo.mat                    # Core metadata read by results_QSR.m
Results/results_QSR/dataT_QSR_fit_*.mat           # fitData outputs: depth, age and none weighting (1.8 GB)
Results/review1_jointSampleCores/*.mat            # Fit files for the rejection-sampler comparison
Results/review1_seeds/dataT_QSR_seed*_fit_*.mat   # Fit files for the reseed comparison
```

With these in place every script runs in its plot-only mode. The reseed `calcData` files (three further 0.7 GB files) and the rejection-only seed-2 `calcData` file are not deposited; the scripts remake them if asked.

---

## Repository Structure

```
nSRdist_QSR/
├── calcData.m                        # Stage 1 main script
├── fitData.m                         # Stage 2 main script
├── Functions/                        # All helper functions, runBchron.R and the calibration curves
├── ResultsScripts/
│   ├── results_QSR.m                 # Stage 3: runs the pipeline and makes the manuscript figures and tables
│   └── ...                           # Sensitivity tests and diagnostics for the response to reviewers
├── BchronFolders/                    # Bchron inputs and outputs (from the data deposit; not in git)
├── BIGMACSdata/                      # The NSR distribution and transition matrix currently used by BIGMACS
├── DataSheets/datasheet.xlsx         # Core metadata and which radiocarbon dates are included or excluded
├── Lin2014Cores/                     # Radiocarbon data of the Lin et al. (2014) cores
├── WA_Foraminiferal_Isotopes_2022/   # Radiocarbon data from the Mulitza et al. (2022) World Atlas
└── Results/                          # .mat output files (from the data deposit or your own runs; not in git
                                      #   except Lin2014_dts.mat, which results_QSR.m reads)
```

---

## Glossary

| Term | Meaning |
|---|---|
| NSR / nSR | Normalised sedimentation rate — sedimentation rate divided by the mean rate for that core |
| SR | Sedimentation rate (cm/kyr) |
| MSPF | Monospecific planktonic foraminifera — the material type used for radiocarbon dating |
| DDD | Doubly-dated depth — a sediment depth with more than one radiocarbon date |
| Scenario | One possible age-depth sequence, constructed to avoid doubly-dated depths and age reversals |
| WA / WA2022 | The Mulitza et al. (2022) World Atlas of foraminiferal isotope data |
| BMode | NSR estimated using the mode of the Bchron age-depth model |
| BMedian | NSR estimated using the median of the Bchron age-depth model |
| BSamp | NSR estimated by sampling from the Bchron posterior |
| RSR | Random Sampling sedimentation Rate — NSR estimated by sampling from calibrated radiocarbon age PDFs |
| MLN | Mixture Lognormal distribution (2-component) |
| LN | Lognormal distribution |
| invGam | Inverse Gamma distribution |
| MSI | Mean Spacing Index — a measure of the resolution of radiocarbon dates in a core |
| BIGMACS | A paleoclimate age model that uses the NSR distribution as a prior |
| CFR | Confirmed (no) reversals — a flag used internally during scenario construction |
| S | The settings structure passed throughout the pipeline |
| fitS | The fitting settings structure used in Stage 2 |

---

## Licence

This code is released under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) licence. You are free to share and adapt the material for any purpose, provided appropriate credit is given.

This repository includes a modified version of **MatCal** (`matcalFast.m`), originally written by Lougheed & Obrochta (2016), which is also licensed under CC BY 4.0. See the `LICENSE` file for full details.

---

## How to cite

Please cite the paper (above) and the code archive: Newall, S. R. (2026). nSRdist_QSR: MATLAB code for estimating normalised sedimentation rate distributions from radiocarbon-dated sediment cores. Zenodo. [Zenodo DOI to be added]

---

## Key References

- Lin, L., Khider, D., Lisiecki, L. E., & Lawrence, C. E. (2014). Probabilistic sequence alignment of stratigraphic records. *Paleoceanography*, 29(6), 976–989. https://doi.org/10.1002/2014PA002713
- Mulitza, S., et al. (2022). World Atlas of late Quaternary foraminiferal oxygen and carbon isotope ratios. *Earth System Science Data*, 14, 2553–2611. https://doi.org/10.5194/essd-14-2553-2022
- Haslett, J., & Parnell, A. (2008). A simple monotone process with application to radiocarbon-dated depth chronologies. *Journal of the Royal Statistical Society: Series C*, 57(4), 399–418. https://doi.org/10.1111/j.1467-9876.2008.00623.x (Bchron)
- Lougheed, B.C. and Obrochta, S.P. (2016). MatCal: Open Source Bayesian 14C Age Calibration in MatLab. *Journal of Open Research Software*, 4, e42. http://dx.doi.org/10.5334/jors.130