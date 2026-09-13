# EEG-fMRI extinction study — imaging quality control

Quality-control analyses for a three-session multi-echo resting-state dataset in
which the middle session (Day 2) was acquired with an EEG cap in place.

**➡️ [Read the results summary](RESULTS.md)** — tSNR and head-motion results with
figures, test statistics, and a conclusion on what the EEG hardware changes.

Headline: grey-matter tSNR is about **6% lower on the EEG day**
(*F*(2, 40) = 6.80, *p* = .003), by the same *relative* amount at every echo, and
head motion is statistically indistinguishable across the three days. The pattern
points to reduced RF receive sensitivity from coil loading rather than to
susceptibility dropout or subject motion.

## Contents

- [`RESULTS.md`](RESULTS.md) — summary of all results
- [`tSNR_results.Rmd`](tSNR_results.Rmd) — full tSNR analysis (day × echo, and the four EEG runs of Day 2)
- [`fd_results.Rmd`](fd_results.Rmd) — full framewise-displacement analysis (mean and max FD)
- [`qc_plot_helpers.R`](qc_plot_helpers.R) — shared plotting and significance-bracket code
- [`extract_fd.py`](extract_fd.py) — extracts FD summaries from fMRIPrep confounds files
- `eeg_fmri_ext_SNR_1_v2.m`, `eeg_fmri_ext_SNR_2_v2.sh` — realignment and tSNR computation
- `figures/` — figures used in the summary

Skin-conductance scripts for the same study are also in this repository
(`SCR_*.m`, `Plot_*.m`, `*_sc_rest_*.m`); see [`instructions.md`](instructions.md).
