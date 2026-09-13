# EEG-fMRI extinction study

Analysis code for a three-session fear-extinction study with simultaneous
EEG-fMRI. The middle session (Day 2) was acquired with an EEG cap in place;
Days 1 and 3 were not.

## Imaging quality control — [`tsnr_fd_analysis/`](tsnr_fd_analysis/)

**➡️ [Read the results summary](tsnr_fd_analysis/RESULTS.md)** — temporal SNR and
head-motion results with figures, test statistics, and a conclusion on what the
EEG hardware changes.

Headline: grey-matter tSNR is about **6% lower on the EEG day**
(*F*(2, 40) = 6.80, *p* = .003), by the same *relative* amount at every echo, and
head motion is statistically indistinguishable across the three days. The pattern
points to reduced RF receive sensitivity from coil loading rather than to
susceptibility dropout or subject motion.

| File | Contents |
|---|---|
| [`RESULTS.md`](tsnr_fd_analysis/RESULTS.md) | Summary of all tSNR and FD results |
| [`tSNR_results.Rmd`](tsnr_fd_analysis/tSNR_results.Rmd) | Full tSNR analysis (day × echo, and the four EEG runs of Day 2) |
| [`fd_results.Rmd`](tsnr_fd_analysis/fd_results.Rmd) | Full framewise-displacement analysis (mean and max FD) |
| [`qc_plot_helpers.R`](tsnr_fd_analysis/qc_plot_helpers.R) | Shared plotting and significance-bracket code |
| [`extract_fd.py`](tsnr_fd_analysis/extract_fd.py) | Extracts FD summaries from fMRIPrep confounds files |
| [`eeg_fmri_ext_SNR_1_v2.m`](tsnr_fd_analysis/eeg_fmri_ext_SNR_1_v2.m), [`eeg_fmri_ext_SNR_2_v2.sh`](tsnr_fd_analysis/eeg_fmri_ext_SNR_2_v2.sh) | Realignment and tSNR computation |
| [`figures/`](tsnr_fd_analysis/figures/) | Figures used in the summary |

## Skin conductance

Scripts for the skin-conductance analyses of the same study are in the
repository root (`SCR_*.m`, `Plot_*.m`, `*_sc_rest_*.m`, `Mean_tonic_Data.m`,
`TTP_onsets.m`); the analysis plan is in [`instructions.md`](instructions.md).
