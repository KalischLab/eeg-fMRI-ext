#!/bin/bash
###############################################################################
#
# eeg_fmri_ext_SNR_2_v2.sh
#
# written by B.Meyer, Neuroimaging Center Mainz, 2026
# corrected version of eeg_fmri_ext_SNR_2.sh
#
# Temporal SNR check for the resting-state runs of all three acquisition days
# (Day 1: rest1, no EEG;  Day 2: rest2-5, with EEG;  Day 3: rest6, no EEG).
# Computes a GM-masked tSNR map per run and echo and writes the grey-matter
# mean to a wide CSV, plus a long-format QC file.
#
# Input comes from eeg_fmri_ext_SNR_1_v2.m: motion-corrected 4D BOLD per echo
# and a GM tissue map already coregistered to that run's mean EPI.
#
# Changes vs. v1:
#   1. Subjects 10-34, matching script 1. v1 looped seq 10 15.
#   2. Rows are assembled in an array and printed once, so every row is
#      guaranteed to have 1 + n_runs*n_echoes fields. v1 emitted commas and
#      newlines from inside the loop based on which cell it was on, which
#      breaks silently if the script is interrupted.
#   3. Results are written to a temp file and moved into place, and the header
#      is written exactly once. v1 appended with >> and had its rm commented
#      out, so rerunning duplicated the header and all rows.
#   4. The time series is detrended before the temporal SD is computed. tSNR
#      on an undetrended series is a drift-plus-noise ratio, and drift depends
#      on scanner state, which is exactly what differs between days.
#   5. The GM mask is resampled using the coregistration script 1 already
#      established (-applyxfm -usesqform). v1 re-estimated a fresh 12-DOF
#      registration between a tissue probability map and a mean EPI, which is
#      poorly posed (inverted GM/WM contrast) and can scale or shear the mask.
#   6. The mask is explicit: GM > 0.5, intersected with voxels that have a
#      non-zero temporal SD. v1 used the tSNR map as its own mask, so the
#      voxel count varied between cells for uncontrolled reasons.
#   7. Per-cell intermediate filenames, so nothing is clobbered and a parallel
#      run is safe. v1 reused one fixed in2ref.mat/in2ref.nii.gz everywhere.
#   8. tSNR is computed once, on the masked images only. v1 also computed an
#      unmasked map that divided by zero outside the brain and was never used.
#   9. Voxel count, median and SD are written to a QC file, so a misregistered
#      or shrunken mask is visible instead of silently shifting the mean.
#  10. Log lines report the run being processed. v1 printed "$run", unset.
#
# Non-steady-state volumes are handled in script 1 (n_dummy), before
# realignment, which is where they have to be dropped.
#
# Requires: FSL (fslmaths, fslstats, fslval, flirt), bc.
#
###############################################################################

set -uo pipefail
export FSLOUTPUTTYPE=NIFTI_GZ

#%% ---------------------------- configuration -----------------------------

input_dir="/home/mey1b/EEG_fmri_ext/data/EEG_fMRI_ext_realigned_for_tSNR"
output_dir="/home/mey1b/EEG_fmri_ext/data/EEG_fMRI_ext_tSNR"

subjects=$(seq 10 34)
runs=(1 2 3 4 5 6)
days=(1 2 2 2 2 3)          # runs[i] was acquired on days[i]
num_echoes=4

gm_threshold=0.5            # GM probability threshold for the mask
hp_cutoff_s=128             # high-pass cutoff in seconds for detrending
keep_intermediates=0        # 1 to keep MEAN/STD/tSNR maps and masks

results="${output_dir}/tSNR_results.txt"
qc="${output_dir}/tSNR_qc.csv"

#%% -------------------------------- setup ---------------------------------

mkdir -p "$output_dir"

for cmd in fslmaths fslstats fslval flirt bc; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: $cmd not found on PATH" >&2
        exit 1
    fi
done

results_tmp=$(mktemp "${output_dir}/.tSNR_results.XXXXXX")
qc_tmp=$(mktemp "${output_dir}/.tSNR_qc.XXXXXX")
trap 'rm -f "$results_tmp" "$qc_tmp"' EXIT

# Header, built from the loop variables so it cannot drift out of sync
# with the columns actually written.
hdr="ID"
for r in "${runs[@]}"; do
    for e in $(seq 1 "$num_echoes"); do
        hdr+=",r${r}e${e}"
    done
done
printf '%s\n' "$hdr" > "$results_tmp"
printf 'ID,run,day,echo,mean_tSNR,median_tSNR,sd_tSNR,n_voxels\n' > "$qc_tmp"

#%% ------------------------------ functions -------------------------------

# process_cell <subject> <run> <day> <echo>
#
# Prints "mean median sd nvox" on stdout and returns 0 on success.
# Returns non-zero if the cell cannot be computed. All logging goes to
# stderr so that it never contaminates the captured value.
process_cell() {
    local sub="$1" run="$2" day="$3" echonr="$4"

    local tag="sub-S${sub}_ses-Day${day}_task-rest${run}"
    local cell="${tag}_acq-mb3me4_echo-${echonr}"

    local bold="${input_dir}/${cell}_bold_realigned.nii"
    [[ -f "$bold" ]] || bold="${bold}.gz"
    local gm="${input_dir}/c1${tag}_T1w.nii"

    if [[ ! -f "$bold" ]]; then
        echo "  MISSING BOLD: ${cell}_bold_realigned.nii" >&2
        return 1
    fi
    if [[ ! -f "$gm" ]]; then
        echo "  MISSING GM MASK: c1${tag}_T1w.nii" >&2
        return 1
    fi

    # Per-cell intermediates: nothing is shared between cells.
    local pfx="${output_dir}/${cell}"
    local mean_file="${pfx}_MEAN.nii.gz"
    local detr_file="${pfx}_DETREND.nii.gz"
    local std_file="${pfx}_STD.nii.gz"
    local tsnr_file="${pfx}_tSNR.nii.gz"
    local gm_resl="${pfx}_GM_resl.nii.gz"
    local gm_bin="${pfx}_GM_bin.nii.gz"
    local nz_file="${pfx}_NONZERO.nii.gz"
    local mask_file="${pfx}_GMmask.nii.gz"

    # --- numerator: temporal mean of the uncorrected series ---------------
    fslmaths "$bold" -Tmean "$mean_file" || return 1

    # --- denominator: temporal SD AFTER removing slow drift ---------------
    # Drift is not thermal noise. Leaving it in makes tSNR depend on scanner
    # state, which is confounded with acquisition day.
    local tr sigma
    tr=$(fslval "$bold" pixdim4 2>/dev/null | tr -d '[:space:]')

    # Validate as a positive decimal before handing it to bc: an unset
    # pixdim4 can come back as "0", "" or "nan" depending on how the file
    # was written, and bc's handling of "nan" is not portable.
    if [[ ! "$tr" =~ ^[0-9]*\.?[0-9]+$ ]] || (( $(echo "$tr <= 0" | bc -l) )); then
        echo "  WARNING: no valid TR in header of ${cell}; NOT detrending" >&2
        fslmaths "$bold" -Tstd "$std_file" || return 1
    else
        # fslmaths -bptf takes sigma in volumes: (cutoff_s / 2) / TR
        sigma=$(echo "scale=6; ($hp_cutoff_s / 2) / $tr" | bc -l)
        fslmaths "$bold" -bptf "$sigma" -1 "$detr_file" || return 1
        # -bptf demeans; the SD is unaffected by that, and the numerator
        # comes from the unfiltered mean above.
        fslmaths "$detr_file" -Tstd "$std_file" || return 1
    fi

    # --- GM mask in this run's EPI space ---------------------------------
    # Script 1 already coregistered the tissue map to this run's mean EPI,
    # so this is a pure resampling that respects the headers. Do NOT let
    # flirt estimate a new transform here.
    flirt -in "$gm" -ref "$mean_file" -applyxfm -usesqform -out "$gm_resl" || return 1

    fslmaths "$gm_resl" -thr "$gm_threshold" -bin "$gm_bin" || return 1

    # Voxels with zero temporal SD carry no signal (outside the FOV for part
    # of the series, or dead). Excluding them explicitly means the mask, and
    # therefore the denominator of the average, is well defined.
    fslmaths "$std_file" -bin "$nz_file" || return 1
    fslmaths "$gm_bin" -mas "$nz_file" "$mask_file" || return 1

    # --- tSNR -------------------------------------------------------------
    fslmaths "$mean_file" -div "$std_file" -mas "$mask_file" "$tsnr_file" || return 1

    local stats m sd med nvox vol
    stats=$(fslstats "$tsnr_file" -k "$mask_file" -m -s -p 50 -v) || return 1
    read -r m sd med nvox vol <<< "$stats"

    if [[ -z "${m:-}" || -z "${nvox:-}" ]] || (( $(echo "$nvox < 1" | bc -l) )); then
        echo "  EMPTY MASK for ${cell}" >&2
        return 1
    fi

    if (( keep_intermediates == 0 )); then
        rm -f "$mean_file" "$detr_file" "$std_file" "$tsnr_file" \
              "$gm_resl" "$gm_bin" "$nz_file" "$mask_file"
    fi

    printf '%s %s %s %s\n' "$m" "$med" "$sd" "$nvox"
    return 0
}

#%% ------------------------------ main loop -------------------------------

for sub in $subjects; do

    row="sub-${sub}"
    echo "=== sub-${sub} ===" >&2

    for i in "${!runs[@]}"; do
        run="${runs[$i]}"
        day="${days[$i]}"

        for echonr in $(seq 1 "$num_echoes"); do

            echo "Processing sub-${sub}, Day ${day}, rest${run}, echo ${echonr}" >&2

            if out=$(process_cell "$sub" "$run" "$day" "$echonr"); then
                read -r m med sd nvox <<< "$out"
                row+=",${m}"
                printf '%s,%s,%s,%s,%s,%s,%s,%s\n' \
                    "sub-${sub}" "$run" "$day" "$echonr" "$m" "$med" "$sd" "$nvox" >> "$qc_tmp"
            else
                row+=",NaN"
                printf '%s,%s,%s,%s,NaN,NaN,NaN,NaN\n' \
                    "sub-${sub}" "$run" "$day" "$echonr" >> "$qc_tmp"
            fi
        done
    done

    printf '%s\n' "$row" >> "$results_tmp"
done

#%% ------------------------------- finish ---------------------------------

mv "$results_tmp" "$results"
mv "$qc_tmp" "$qc"
trap - EXIT

echo >&2
echo "wrote $results" >&2
echo "wrote $qc" >&2
