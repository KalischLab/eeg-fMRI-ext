#!/bin/bash

###############################################################################
#                                                                               
# eeg_fmri_ext_SNR_2.sh                                                     
#                                                                              
# written by B.Meyer, Neuroimaging Center Mainz, 2026                          
#                                                                                
# Temporal SNR check for resting state data acquired at Day 1 (without EEG)     
# and at Day 2 (with EEG). This code calculates tSNR maps for BOLD data from
# Day 1 and 2 and applies GM masks to the results.
#
###############################################################################


# Define variables for subject, session, directories, and runs
input_dir="/home/mey1b/EEG_fmri_ext/data/EEG_fMRI_ext_realigned_for_tSNR"
output_dir="/home/mey1b/EEG_fmri_ext/data/EEG_fMRI_ext_tSNR"
num_echos=4 # Number of echoes

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

results="/home/mey1b/EEG_fmri_ext/data/EEG_fMRI_ext_tSNR/results.txt"

#if [ -f ${results} ]; then
#    rm ${results}
#fi

echo "ID,r1e1,r1e2,r1e3,r1e4,r2e1,r2e2,r2e3,r2e4,\
r3e1,r3e2,r3e3,r3e4,r4e1,r4e2,r4e3,r4e4,\
r5e1,r5e2,r5e3,r5e4,r6e1,r6e2,r6e3,r6e4" >> $results

for subject in $(seq 10 15); do

    for ses_task in "ses-Day1_task-rest1" "ses-Day2_task-rest2" "ses-Day2_task-rest3" "ses-Day2_task-rest4" "ses-Day2_task-rest5" "ses-Day3_task-rest6"; do	   

        # Loop over each echo to calculate tSNR for each echo separately
        for echonr in $(seq 1 $num_echos); do
            echo "Processing echo $echonr for run $run..."

            # Define the echo-specific input file
            bold_file="${input_dir}/sub-S${subject}_${ses_task}_acq-mb3me4_echo-${echonr}_bold_realigned.nii"

	    if [ -f ${bold_file} ]
	    then
		echo "Processing BOLD file ${bold_file}...!"
	    else
		echo "BOLD file ${bold_file} does not exist!"
		
	    	if [[ "$echonr" -eq 1 && "$ses_task" == "ses-Day1_task-rest1" ]]; then
		    echo -n "sub-${subject},NaN," >> $results
		    continue
	    	elif [[ "$echonr" -eq 4 && "$ses_task" == "ses-Day3_task-rest6" ]]; then
		    echo "NaN" >> $results
		    continue
	        else
		    echo -n "NaN," >> $results
		    continue
		fi
	    fi

	    # Define output file names for mean, std, and tSNR for this echo
            mean_file="${output_dir}/sub-S${subject}_${ses_task}_acq-mb3me4_echo-${echonr}_bold_MEAN.nii.gz"
            mean_file_masked="${output_dir}/sub-S${subject}_${ses_task}_acq-mb3me4_echo-${echonr}_bold_MEAN_masked.nii.gz"
	   
	    std_file="${output_dir}/sub-S${subject}_${ses_task}_acq-mb3me4_echo-${echonr}_bold_STD.nii.gz"
            std_file_masked="${output_dir}/sub-S${subject}_${ses_task}_acq-mb3me4_echo-${echonr}_bold_STD_masked.nii.gz"

            tsnr_file="${output_dir}/sub-S${subject}_${ses_task}_acq-mb3me4_echo-${echonr}_bold_tSNR.nii.gz"
            tsnr_file_masked="${output_dir}/sub-S${subject}_${ses_task}_acq-mb3me4_echo-${echonr}_bold_tSNR_masked.nii.gz"

            gm_mask="${input_dir}/c1sub-S${subject}_ses-${ses_task:14}_T1w.nii"
            gm_mask_resl="${input_dir}/c1sub-S${subject}_ses-${ses_task:14}_T1w_resl.nii"
            gm_mask_resl_bin="${input_dir}/c1sub-S${subject}_ses-${ses_task:14}_T1w_resl_bin.nii"

            in2ref_mat="${input_dir}/in2ref.mat"
            in2ref_im="${input_dir}/in2ref.nii.gz"

             # Calculate the mean BOLD image across time for this echo
            fslmaths "$bold_file" -Tmean "$mean_file"

	    # Calculate the standard deviation image across time for this echo
            fslmaths "$bold_file" -Tstd "$std_file"
    
	    # Reslice GM mask with respect to BOLD image
            flirt -in "$gm_mask" -ref "$mean_file" -omat "$in2ref_mat" -out "$in2ref_im"
            flirt -in "$gm_mask" -ref "$mean_file" -applyxfm -init "$in2ref_mat" -out "$gm_mask_resl"           

            # Binarize GM mask thresholded at 0.5
            fslmaths $gm_mask_resl -thr 0.5 -bin $gm_mask_resl_bin

            # Mask mean BOLD image with GM mask
            fslmaths "$mean_file" -mas "$gm_mask_resl_bin" "$mean_file_masked"

            # Mask STD BOLD image with GM mask
            fslmaths "$std_file" -mas "$gm_mask_resl_bin" "$std_file_masked"

            # Calculate tSNR by dividing the mean image by the standard deviation
            fslmaths "$mean_file" -div "$std_file" "$tsnr_file"

            # Mask TSNR BOLD image with GM mask
            fslmaths "$mean_file_masked" -div "$std_file_masked" "$tsnr_file_masked"
	    v=$(fslstats "$tsnr_file_masked" -k "$tsnr_file_masked" -M)
	    v="${v:0:-1}"

	    if [[ "$echonr" -eq 1 && "$ses_task"  == "ses-Day1_task-rest1" ]]; then
		echo -n "sub-${subject},${v}," >> $results
	    elif [[ "$echonr" -eq 4 && "$ses_task" == "ses-Day3_task-rest6" ]]; then
		echo "$v" >> $results
	    else
	        echo -n "$v," >> $results
	    fi

            echo "tSNR calculation completed for $run, echo $echonr."
        done
    	echo "Processing BOLD file ${bold_file} completed!"
    done
    
done


