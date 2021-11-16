#!/bin/bash

# This script is to set up SEQ data on talapas in BIDS (ish) format
# Run 01_seq_setup.sh first, then create TSV files on local computer and transfer to 'bids_data/$SSID/func/' (or figure out a  way to create json files instead)

#SBATCH --account=mayrlab --time=0-2:00:00 --output=logs/01_setup_%j.txt


BIDS_DIR=/gpfs/projects/mayrlab/shared/SEQ/bids_data/
OUTPUT_DIR=/gpfs/projects/mayrlab/shared/SEQ/bids_data_preprocessed/
WORKING_DIR=/gpfs/projects/mayrlab/shared/SEQ/bids_data_tmp/
FS_FILE=/gpfs/projects/mayrlab/shared/SEQ_prep/fMRI_prep/fs_license_mmoss.txt

fmriprep -v --participant-label $SUB -t $TASK_NAME --skull-strip-t1w force –fs-license-file $FS_FILE  -w $WORKING_DIR $BIDS_DIR $OUTPUT_DIR $SUB



fmriprep --participant-label 301 -t Run1 NAME --skull-strip-t1w force –fs-license-file /gpfs/projects/mayrlab/shared/SEQ_prep/fMRI_prep/fs_license_mmoss.txt  -w /gpfs/projects/mayrlab/shared/SEQ/bids_data_tmp/ /gpfs/projects/mayrlab/shared/SEQ/bids_data/ /gpfs/projects/mayrlab/shared/SEQ/bids_data_preprocessed/ 301

