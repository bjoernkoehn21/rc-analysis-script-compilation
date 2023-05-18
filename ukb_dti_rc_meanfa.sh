#!/bin/bash

#subs=$(<dti_surface_20210802.txt) #get all participant IDs
subs=$(</slow/projects/01_UKB/00_scripts/ukb_rc_ids.txt)
# define function called "task"

export FREESURFER_HOME=/opt/freesurfer
source $FREESURFER_HOME/SetUpFreeSurfer.sh

task () {
    local run=$1

cd /slow/projects/01_UKB/surface/00_batch1/${j} # cd to results directory

unzip -qq *zip FreeSurfer/mri/aparc+aseg.mgz

mri_binarize --i FreeSurfer/mri/aparc+aseg.mgz --wm --o wm.mask.mgz

/opt/freesurfer/bin/mri_convert -rt nearest -rl /slow/projects/01_UKB/dti/00_batch1/${j}/DWI_processed_v311/*_fractional_anisotropy.nii.gz wm.mask.mgz wm.nii.gz

fslmeants -i /slow/projects/01_UKB/dti/00_batch1/${j}/DWI_processed_v311/*_fractional_anisotropy.nii.gz -m wm.nii.gz -o meanfa.txt

rm -rf FreeSurfer/
rm -f wm*
 
}

N=100
(
for j in $subs; do
   ((i=i%N)); ((i++==0)) && wait
   task "$j" & 
done
)
