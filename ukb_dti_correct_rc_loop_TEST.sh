#!/bin/bash

#subs=$(<dti_surface_20210802.txt) #get all participant IDs
subs=$(</slow/projects/01_UKB/00_scripts/ukb_subjIDs.txt)
# define function called "task"

task () {
    local run=$1

#cd /slow/projects/01_UKB/dti/00_batch1/${j}/DWI_processed_v311 # cd to results directory

# yes | cp -rf /slow/projects/01_UKB/00_scripts/ukb_correct_rc_dti.m ukb_correct_rc_dti.m # copy matlab function to directory

/opt/matlab/bin/matlab -nojvm -nodesktop -r "run /slow/projects/01_UKB/00_scripts/ukb_correct_rc_dti_TEST($1); exit" # invoke matlab and execute script

# rm -f slow/projects/01_UKB/dti/00_batch1/${j}/DWI_processed_v311/ukb_correct_rc_dti.m # no more need for matlab script
 
}

N=100
(
for j in $subs; do
   ((i=i%N)); ((i++==0)) && wait
   task "$j" & 
done
)
