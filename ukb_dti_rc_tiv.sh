#!/bin/bash

#subs=$(<dti_surface_20210802.txt) #get all participant IDs
subs=$(</slow/projects/01_UKB/00_scripts/ukb_rc_ids.txt)
# define function called "task"

task () {
    local run=$1

cd /slow/projects/01_UKB/surface/00_batch1/${j} # cd to results directory

unzip -qq *zip FreeSurfer/stats/aseg.stats
awk 'NR==35 {print $9; exit}' FreeSurfer/stats/aseg.stats>s.txt
sed 's/,//g' s.txt >tiv.txt
rm -rf FreeSurfer/
rm -f s.txt
 
}

N=100
(
for j in $subs; do
   ((i=i%N)); ((i++==0)) && wait
   task "$j" & 
done
)
