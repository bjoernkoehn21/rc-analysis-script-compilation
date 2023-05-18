# Explanation of the initial individual scripts
The current explanations only include the scipts constituting the starting point. These scripts are to be altered, concatenated, generalized, etc. to deliver a comprehensive script compilation fit for public use.  

## Overview
- Compute rc-coefficient for all subjects
  - ukb_compute_rc_dti.m
  - ukb_dti_rc_loop.sh
- Correct and normalize
  - ukb_correct_rc_dti.m
  - ukb_dti_correct_rc_loop.sh
- Calculate density
  - ukb_dti_rc_density.m
- Summarize/assemble relevant meassures
  - ukb_assemble_rc_dti.m
  - ukb_dti_rc_tiv.sh
  - ukb_dti_rc_meanfa.sh

## More detailed information
First of all: The Matlab script, that does the actual RC analysis: So loading the connectivity mats and calculating the RC curves: ukb_compute_rc_dti.m
The script always operates on a person. Therefore you need a wrapper.

Here is my wrapper script. The wrapper is a shell script, which is the easiest to parallelize. It always copies the above Matlab script into the subject folder, executes it, and deletes it again. This way you get the results for each subject saved: ukb_dti_rc_loop.sh

Now we have the RC curves. But: We have to correct the range if necessary. We do this if there is no unique RC regime. For this we have a bash wrapper and a Matlab script: ukb_correct_rc_dti.m and ukb_dti_correct_rc_loop.sh.

And finally: Summarize/assemble the actual measures! For this purpose a Matlab script: ukb_assemble_rc_dti.m

For this purpose it is useful to read out covariates. I have done this somewhat unsystematically so far (i.e. all covariates individually).

Density! So: How many of the theoretically possible connections are there? For this purpose a Matlab script running with parfor in parallel. Here we don't need nested loops, so directly in Matlab: ukb_dti_rc_density.m

Then total intracranial volume. For this, I use the Freesurfer output. The script unpacks the Freesurfer folder of the UKB (the relevant files) and reads the corresponding line. The whole thing in bash: ukb_dti_rc_tiv.sh
For the HCP data it should be sufficient to download the corresponding file from the Freesurfer folder (if not already done).

The same applies to the middle FA. CATO stores an FA image (with the FA value in each voxel), but this is for all tissue classes, and we actually only want white matter. For this I use the white matter mask from FreeSurfer, which I first create and then bring into register with the CATO data. Also a bash script: ukb_dti_rc_meanfa.sh
