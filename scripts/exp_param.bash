#!/bin/bash
#############################################################
# <next few lines under version control, D O  N O T  E D I T>
# $Date$
# $Author$
# $Revision$
# $Id$
#############################################################

#############  U S E R  I N P U T  -  B E G I N  ############

# define experiment names and run types
exp_name="shear1.0_amu1.0"   # Name of the experiment
run_type="control"           # Type of run [ truth / control / perturbed ]
let run_number=1             # Number of the run [ 1,2,3, ... ]
let Ne=5000                  # Ensemble size

# define "local source" and "remote storage" experiment directories 
dir_expt_src="/home/disk/pvort/rahulm/sqg_experiments"
dir_expt_run="/home/disk/p/nobackup/rahulm/ALL_ETC/sqg_experiments"

# Based on the choice of ${run_type} above, set below:
if [ ${run_type} == "preptruth" ]; then
	
	echo "THIS IS A PREPTRUTH RUN"

elif [ ${run_type} == "truth" ]; then
	
	echo "THIS IS A TRUTH RUN"

elif [ ${run_type} == "control" ]; then

	echo "THIS IS A CONTROL RUN"
	
	let tint=21                # no. of times in each ensemble member
	let truth_run_number=1     # Run Number of TRUTH

elif [ ${run_type} == "perturbed" ]; then
	
	echo "THIS IS A PERTURBED RUN"
	
	let control_run_number=1   # Run Number of CONTROL
	exp_type="THETA"           # Type of experiment [ THETA / PHI ] 

fi

###############  U S E R  I N P U T  -  E N D  ##############
