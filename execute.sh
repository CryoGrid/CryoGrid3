#!/bin/bash

#change to model directory
cd $HOME/CryoGrid/CryoGrid3
printf "Switched to CryoGrid3 directory...\n"

# run setup script which generates setup files
matlab -nosplash -nodisplay -nodesktop < setup_CryoGrid3.m > MATLAB.temp
#rm MATLAB.temp
printf "Configured run and created output directory...\n"

# generate run and submit scripts from templates
sed -e "s@SAVEDIR@$(< SAVEDIR.temp)@g" -e "s/RUN/$(< RUN.temp)/g" submit_SLURM.template > submit_SLURM.sh
sed -e "s@SAVEDIR@$(< SAVEDIR.temp)@g" -e "s/RUN/$(< RUN.temp)/g" run_CryoGrid3.template > run_$(< RUN.temp).m
#rm SAVEDIR.temp
#rm RUN.temp
printf "Created run and submit scripts from templates...\n"

# submit the batch job to the slurm queue
sbatch submit_SLURM.sh
printf "Submitted job to SLURM queue!\n"
