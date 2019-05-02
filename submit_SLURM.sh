#!/bin/bash

#SBATCH --job-name="SCENARIO_REV4_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_yedomaREV2_DarcyReservoirNoInflow_eRes-10.00_snowDens250_DeltaXice0.00_diff0.0_adv1.0_Kland3.0e-10_Kwater3.0e-08"
#SBATCH --qos=medium
#SBATCH --account=permarisk
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --output=/data/scratch/nitzbon/CryoGrid/CryoGrid3_infiltration_xice_mpi_polygon/runs/SCENARIO_REV4_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_yedomaREV2_DarcyReservoirNoInflow_eRes-10.00_snowDens250_DeltaXice0.00_diff0.0_adv1.0_Kland3.0e-10_Kwater3.0e-08/CryoGrid3.out
#SBATCH --error=/data/scratch/nitzbon/CryoGrid/CryoGrid3_infiltration_xice_mpi_polygon/runs/SCENARIO_REV4_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_yedomaREV2_DarcyReservoirNoInflow_eRes-10.00_snowDens250_DeltaXice0.00_diff0.0_adv1.0_Kland3.0e-10_Kwater3.0e-08/CryoGrid3.err
#SBATCH --workdir=/home/nitzbon/CryoGrid/CryoGrid3

matlab -nosplash -nodisplay -nodesktop < run_SCENARIO_REV4_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_yedomaREV2_DarcyReservoirNoInflow_eRes-10.00_snowDens250_DeltaXice0.00_diff0.0_adv1.0_Kland3.0e-10_Kwater3.0e-08.m
rm run_SCENARIO_REV4_199910-209912_rcp85_xice1_xE1_xH1_xW1_xS1_yedomaREV2_DarcyReservoirNoInflow_eRes-10.00_snowDens250_DeltaXice0.00_diff0.0_adv1.0_Kland3.0e-10_Kwater3.0e-08.m