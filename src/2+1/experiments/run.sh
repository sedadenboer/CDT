#!/bin/bash

#SBATCH --partition=gpu_mig
#SBATCH --gpus=1
#SBATCH --job-name=ConvergenceT5
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --time=48:00:00
#SBATCH --output=slurm_output_%A.out
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=seda.denboer@student.uva.nl

module purge
module load 2023
module load Anaconda3/2023.07-2

cd $HOME/CDT/

# Define the values for k0 and k3
k0_values=(3 4 5 6 7 8)
k3_init_guesses=(1.3 1.5 1.7 1.9 2.1 2.3)

# Loop through each combination of k0 and k3 values
for ((i = 0; i < ${#k0_values[@]}; i++)); do
    k0=${k0_values[$i]}
    k3=${k3_init_guesses[$i]}
    echo "Running main.py with k0=$k0 and k3=$k3"
    python main.py\
		--seed 0\
		--k0 $k0\
		--k3 $k3\
		--tune_flag True\
		--thermal_sweeps 1000\
		--sweeps 0\
		--k_steps 500000\
		--target_volume 5000\
		--include_mcmc_data True\
		--measuring_interval 1\
		--measuring_thermal True\
		--save_thermal true\
		--saving_interval 100
done