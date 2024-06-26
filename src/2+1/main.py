# main.py
# 
# Author: Seda den Boer
# Date: 15-05-2024
#
# Description: Main file for running the 2+1D CDT Monte Carlo simulation.

import sys
sys.path.append('..')
import argparse
from classes.simulation import Simulation
from classes.universe import Universe
import numpy as np


def main(args):
    # Create an instance of the Universe class
    if args.intial_geometry:
        # Create an initial universe with the provided number of time slices
        universe = Universe(geometry_infilename=f'classes/initial_universes/initial_g=0_T={args.T}.CDT')
        k3 = args.k3
    else:
        # Load the universe and k3 value from the provided files (has to be adjusted dependent on the path and file names)
        geometry_infile = 'saved_universes_thermal/k0={args.k0}/' \
                         f'T{args.T}_k0={args.k0}_tswps=1000_swps=0_kstps={args.k_steps}_chain={args.seed}_thermal_1000.txt'
        k3_infile = 'measurements_thermal/k0={args.k0}/' \
                         f'T{args.T}_k0={args.k0}_tswps=1000_swps=0_kstps={args.k_steps}_chain={args.seed}_k3_values.npy'
        
        # Load the universe and k3 value
        universe = Universe(geometry_infilename=geometry_infile)
        k3_values = np.load(k3_infile)
        k3 = k3_values[-1]

        print(f'Loaded universe with k0={args.k0} and k3={k3}')

    # Create an instance of the Simulation class with the provided arguments
    simulation = Simulation(
        universe=universe,
        seed=args.seed,
        k0=args.k0,
        k3=k3,
        tune_flag=args.tune_flag,
        thermal_sweeps=args.thermal_sweeps,
        sweeps=args.sweeps,
        k_steps=args.k_steps,
        v1=args.v1,
        v2=args.v2,
        v3=args.v3,
        volfix_switch=args.volfix_switch,
        target_volume=args.target_volume,
        target2_volume=args.target2_volume,
        epsilon=args.epsilon,
        observables=args.observables,
        include_mcmc_data=args.include_mcmc_data,
        measuring_interval=args.measuring_interval,
        measuring_thermal=args.measuring_thermal,
        measuring_main=args.measuring_main,
        save_main=args.save_main,
        save_thermal=args.save_thermal,
        saving_interval=args.saving_interval,
        validity_check=args.validity_check,
    )

    # Start the simulation
    simulation.start(
        outfile=f'T{args.T}_k0={simulation.k0}_tswps={args.thermal_sweeps}_swps={args.sweeps}_kstps={args.k_steps}_chain={args.seed}'
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Monte Carlo simulation")
    parser.add_argument("--intial_geometry", action="store_true", help="Flag to use initial universe")
    parser.add_argument("--T", type=int, default=3, help="Number of time slices of the universe")
    parser.add_argument("--seed", type=int, default=42, help="Seed for random number generator")
    parser.add_argument("--k0", type=float, default=1, help="Number of k0 moves to perform")
    parser.add_argument("--k3", type=float, default=1, help="Number of k3 moves to perform")
    parser.add_argument("--tune_flag", action="store_true", help="Flag to tune the k3 parameter")
    parser.add_argument("--thermal_sweeps", type=int, default=10, help="Number of thermal sweeps to perform")
    parser.add_argument("--sweeps", type=int, default=10, help="Number of sweeps to perform")
    parser.add_argument("--k_steps", type=int, default=1000, help="Number of k steps to perform")
    parser.add_argument("--v1", type=int, default=1, help="Frequency of the v1 move")
    parser.add_argument("--v2", type=int, default=1, help="Frequency of the v2 move")
    parser.add_argument("--v3", type=int, default=1, help="Frequency of the v3 move")
    parser.add_argument("--volfix_switch", type=int, default=0, help="Volfix switch")
    parser.add_argument("--target_volume", type=int, default=0, help="Target volume")
    parser.add_argument("--target2_volume", type=int, default=0, help="Second target volume")
    parser.add_argument("--epsilon", type=float, default=0.00005, help="Epsilon parameter")
    parser.add_argument("--observables", nargs="+", default=[
        'n_vertices', 'n_tetras', 'n_tetras_31', 'n_tetras_22', 'slice_sizes', 'slab_sizes', 'curvature', 'connections'
        ], help="List of observables to measure")
    parser.add_argument("--include_mcmc_data", action="store_true", help="Flag to include MCMC data in the observables)")
    parser.add_argument("--measuring_interval", type=int, default=1, help="Measuring interval")
    parser.add_argument("--measuring_thermal", action="store_true",help="Flag to measure thermal data")
    parser.add_argument("--measuring_main", action="store_true", help="Flag to measure main data")
    parser.add_argument("--save_main", action="store_true", help="Flag to save main data")
    parser.add_argument("--save_thermal", action="store_true",help="Flag to save thermal data")
    parser.add_argument("--saving_interval", type=int, default=1, help="Saving interval")
    parser.add_argument("--validity_check", action="store_true", help="Flag to perform validity check")
  
    args = parser.parse_args()
    main(args)
