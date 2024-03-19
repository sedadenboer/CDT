# main.py
#
# Author: Seda den Boer
# Date: 29-02-2024
#
# Description: The main file for the 2+1D CDT simulation.

from classes.simulation import Simulation
from classes.universe import Universe
from classes.observable import Observable
import random
import argparse

# Initialize random number generator
rng = random.Random(1)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k0", type=float, help="Value for k0")
    parser.add_argument("--k3", type=float, help="Initial value for k3")
    parser.add_argument("--targetvolume", type=int, help="Value for target volume")
    parser.add_argument("--target2volume", type=int, help="Value for target2 volume")
    parser.add_argument("--volfixswitch", type=int, help="Value for volfix switch")
    parser.add_argument("--seed", type=int, help="Value for seed")
    parser.add_argument("--outputdir", type=str, help="Value for output directory")
    parser.add_argument("--fileid", type=str, help="Value for file ID")
    parser.add_argument("--thermalsweeps", type=int, help="Number of thermalization sweeps")
    parser.add_argument("--measuresweeps", type=int, help="Number of measurement sweeps")
    parser.add_argument("--ksteps", type=int, help="Value for k steps")
    parser.add_argument("--strictness", type=int, help="Value for strictness")
    parser.add_argument("--v1", type=int, help="Value for v1")
    parser.add_argument("--v2", type=int, help="Value for v2")
    parser.add_argument("--v3", type=int, help="Value for v3")
    parser.add_argument("--infile", type=str, help="Value for input file")
    parser.add_argument("--outfile", type=str, help="Value for output file")

    args = parser.parse_args()

    k0 = args.k0 if args.k0 is not None else 0.0
    k3_i = args.k3 if args.k3 is not None else 0.0
    target_volume = args.targetvolume if args.targetvolume is not None else 0
    target2_volume = args.target2volume if args.target2volume is not None else 0
    volfix_switch = args.volfixswitch if args.volfixswitch is not None else 0
    seed = args.seed if args.seed is not None else 0
    output_dir = args.outputdir if args.outputdir is not None else ""
    thermal_sweeps = args.thermalsweeps if args.thermalsweeps is not None else 0
    sweeps = args.measuresweeps if args.measuresweeps is not None else 0
    k_steps = args.ksteps if args.ksteps is not None else 0
    strictness = args.strictness if args.strictness is not None else 0
    v1 = args.v1 if args.v1 is not None else 1
    v2 = args.v2 if args.v2 is not None else 1
    v3 = args.v3 if args.v3 is not None else 1
    in_file = args.infile if args.infile is not None else ""
    out_file = args.outfile if args.outfile is not None else ""

    # Set data directory for observables
    Observable.set_data_dir(output_dir)

    # Initialize Universe
    universe = Universe(geometry_infilename=in_file, strictness=strictness, volfix_switch=volfix_switch)

    print("\n\n#######################")
    print("* * * Initialized * * *")
    print("#######################\n\n")

    #TODO Add observables

    # Start simulation
    simulation = Simulation(universe)
    simulation.start(k0, k3_i, sweeps, thermal_sweeps, k_steps, target_volume, target2_volume, seed, out_file, v1, v2, v3)

    print("\n\n####################")
    print("* * * Finished * * *")
    print("####################\n\n")

    print(f"Number of tetrahedra: {simulation.universe.tetrahedron_pool.get_number_occupied()}\n")
    print(f"Number of (3,1)-tetrahedra: {simulation.universe.tetras_31.get_number_occupied()}\n")
    print(f"Number of vertices: {simulation.universe.vertex_pool.get_number_occupied()}\n")


if __name__ == "__main__":
    main()
