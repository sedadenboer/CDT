import sys
sys.path.append('..')
import argparse
from classes.simulation import Simulation
from classes.universe import Universe
import numpy as np


def main(args):
    # Set up the universe
    universe = Universe(total_time=args.T, initial_slice_size=args.slice_size)

    # Set up the simulation
    simulation = Simulation(
        universe=universe,
        lambd=args.lambd,
        seed=args.seed,
        steps=args.steps,
        weighted_moves=args.weighted_moves,      
    )

    # Run the simulation
    simulation.progress_universe(silence=args.silence, save_data=args.save)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Monte Carlo simulation")
    parser.add_argument('--T', type=int, help='Total time', default=20)
    parser.add_argument('--slice_size', type=int, help='Initial slice size', default=20)
    parser.add_argument('--lambd', type=float, help='Lambda', default=np.log(2))
    parser.add_argument('--seed', type=int, help='Seed', default=42)
    parser.add_argument('--steps', type=int, help='Steps', default=100)
    parser.add_argument('--weighted_moves', action='store_true', help='Weighted moves')
    parser.add_argument('--silence', action='store_true', help='Silence')
    parser.add_argument('--save', action='store_true', help='Save data')
    args = parser.parse_args()
    main(args)