import sys
sys.path.append('..')
from typing import List
import numpy as np
from classes.universe import Universe
from classes.helper_functions.helpers import get_spatial_neighbours, get_vertices_in_slice
import argparse
import os


def pick_vertex_dense(universe: Universe) -> int:
    """
    Picks a vertex in the universe at a timeslice with the maximum spatial volume.
    """
    slice_sizes = universe.slice_sizes
    vertices_in_slice = get_vertices_in_slice(universe)
    
    # Find which timeslice has the maximum spatial volume
    max_spatial_volume_timeslice = np.argmax(slice_sizes)
    
    # Pick a random vertex
    vertices = vertices_in_slice[max_spatial_volume_timeslice]
    vertex = np.random.choice(vertices)

    return vertex

def random_walk(universe: Universe, sigma_max: int, n_walkers: int, spatial_flag: int) -> np.ndarray:
    """
    Perform a random walk on the universe for a given number of walkers and sigma_max steps.
    Also counts the number of visits to the starting vertex for each sigma.
    """
    universe.update_vertices()
    
    # Get the (spatial) connections of the universe and the starting vertex
    if spatial_flag == 0:
        connections = get_spatial_neighbours(universe)
    else:
        connections = universe.vertex_neighbours

    v_i = pick_vertex_dense(universe)
   
    # Make a list to count for each sigma how many walkers visited the starting vertes
    n_visits = np.zeros(sigma_max + 1)

    # Start n random walkers from initial vertex for sigma_max steps
    for _ in range(n_walkers):
        v = v_i
        for sigma in range(1, sigma_max + 1):
            # Get the neighbours of the current vertex and move to a random one
            neighbours = connections[v]
            v = np.random.choice(neighbours)

            # If the walker is back at the start, increase the counter for the current sigma
            if v == v_i:
                n_visits[sigma] += 1

    # Compute return probabilities
    return_probabilities = n_visits / n_walkers
    # At sigma = 0, the probability is always 1
    return_probabilities[0] = 1.0
    
    return return_probabilities

def diffusion(universe: Universe, sigma_max: int, n_walkers: int, n_starts: int, spatial_flag: int) -> List[np.ndarray]:
    """
    Perform a diffusion simulation on the universe for a given number of walkers and sigma_max steps.
    Gets the average probability of returning to a starting vertex for each sigma.
    """
    total = []

    # Perform the random walk for multiple starting points and sum the results
    for i in range(n_starts):
        print(f'START ITERATION {i}!')
        return_probabilities = random_walk(universe, sigma_max, n_walkers, spatial_flag)
        total.append(return_probabilities)
    
    return total

def main(T: int, target_volume: int, k0: float, thermal_sweeps: int, main_sweeps: int, chain: int, sigma_max: int, n_walkers: int, n_starts: int, spatial_flag: int) -> None:
    path = f'saved_universes_thermal/k0={k0}/T{T}_k0={k0}_tswps=1000_swps=0_kstps={target_volume * 100}_chain={chain}_thermal_1000.txt'
    universe = Universe(path)

    result = diffusion(universe, sigma_max, n_walkers, n_starts, spatial_flag)

    # Save result
    filename = f'prob_sptl={spatial_flag}_T{T}_N{target_volume}_k0={k0}_tswps={thermal_sweeps}_swps={main_sweeps}_chain={chain}'
    path = f'diffusion/k0={k0}'
    if not os.path.exists(path):
        os.makedirs(path)
    
    np.save(path + '/' + filename, result)
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--T', type=int, default=3)
    parser.add_argument('--target_volume', type=int, default=3000)
    parser.add_argument('--k0', type=float, default=0.0)
    parser.add_argument('--thermal_sweeps', type=int, default=1000)
    parser.add_argument('--main_sweeps', type=int, default=0)
    parser.add_argument('--chain', type=int, default=0)
    parser.add_argument('--sigma_max', type=int, default=200)
    parser.add_argument('--n_walkers', type=int, default=1000)
    parser.add_argument('--n_starts', type=int, default=100)
    parser.add_argument('--spatial_flag', type=int, default=0)
    args = parser.parse_args()

    main(args.T, args.target_volume, args.k0, args.thermal_sweeps, args.main_sweeps, args.chain, args.sigma_max, args.n_walkers, args.n_starts, args.spatial_flag)
