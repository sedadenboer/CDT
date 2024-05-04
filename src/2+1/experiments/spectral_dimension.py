# Load packages
import sys
sys.path.append('..')
from typing import Tuple
import numpy as np
from classes.universe import Universe
from classes.helper_functions.helpers import get_spatial_neighbours, get_vertices_in_slice
import argparse
import os


def pick_vertex_dense(universe: Universe):
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

def measure_spectral_dimension(
        universe: Universe,
        sigma_max: int,
        diffusion_constant: float,
        only_spatial: bool,
        start_flag: int,
        T: int,
        k0: float,
        target_volume: int,
        thermal_sweeps: int,
        main_sweeps: int,
        chain: int,
        run: int
        ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Measure the spectral dimension of the universe.

    Args:
        universe (Universe): Universe object
        sigma_max (int): Maximum sigma value
        diffusion_constant (float): Diffusion constant
        only_spatial (bool): If True, only consider spatial neighbors
        start_flag (int): 0 if start vertex should be first vertex, 1 if start vertex should be picked randomly, \
            2 if start vertex should be picked randomly from the densest timeslice
        T (int): Number of timeslices
        k0 (float): Coupling constant
        target_volume (int): Target volume
        thermal_sweeps (int): Number of thermalization sweeps
        main_sweeps (int): Number of main sweeps
        chain (int): Chain number
        run (int): Run number

    Returns:
        prob (np.ndarray): Array of probabilities
        spec_dim (np.ndarray): Array of spectral dimensions
        spec_dim1 (np.ndarray): Array of spectral dimensions
    """
    # Initialize variables
    if only_spatial:
        connectivity = get_spatial_neighbours(universe)
    else:
        connectivity = universe.vertex_neighbours

    if start_flag == 0:
        start = connectivity.keys()[0]
    elif start_flag == 1:
        start = np.random.choice(list(connectivity.keys()))
    else:
        start = pick_vertex_dense(universe)

    spec_dim = np.zeros(sigma_max)
    spec_dim1 = np.zeros(sigma_max)
    prob_buffers = [np.zeros(len(connectivity)) for _ in range(2)] 
    cur = 0
    prob_buffers[cur][start] = 1.0 
    epsilon = 1e-10  # small value to avoid numerical instability
    prob = np.zeros(sigma_max)

    # Main loop over sigma values
    for sigma in range(sigma_max):
        prob[sigma] = prob_buffers[cur][start]  # Store probability for current sigma

        # Update probability buffers for each vertex
        for vertex_id, neighbors in connectivity.items():
            # Make sure that the probability is not too small
            if prob_buffers[cur][vertex_id] > epsilon:
                # Update probability for each neighbor
                for neighbor_id in neighbors:
                    # If this vertex could not have been reached by diffusion in half of the time,
                    # there is no probability to return to the starting vertex
                    if sigma > sigma_max / 2 and prob_buffers[cur][neighbor_id] < epsilon:
                        continue
                    
                    # Update probability
                    prob_buffers[(cur + 1) % 2][neighbor_id] += (diffusion_constant *
                                                                prob_buffers[cur][vertex_id] / len(neighbors))
            
            # Prevent oscillations at low sigma (diffusion particles can sometimes remain where they are)
            prob_buffers[(cur + 1) % 2][vertex_id] += ((1.0 - diffusion_constant) *
                                                        prob_buffers[cur][vertex_id])

        # Reset current buffer
        for i in range(len(connectivity)):
            prob_buffers[cur][i] = 0.0

        cur = (cur + 1) % 2

    # Update spectral dimension
    for sigma in range(2, sigma_max - 1):
        spec_dim[sigma] = -2.0 * sigma * (prob[sigma + 1] / prob[sigma] - 1.0)
        spec_dim1[sigma] = -2.0 * np.log(prob[sigma + 1] / prob[sigma]) / np.log1p(1.0 / sigma)

    # Check if directory exists
    if not os.path.exists(f'thermal_{target_volume}/T{T}/prob'):
        os.makedirs(f'thermal_{target_volume}/T{T}/prob')
    if not os.path.exists(f'thermal_{target_volume}/T{T}/specdim'):
        os.makedirs(f'thermal_{target_volume}/T{T}/specdim')
    if not os.path.exists(f'thermal_{target_volume}/T{T}/specdim1'):
        os.makedirs(f'thermal_{target_volume}/T{T}/specdim1')

    if only_spatial:
        sptl = 0
    else:
        sptl = 1

    output_prob = f'thermal_{target_volume}/T{T}/prob/prob_sptl={sptl}_strt={start_flag}_T{T}_N{target_volume}_k0={k0}_tswps={thermal_sweeps}_swps={main_sweeps}_chain={chain}_run_{run}.txt'
    output_spec_dim = f'thermal_{target_volume}/T{T}/specdim/specdim_sptl={sptl}_strt={start_flag}_T{T}_N{target_volume}_k0={k0}_tswps={thermal_sweeps}_swps={main_sweeps}_chain={chain}_run_{run}.txt'
    output_spec_dim1 = f'thermal_{target_volume}/T{T}/specdim1/specdim1_sptl={sptl}_strt={start_flag}_T{T}_N{target_volume}_k0={k0}_tswps={thermal_sweeps}_swps={main_sweeps}_chain={chain}_run_{run}.txt'

    # Save results
    np.savetxt(output_prob, prob)
    np.savetxt(output_spec_dim, spec_dim)
    np.savetxt(output_spec_dim1, spec_dim1)

    return prob, spec_dim, spec_dim1


def main(T: int,
         k0: float,
         target_volume: int,
         thermal_sweeps: int,
         main_sweeps: int,
         chain: int,
         sigma_max: int,
         diffusion_constant:
         bool,
         start_flag: int,
         runs: int
         ):
    # Load universe
    # filename = f'../experiments/thermal_{target_volume}/T{T}/saved_universes/k0={k0}/T{T}_k0={k0}_tswps=1000_swps=0_kstps={target_volume * 100}_chain={int(k0)}_thermal_1000.txt'
    filename = f'saved_universes/k0={k0}/T{T}_k0={k0}_tswps={thermal_sweeps}_swps={main_sweeps}_kstps={target_volume * 100}_chain={chain}_thermal_1000.txt'
    universe = Universe(geometry_infilename=filename)
    universe.update_vertices()
    
    # Measure spectral dimension
    for run in range(runs):
        measure_spectral_dimension(
            universe=universe,
            sigma_max=sigma_max,
            diffusion_constant=diffusion_constant,
            only_spatial=True,
            start_flag=start_flag,
            T=T,
            k0=k0,
            target_volume=target_volume,
            thermal_sweeps=thermal_sweeps,
            main_sweeps=main_sweeps,
            chain=chain,
            run=run
        )

        measure_spectral_dimension(
            universe=universe,
            sigma_max=sigma_max,
            diffusion_constant=diffusion_constant,
            only_spatial=False,
            start_flag=start_flag,
            T=T,
            k0=k0,
            target_volume=target_volume,
            thermal_sweeps=thermal_sweeps,
            main_sweeps=main_sweeps,
            chain=chain,
            run=run
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--T', type=int, default=3)
    parser.add_argument('--k0', type=float, default=0.0)
    parser.add_argument('--target_volume', type=int, default=3000)
    parser.add_argument('--thermal_sweeps', type=int, default=1000)
    parser.add_argument('--main_sweeps', type=int, default=0)
    parser.add_argument('--chain', type=int, default=0)
    parser.add_argument('--sigma_max', type=int, default=1000)
    parser.add_argument('--diffusion_constant', type=float, default=0.8)
    parser.add_argument('--start_flag', type=int, default=0)
    parser.add_argument('--runs', type=int, default=1)
    args = parser.parse_args()

    main(args.T, args.k0, args.target_volume, args.thermal_sweeps, args.main_sweeps, args.chain, args.sigma_max, args.diffusion_constant, args.start_flag, args.runs)
    