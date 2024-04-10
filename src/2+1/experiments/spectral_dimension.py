# Load packages
import sys
sys.path.append('..')
from typing import Dict, List
import numpy as np
import gzip
import pickle
from multiprocessing import Pool
from classes.universe import Universe
from classes.helper_functions.helpers import get_spatial_neighbours, get_vertices_in_slice


def pick_from_max_spatial_volume(universe: Universe):
    """
    Picks n tetrahedra with the maximum spatial volume.
    """
    # Find which timeslice has the maximum spatial volume
    slice_sizes = universe.slice_sizes
    max_spatial_volume_timeslice = max(slice_sizes, key=lambda k: len(slice_sizes[k]) if isinstance(slice_sizes[k], list) else -1)

    # Pick a random tetrahedron
    tetrahedra_id = universe.tetrahedron_pool.pick()
    tetrahedra = universe.tetrahedron_pool.get(tetrahedra_id)

    # Keep picking until we have a tetrahedron at the timeslice with the maximum spatial volume
    while tetrahedra.time != max_spatial_volume_timeslice:
        tetrahedra_id = universe.tetrahedron_pool.pick()
        tetrahedra = universe.tetrahedron_pool.get(tetrahedra_id)

    return tetrahedra

def spectral_dimension(universe: Universe, diffusion_times: List[int], n_walkers: int = 1000):
    """
    Calculates the spectral dimension of the universe. Does 
    this by performing a random walk on the universe and
    counting the number of times the walker returns to the
    original tetrahedron.
    """
    n_returned = 0
    spectral_dimensions = []

    # Start the random walkers
    for diffusion_time in diffusion_times:
        for _ in range(n_walkers):
            # Start the walker at a random tetrahedron in a dense region
            walker = pick_from_max_spatial_volume(universe)
            walker_id = walker.ID

            # Perform the random walk
            for _ in range(diffusion_time):
                neighbours = walker.get_tetras()
                walker = np.random.choice(neighbours)

            # Check if the walker returned to the original tetrahedron
            if walker.ID == walker_id:
                n_returned += 1
    
        # Calculate the spectral dimension
        return_probability = np.float64(n_returned) / np.float64(n_walkers)
        spectral_dimension = -2 * np.log(return_probability) / np.log(diffusion_time)
        spectral_dimensions.append(spectral_dimension)
        print(f'Time: {diffusion_time}, Return probability: {return_probability}, Spectral dimension: {spectral_dimension}\n')

    return spectral_dimension

def spectral_dimension_per_timeslice(universe: Universe, diffusion_times: List[int], n_walkers: int = 1000) -> List[float]:
    """
    Calculates the spectral dimension of the universe per time slice. Does this by performing
    a random walk on the timeslice and counting the number of times the walker returns to the
    original tetrahedron.
    """
    T = universe.n_slices
    spectral_dimensions = {t: [] for t in range(T)}
    vertices_in_slice = get_vertices_in_slice(universe)
    spatial_neighbours = get_spatial_neighbours(universe)

    # For each timeslice
    for t in range(T):
        for diffusion_time in diffusion_times:
            n_returned = 0
            for _ in range(n_walkers):
                # Pick random initial walker in each timeslice
                walker = np.random.choice(vertices_in_slice[t])
                walker_i = walker

                # Perform the random walk
                for _ in range(diffusion_time):
                    neighbours = spatial_neighbours[walker]
                    walker = np.random.choice(neighbours)

                # Check if the walker returned to the original vertex
                if walker == walker_i:
                    n_returned += 1
            
            return_probability = np.float64(n_returned) / np.float64(n_walkers)
            spectral_dimension = -2 * np.log(return_probability) / np.log(diffusion_time)
            spectral_dimensions[t].append(spectral_dimension)
            print(f'Timeslice: {t}, Time: {diffusion_time}, Return probability: {return_probability}, Spectral dimension: {spectral_dimension}\n')
    
    # Save spectral dimensions
    output_filename = f'measurements/spectral_dimensions_per_timeslice_T={T}'
    with gzip.open(output_filename + '.pkl.gz', 'wb', compresslevel=6) as f:
        pickle.dump(spectral_dimensions, f)

    return spectral_dimensions


if __name__ == "__main__":
    # Load files
    k0_values = [0, 1, 2, 3, 4, 5, 6, 7]
    k3_init_guess = [1.0417799999999955, 1.1760799999999827, 1.3237199999999882, 1.4782399999999793, 1.6480799999999842, 1.8284799999999846, 2.0621400000000256, 2.3117000000000187]
    filenames = [f'saved_universes/thermal_1000/universe_k0={k0}_tswps=1000_swps=0_kstps=1000000_trgtv=10000_trgtv2=0_fx=0_chn=0_thermal_1000.pkl.gz' for k0 in k0_values]

    universe_0 = Universe(geometry_infilename=filenames[0])
    step = 50
    diffusion_times = np.arange(step, 400 + step, step)
    n_walkers = 10000

    result_per_timeslice = spectral_dimension_per_timeslice(universe_0, diffusion_times, n_walkers)
    # result = spectral_dimension(universe_0, diffusion_times, n_walkers)
