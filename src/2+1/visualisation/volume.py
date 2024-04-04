import sys
sys.path.append('..')
from classes.universe import Universe
from classes.simulation import Simulation
from typing import Dict
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 16
import gzip
import pickle
import numpy as np
import seaborn as sns
import pandas as pd

def volume_profile(universe: Universe, k0: int, savefig=False, filename='volume_profile.png'):
    """
    Plots the volume profile of the universe.
    """
    # universe.log()

    # Make a dict for each time slice with the number of vertices
    volumes = universe.slice_sizes
    T = universe.n_slices
    times = list(range(T))
    print(volumes)
    print(times)
    # Plot the volume profile
    plt.figure(figsize=(10, 7))
    plt.bar(times, volumes, color='b', alpha=0.7)
    plt.xlabel('Time')
    plt.ylabel('$n(t)$')
    plt.title(f'Volume profile of a universe with $k_0 = {k0}$')

    # TODO: make the plot look better (cocoon triangulation shape)

    if savefig:
        plt.savefig(filename, dpi=400)

def compare_volume_profile(savefig=False, filename='volume_profile_comparison.png'):
    """
    Compares the volume profiles of multiple universes.
    """
    # Load the universes
    chain = 0
    filename1 = f'saved_universes/N/measurements/50_50/swps=50_tswps=50_kstps=100000_trgtvol=10000_k0=0_chain={chain}_final.pkl.gz'
    filename2 = f'saved_universes/N/measurements/50_50/swps=50_tswps=50_kstps=100000_trgtvol=10000_k0=3_chain={chain}_final.pkl.gz'
    filename3 = f'saved_universes/N/measurements/50_50/swps=50_tswps=50_kstps=100000_trgtvol=10000_k0=7_chain={chain}_final.pkl.gz'
    universe_1 = Universe(geometry_infilename=filename1, strictness=3)
    universe_2 = Universe(geometry_infilename=filename2, strictness=3)
    universe_3 = Universe(geometry_infilename=filename3, strictness=3)

    # Get volume data
    volumes_1 = universe_1.slice_sizes
    volumes_2 = universe_2.slice_sizes
    volumes_3 = universe_3.slice_sizes

    # Combine the datasets
    data = {'time': [], 'volume': [], 'k0': []}
    volumes = [volumes_1, volumes_2, volumes_3]
    k0_values = [0, 3, 7]
    for i, volumes_i in enumerate(volumes):
        k0 = k0_values[i]
        for t, volume in enumerate(volumes_i):
            data['time'].append(t)
            data['volume'].append(volume)
            data['k0'].append(k0)

    # Plot the data
    plt.figure(figsize=(10, 7))
    sns.barplot(x='time', y='volume', hue='k0', data=pd.DataFrame(data), palette='viridis')
    plt.xlabel('Time')
    plt.ylabel('$n(t)$')
    plt.title('Volume profile comparison')
    plt.legend(title='$k_0$')
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.gca().set_axisbelow(True)
    if savefig:
        plt.savefig(filename, dpi=400)
    plt.show()

def curvature_profile(universe: Universe):
    """
    Gets the curvature profile (spatial coordination number) of a universe per time slice.
    """
    T = universe.n_slices
    times = list(range(T))

    # Make a dict for each time slice with the vertices
    scnums_per_slice = {t: [] for t in times}

    for v in universe.vertex_pool.get_objects():
        scnums_per_slice[v.time].append(v.scnum)

    return scnums_per_slice

def plot_curvature_profile(universe: Universe, k0: int, savefig=False, filename='curvature_profile.png'):
    """
    Plots the curvature profile of the universe per time slice,
    as a violin plot.
    """
    # Get curvature data
    scnums_per_slice = curvature_profile(universe)
    T = universe.n_slices
    times = list(range(T))

    # Plot
    plt.figure(figsize=(10, 7))
    sns.violinplot(data=[scnums_per_slice[t] for t in times], palette='viridis')
    plt.xlabel('Time')
    plt.ylabel('Spatial coordination number')
    plt.title(f'Curvature profile ($k_0 = {k0}$)')
    plt.xticks(times)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.gca().set_axisbelow(True)
    if savefig:
        plt.savefig(filename, dpi=400)
    plt.show()

def compare_curvature_profile(savefig=False, filename='curvature_profile_comparison.png'):
    """
    Compares the curvature profiles of multiple universes.
    """
    # Load the universes
    chain = 0
    filename1 = f'saved_universes/N/measurements/50_50/swps=50_tswps=50_kstps=100000_trgtvol=10000_k0=0_chain={chain}_final.pkl.gz'
    filename2 = f'saved_universes/N/measurements/50_50/swps=50_tswps=50_kstps=100000_trgtvol=10000_k0=3_chain={chain}_final.pkl.gz'
    filename3 = f'saved_universes/N/measurements/50_50/swps=50_tswps=50_kstps=100000_trgtvol=10000_k0=7_chain={chain}_final.pkl.gz'
    universe_1 = Universe(geometry_infilename=filename1, strictness=3)
    universe_2 = Universe(geometry_infilename=filename2, strictness=3)
    universe_3 = Universe(geometry_infilename=filename3, strictness=3)

    # Get curvature data
    scnums_per_slice_1 = curvature_profile(universe_1)
    scnums_per_slice_2 = curvature_profile(universe_2)
    scnums_per_slice_3 = curvature_profile(universe_3)

    # Combine the datasets
    data = {'time': [], 'scnum': [], 'k0': []}
    scnums_per_slice = [scnums_per_slice_1, scnums_per_slice_2, scnums_per_slice_3]
    k0_values = [0, 3, 7]
    for i, scnums_per_slice_i in enumerate(scnums_per_slice):
        k0 = k0_values[i]
        for t, scnums in scnums_per_slice_i.items():
            for scnum in scnums:
                data['time'].append(t)
                data['scnum'].append(scnum)
                data['k0'].append(k0)

    # Plot the data
    plt.figure(figsize=(10, 7))
    sns.violinplot(x='time', y='scnum', hue='k0', data=pd.DataFrame(data), palette='viridis')
    plt.xlabel('Time')
    plt.ylabel('Spatial coordination number')
    plt.title('Curvature profile comparison')
    plt.legend(title='$k_0$')
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.gca().set_axisbelow(True)
    if savefig:
        plt.savefig(filename, dpi=400)
    plt.show()
   

if __name__ == '__main__':
    k0 = 0
    chain = 1
    initial_universe = '../classes/initial_universes/sample-g0-T3.cdt'
    filename = f'saved_universes/N/measurements/50_50/swps=50_tswps=50_kstps=100000_trgtvol=10000_k0={k0}_chain={chain}_final.pkl.gz'
    # universe = Universe(geometry_infilename=filename , strictness=3)
    
    # plot_curvature_profile(universe=universe, k0=k0, savefig=True, filename=f'plots/curvature_profile_k0={k0}_chain={chain}.png')

    compare_curvature_profile(savefig=True, filename='plots/curvature_profile_comparison_N31_fixed.png')
    compare_volume_profile(savefig=True, filename='plots/volume_profile_comparison_N31_fixed.png')