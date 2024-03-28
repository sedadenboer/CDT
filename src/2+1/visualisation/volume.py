import sys
sys.path.append('..')
from classes.universe import Universe
from classes.simulation import Simulation
import matplotlib.pyplot as plt

def volume_profile(universe: Universe, save=False, filename='volume_profile.png'):
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
    plt.title('Volume profile of the universe')

    if save:
        plt.savefig(filename, dpi=400)
    plt.show()


if __name__ == '__main__':
    # k0 = 0
    # universe = Universe(geometry_infilename=f'saved_universes/N/measurements/swps=10_tswps=50_kstps=100000_trgtvol=10000_k0={k0}_chain=1_final.pkl.gz')
    # initial_universe = Universe(geometry_infilename='/home/seda2102/epic/CDT/src/2+1/classes/initial_universes/sample-g0-T3.cdt')

    universe = Universe(geometry_infilename='/home/seda2102/epic/CDT/src/2+1/classes/initial_universes/sample-g0-T3.cdt', strictness=3)
    # universe = Universe(geometry_infilename='/home/seda2102/epic/CDT/src/2+1/classes/initial_universes/output_g=0_T=32.txt', strictness=3)
    
    simulation = Simulation(universe)
    simulation.save_process = False
    simulation.save_final = False
    simulation.as_pickle = True
    simulation.start( 
        k0=7, k3=2.1, sweeps=0, thermal_sweeps=10, k_steps=100000,
        volfix_switch=0, target_volume=16000, target2_volume=0,
        seed=0, outfile="volume/swps=10_tswps=50_kstps=10000_trgtvol=16000_k0=7", validity_check=True,
        v1=1, v2=1, v3=1
    )
    # simulation.universe.make_network_graph()
    N3 = simulation.universe.tetrahedron_pool.size
    N31 = simulation.universe.tetras_31.size
    N22 = N3 - N31

    print(f"N3: {N3}, N31: {N31}, N22: {N22}")
    print("Done")
    volume_profile(universe, save=True, filename='plots/volume_profile')