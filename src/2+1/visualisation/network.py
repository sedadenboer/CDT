import sys
sys.path.append('..')
from classes.universe import Universe
import networkx as nx
from typing import List
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import plotly.graph_objects as go
from pyvis.network import Network
import numpy as np

def generate_network_graph(universe: Universe):
    # Create the network graph
    G = nx.Graph()
    for vertex in universe.vertex_pool.get_objects():
        G.add_node(
            vertex.ID,
            time=vertex.time,
            tetrahedron=vertex.tetra.ID,
            degree=len(universe.vertex_neighbours[vertex.ID]),
            cnum=vertex.cnum,
            scnum=vertex.scnum,
            )

    for vertex in universe.vertex_pool.get_objects():
        for neighbour_vertex in universe.vertex_neighbours[vertex.ID]:
            neighbour_vertex = universe.vertex_pool.get(neighbour_vertex)

            if vertex.time == neighbour_vertex.time:
                G.add_edge(vertex.ID, neighbour_vertex.ID)
            else:
                G.add_edge(vertex.ID, neighbour_vertex.ID)
    
    return G

def plot_network_graph(universe, save=False, filename='network_graph.png'):
    """
    Plots a network graph of the triangulation.
    """
    universe.log()
    
    node_colors_dict = {0: 'gold', 1: 'indigo', 2: 'cyan'}

    # Create the network graph
    G = nx.Graph()
    for vertex in universe.vertex_pool.get_objects():
        G.add_node(
            vertex.ID,
            color=node_colors_dict[vertex.time],
            time=vertex.time,
            tetrahedron=vertex.tetra.ID,
            degree=len(universe.vertex_neighbours[vertex.ID]),
            cnum=vertex.cnum,
            scnum=vertex.scnum,
            )

    for vertex in universe.vertex_pool.get_objects():
        for neighbour_vertex in universe.vertex_neighbours[vertex.ID]:
            neighbour_vertex = universe.vertex_pool.get(neighbour_vertex)

            if vertex.time == neighbour_vertex.time:
                G.add_edge(vertex.ID, neighbour_vertex.ID, color='b')
            else:
                G.add_edge(vertex.ID, neighbour_vertex.ID, color='r')

    # Draw the network graph
    pos = nx.spring_layout(G)
    edges = G.edges()
    edge_colors = [G.edges[edge]['color'] for edge in edges]
    node_colors = [G.nodes[node]['color'] for node in G.nodes()]
    node_sizes = [G.nodes[node]['degree'] for node in G.nodes()]
    plt.figure(figsize=(10, 7))
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, edgecolors='black', alpha=0.7)
    nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=1, alpha=0.7)

    # Add legend for nodes and edges
    edge_elements = [Line2D([0], [0], color='b', lw=1, label='Spacelike'),
                       Line2D([0], [0], color='r', lw=1, label='Timelike')]
    node_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=f't = {time}') for time, color in node_colors_dict.items()]

    plt.legend(handles=node_elements + edge_elements)

    if save:
        plt.savefig(filename + 'png', dpi=400)
        # Export as graph for gephi
        nx.write_gexf(G, f'{filename}.gexf')
    plt.show()

    return G

def make_network_graph_plotly(universe: Universe):
    """
    Plots an interactive network graph of the triangulation.
    """
    universe.log()
    
    node_colors_dict = {0: 'gold', 1: 'indigo', 2: 'cyan'}

    # Create the network graph
    G = nx.Graph()
    for vertex in universe.vertex_pool.get_objects():
        G.add_node(
            vertex.ID,
            color=node_colors_dict[vertex.time],
            time=vertex.time,
            tetrahedron=vertex.tetra.ID,
            degree=len(universe.vertex_neighbours[vertex.ID]),
            cnum=vertex.cnum,
            scnum=vertex.scnum
            )

    for vertex in universe.vertex_pool.get_objects():
        for neighbour_vertex in universe.vertex_neighbours[vertex.ID]:
            neighbour_vertex = universe.vertex_pool.get(neighbour_vertex)

            if vertex.time == neighbour_vertex.time:
                G.add_edge(vertex.ID, neighbour_vertex.ID, color='blue')
            else:
                G.add_edge(vertex.ID, neighbour_vertex.ID, color='red')

    pos = nx.spring_layout(G)
    
    # Make plotly graph
    colors = [G.edges[edge]['color'] for edge in G.edges()]
    edge_traces = []
    for edge_set, color in zip(G.edges(), colors):
        edge_x = []
        edge_y = []
        for edge in edge_set:
            x0, y0 = G.nodes[edge[0]]['pos']
            x1, y1 = G.nodes[edge[1]]['pos']
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])

        edge_traces.append(go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color=color),
            hoverinfo='none',
            mode='lines'))
        
    # Add node info to hovertext
    node_hovertext = []
    for node in G.nodes():
        vertex = universe.vertex_pool.get(node)
        hovertext = f'ID: {node}\n'
        hovertext += f'time: {vertex.time}\n'
        hovertext += f'tetrahedron: {vertex.tetra.ID}\n'
        hovertext += f'degree: {len(universe.vertex_neighbours[node])}\n'
        hovertext += f'cnum: {vertex.cnum}\n'
        hovertext += f'scnum: {vertex.scnum}'
        node_hovertext.append(hovertext)
    
    node_x = []
    node_y = []

    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        hovertext=node_hovertext,
        marker=dict(
            showscale=True,
            colorscale='viridis',
            reversescale=True,
            color=[],
            size=10,
            colorbar=dict(
                thickness=15,
                title='Time',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))

    node_adjacencies = []
    node_text = []
    for node in G.nodes():
        node_text.append(f'ID: {node}')
        node_adjacencies.append(universe.vertex_pool.get(node).time)

    node_trace.marker.color = node_adjacencies
    node_trace.text = node_text

    fig = go.Figure(data=edge_traces + [node_trace],
                layout=go.Layout(
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=0,l=0,r=0,t=0),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)
                ))
    
    fig.show()

def make_network_pyvis(universe: Universe):
    universe.log()

    # Create a pyvis network
    nt = Network(height='500px', width='100%')

    # Generate n colors based on the times
    colors = ['red', 'blue', 'green', 'purple', 'orange', 'yellow', 'pink', 'brown', 'black', 'grey']
    color_map = {time: colors[time] for time in range(10)}
    
    # Add nodes
    for vertex in universe.vertex_pool.get_objects():
        nt.add_node(
                vertex.ID,
                label=vertex.ID,
                title=f'ID: {vertex.ID},\n\
                        time: {vertex.time},\n\
                        tetrahedron: {vertex.tetra.ID},\n\
                        degree: {len(universe.vertex_neighbours[vertex.ID])},\n\
                        cnum: {vertex.cnum},\n\
                        scnum: {vertex.scnum}',
                color=color_map[vertex.time],
                size=10,
                )
    
    # Add edges
    for vertex in universe.vertex_pool.get_objects():
        for neighbour_vertex in universe.vertex_neighbours[vertex.ID]:
            nt.add_edge(
                vertex.ID,
                neighbour_vertex,
                color='blue' if universe.vertex_pool.get(neighbour_vertex).time == vertex.time else 'red'
                )

    nt.show_buttons()
    nt.show('network.html', notebook=False)

def pick_from_max_spatial_volume(universe: Universe):
    """
    Picks n tetrahedra with the maximum spatial volume.
    """
    # Find which timeslice has the maximum spatial volume
    T = universe.slice_sizes
    max_spatial_volume = max(T)
    max_spatial_volume_timeslice = T.index(max_spatial_volume)

    # Pick a random tetrahedron
    tetrahedra_id = universe.tetrahedron_pool.pick()
    tetrahedra = universe.tetrahedron_pool.get(tetrahedra_id)

    # Keep picking until we have a tetrahedron at the timeslice with the maximum spatial volume
    while tetrahedra.time != max_spatial_volume_timeslice:
        tetrahedra_id = universe.tetrahedron_pool.pick()
        tetrahedra = universe.tetrahedron_pool.get(tetrahedra_id)

    return tetrahedra

def spectral_dimension(universe: Universe, diffusion_time: int = 1000, n_walkers: int = 1000):
    """
    Calculates the spectral dimension of the universe. Does 
    this by performing a random walk on the universe and
    counting the number of times the walker returns to the
    original tetrahedron.
    """
    n_returned = 0

    # # Start the random walkers
    # for _ in range(n_walkers):
    #     # Start the walker at a random tetrahedron in a dense region
    #     walker = pick_from_max_spatial_volume(universe)
    #     walker_id = walker.ID

    #     # Perform the random walk
    #     for _ in range(diffusion_time):
    #         neighbours = walker.get_tetras()
    #         walker = np.random.choice(neighbours)

    #     # Check if the walker returned to the original tetrahedron
    #     if walker.ID == walker_id:
    #         n_returned += 1

    # Start the random walkers
    for _ in range(n_walkers):
        # Start the walker at a random vertex and save it
        walker = universe.vertex_pool.pick()
        walker_i = walker

        # Perform the random walk
        for _ in range(diffusion_time):
            neighbours = universe.vertex_neighbours[walker]
            walker = np.random.choice(neighbours)

        # Check if the walker returned to the original vertex
        if walker == walker_i:
            n_returned += 1
    
    # Calculate the spectral dimension
    return_probability = np.float64(n_returned) / np.float64(n_walkers)
    spectral_dimension = -2 * np.log(return_probability) / np.log(diffusion_time)
    print(f'Time: {diffusion_time}, Return probability: {return_probability}, Spectral dimension: {spectral_dimension}\n')

    return spectral_dimension

def plot_spectral_dimension(universe: Universe, save=False, filename='spectral_dimension.png'):
    """
    Plots the spectral dimension of the universe.
    """
    step = 50
    diffusion_times = np.arange(step, 400 + step, step)
    n_walkers = 10000
    spectral_dimensions = []

    # Retrieve the spectral dimension for each diffusion time
    for diffusion_time in diffusion_times:
        spectral_dimensions.append(spectral_dimension(universe, diffusion_time, n_walkers))
    
    # Save spectral dimensions as a txt 
    np.savetxt('saved_universes/spectral_dimensions.txt', spectral_dimensions)
    
    # Plot the spectral dimension
    plt.figure(figsize=(10, 7))
    plt.plot(diffusion_times, spectral_dimensions, marker='o')
    plt.xlabel('$\sigma$')
    plt.ylabel('$D_s$')
    if save:
        plt.savefig(filename, dpi=400)
    plt.show()
    

if __name__ == '__main__':
    k0 = 0
    universe = Universe(geometry_infilename=f'saved_universes/N/measurements/50_50/swps=50_tswps=50_kstps=100000_trgtvol=10000_k0={k0}_chain=0_final.pkl.gz')
    # initial_universe = Universe(geometry_infilename='/home/seda2102/epic/CDT/src/2+1/classes/initial_universes/sample-g0-T3.cdt')
    # make_network_graph_plotly(initial_universe)
    # plot_network_graph(universe, save=True, filename=f'plots/network_graph_k0={k0}')
    # make_network_graph(initial_universe, save=True, filename=f'plots/network_graph_k0={k0}')
    # make_network_pyvis(universe)

    plot_spectral_dimension(universe, save=False, filename=f'plots/spectral_dimension_k0={k0}.png')
