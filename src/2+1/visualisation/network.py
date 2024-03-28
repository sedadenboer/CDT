import sys
sys.path.append('..')
from classes.universe import Universe
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import plotly.graph_objects as go
from pyvis.network import Network
import bokeh.plotting as bp
from bokeh.models import MultiLine, Scatter, HoverTool, ColumnDataSource

def make_network_graph(universe, save=False, filename='network_graph.png'):
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
            scnum=vertex.scnum
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
    plt.figure(figsize=(10, 7))
    nx.draw_networkx_nodes(G, pos, node_size=20, node_color=node_colors, edgecolors='black', alpha=0.7)
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
    

if __name__ == '__main__':
    k0 = 7
    universe = Universe(geometry_infilename=f'saved_universes/N/measurements/swps=10_tswps=50_kstps=100000_trgtvol=10000_k0={k0}_chain=1_final.pkl.gz')
    initial_universe = Universe(geometry_infilename='/home/seda2102/epic/CDT/src/2+1/classes/initial_universes/sample-g0-T3.cdt')
    # make_network_graph_plotly(initial_universe)
    make_network_graph(universe, save=True, filename=f'plots/network_graph_k0={k0}')
    # make_network_graph(initial_universe, save=True, filename=f'plots/network_graph_k0={k0}')
    # make_network_pyvis(universe)
