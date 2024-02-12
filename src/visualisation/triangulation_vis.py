# triangulation_vis.py
#
# Author: Seda den Boer
# Date: 03-02-2024
#
# Description: This script visualises the triangulation.

import sys
sys.path.append('..')
from classes.universe import Universe
from classes.simulation import Simulation
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


def get_triangulation_matrix_periodic(universe: Universe, silent: bool = True):
    vertex_sheet = universe.get_triangulation_state()
    matrix = list(vertex_sheet.values())

    # Copy the first column item at the back of each row
    for row in matrix:
        row.append(row[0])

    # Copy the first row to the last row
    matrix.append(matrix[0])

    # Print the matrix with their IDs
    if not silent:
        print("Triangulation matrix:")
        for row in matrix:
            for vertex in row:
                print(vertex.ID, end=" ")
            print()

    return matrix

def make_graph(universe: Universe):
    # state = universe.get_triangulation_state()
    # vertex_sheet = list(state.values())
    vertex_sheet = get_triangulation_matrix_periodic(universe)
    G = nx.Graph()
    pos=nx.spring_layout(G)
    for y, row in enumerate(vertex_sheet):
        for x, col in enumerate(row):
            G.add_node(col.ID, pos=(x, y))

    # Add the edges from the space, future, and past neighbors in the vertex objects
    for y, row in enumerate(vertex_sheet):
        for x, col in enumerate(row):
            space_neighbours = [neighbour.ID for neighbour in col.get_space_neighbours()]
            future_neighbours = [neighbour.ID for neighbour in col.get_future_neighbours()]
            past_neighbours = [neighbour.ID for neighbour in col.get_past_neighbours()]
            # print(f"VERTEX {col.ID}: Space neighbours: {space_neighbours}, Future neighbours: {future_neighbours}, Past neighbours: {past_neighbours}")
            for neighbour in col.get_space_neighbours():
                G.add_edge(col.ID, neighbour.ID, color='blue', weight=2)
            for neighbour in col.get_future_neighbours():
                G.add_edge(col.ID, neighbour.ID, color='red', weight=1)
            # for neighbour in col.get_past_neighbours():
            #     G.add_edge(col.ID, neighbour.ID, color='red', weight=1)
         
     # Plot the graph with toroidal coordinates and colored edges
    pos = nx.get_node_attributes(G, 'pos')
    edges = G.edges()
    edge_colors = [G[u][v]['color'] for u, v in edges]
    weights = [G[u][v]['weight'] for u,v in edges]
    pos = nx.spring_layout(G)
    nx.draw(G, pos, with_labels=False, edge_color=edge_colors, node_size=1, node_color='black', width=weights)
    plt.show()

    return G

def plot_triangulation(universe: Universe):
    vertex_sheet = get_triangulation_matrix_periodic(universe, silent=True)
    fig, ax = plt.subplots()
    for y, row in enumerate(vertex_sheet):
        for x, col in enumerate(row):
            ax.plot(x, y, 'o', color='black')
            ax.text(x, y, str(col.ID), fontsize=12)

            # Color the last vertex differently in each row and also the full last row
            if y == len(vertex_sheet) - 1 or x == len(row) - 1:
                ax.plot(x, y, 'o', color='red')
                ax.text(x, y, str(col.ID), fontsize=12)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()

def plot_torus(universe):
    vertex_sheet = get_triangulation_matrix_periodic(universe)

    

if __name__ == "__main__":
    universe = Universe(10, 10)
    make_graph(universe)
    simulation = Simulation(universe, lambd=np.log(2), target_volume=0, epsilon=0)
    simulation.progress_universe(100, silence=False)
    make_graph(simulation.universe)
