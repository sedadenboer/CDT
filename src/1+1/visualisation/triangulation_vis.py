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
import matplotlib.pyplot as plt
import vtk
import networkx as nx
import copy

def torus_vtk(universe: Universe):
    """
    Visualise the triangulation of a 1+1D torus using VTK.

    Args:
        universe (Universe): Universe to visualise.
    """
    state = universe.get_triangulation_state()
    vertex_sheet = list(state.values())

    # Create vtk object from the universe
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    triangles = vtk.vtkCellArray()
    colors = vtk.vtkUnsignedCharArray()
    colors.SetNumberOfComponents(3)
    colors.SetName("Colors")

    # Get radius list 
    n_per_time = [len(row) for row in vertex_sheet]
    print(n_per_time)
    print(len(vertex_sheet))

    # Get a node
    all_nodes_indices = universe.vertex_pool.used_indices
    all_nodes = [universe.vertex_pool.elements[i] for i in all_nodes_indices]
    pivot_node = all_nodes[0]

    # Create the torus points
    for i in range(universe.total_time):
        num_points = n_per_time[i]
        for j in range(num_points):
            # Calculate x, y, z coordinates for each point
            radius = num_points / (2 * np.pi)
            theta = 2 * np.pi * j / num_points
            phi = 2 * np.pi * i / len(n_per_time)
            x = (radius + 10 + radius * np.cos(theta)) * np.cos(phi)
            y = (radius + 10 + radius * np.cos(theta)) * np.sin(phi)
            z = radius * np.sin(theta)

            # Add the point to the vtk object
            points.InsertPoint(pivot_node.ID, (x, y, z))
            pivot_node = pivot_node.get_neighbour_right()
        
        # Move to the next time
        if i < universe.total_time - 1:
            pivot_node = pivot_node.get_future_neighbours()[0]

    # Create the lines
    for node in all_nodes:
        # Add right neighbours
        right_neighbour = node.get_neighbour_right()
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, node.ID)
        line.GetPointIds().SetId(1, right_neighbour.ID)
        lines.InsertNextCell(line)
        colors.InsertNextTuple3(0, 0, 0)

        for future_node in node.get_future_neighbours():
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, future_node.ID)
            line.GetPointIds().SetId(1, node.ID)
            lines.InsertNextCell(line)
            colors.InsertNextTuple3(0, 0, 0)
    
    # Add triangles
    all_triangles_indices = universe.triangle_pool.used_indices
    all_triangles = [universe.triangle_pool.elements[i] for i in all_triangles_indices]
    
    for triangle in all_triangles:
        left, right, center = triangle.get_vertices()
        triangle_vtk = vtk.vtkTriangle()
        triangle_vtk.GetPointIds().SetId(0, left.ID)
        triangle_vtk.GetPointIds().SetId(1, right.ID)
        triangle_vtk.GetPointIds().SetId(2, center.ID)
        triangles.InsertNextCell(triangle_vtk)
        if triangle.is_upwards():
            colors.InsertNextTuple3(255, 0, 0)
        else:
            colors.InsertNextTuple3(0, 0, 255)

    # Create a polydata object to store everything
    linesPolyData = vtk.vtkPolyData()
    linesPolyData.SetPoints(points)
    linesPolyData.SetLines(lines)
    linesPolyData.SetPolys(triangles)
    linesPolyData.GetCellData().SetScalars(colors)
    linesPolyData.Modified()

    # Create a mapper and actor
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(linesPolyData)
    mapper.Update()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    axes_actor = vtk.vtkAxesActor()
    axes = vtk.vtkOrientationMarkerWidget()
    axes.SetOrientationMarker(axes_actor)

    # Create a renderer, render window, and interactor
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(0, 0, 0)
    renderer.ResetCamera()

    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(renderer)
    renWin.SetSize(450, 350)
    renWin.SetWindowName("Causal Dynamical Triangulation 1+1D")

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
    iren.SetRenderWindow(renWin)
    axes.SetInteractor(iren)
    axes.EnabledOn()
    axes.InteractiveOn()

    iren.Initialize()
    iren.Start()

def get_triangulation_matrix_periodic(universe: Universe, silent: bool = True) -> list:
    """
    Get the triangulation matrix for a periodic universe.

    Args:
        universe (Universe): Universe to get the triangulation matrix from.
        silent (bool, optional): If matrix should be printed. Defaults to True.

    Returns:
        list: periodic triangulation matrix
    """
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
    
def plot_triangulation_flat(universe: Universe):
    """
    Plot the triangulation network.

    Args:
        universe (Universe): Universe to plot.
    """
    state = universe.get_triangulation_state()
    vertex_sheet = list(state.values())

    # Al triangles
    all_triangles_indices = universe.triangle_pool.used_indices
    all_triangles = [universe.triangle_pool.elements[i] for i in all_triangles_indices]

    # Get a triangle at time = 0
    row_0 = vertex_sheet[0]
    vertex_t0 = row_0[0]
    triangle_t0 = vertex_t0.get_triangle_right()

    # Add to graph the first triangle
    G = nx.Graph()
    for vertex in triangle_t0.get_vertices():
        G.add_node(vertex.ID)
    G.add_edge(triangle_t0.vl_.ID, triangle_t0.vr_.ID)
    G.add_edge(triangle_t0.vr_.ID, triangle_t0.vc_.ID)
    G.add_edge(triangle_t0.vc_.ID, triangle_t0.vl_.ID)
    
    for i in range(len(all_triangles)):
        triangle = all_triangles[i]
        for vertex in triangle.get_vertices():
            G.add_node(vertex.ID)
        G.add_edge(triangle.vl_.ID, triangle.vr_.ID)
        G.add_edge(triangle.vr_.ID, triangle.vc_.ID)
        G.add_edge(triangle.vc_.ID, triangle.vl_.ID)
        
    # Plot the graph
    pos = nx.spring_layout(G, dim=2)
    nx.draw(G, pos, with_labels=False, node_size=5, node_color='skyblue', font_size=8, font_weight='bold', edge_color='black', linewidths=1, width=1, alpha=0.7)
    plt.show()


if __name__ == "__main__":
    universe = Universe(40, 50)

    simulation = Simulation(universe, lambd=np.log(2))
    simulation.progress_universe(1000, silence=True)

    universe.print_state()

    torus_vtk(universe)

