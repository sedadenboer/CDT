# testing_plot.py
#
# Author: Seda den Boer
# Date: 03-01-2024
#
# Description:

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import vtk


def visualise_universe(grid):
    """
    Visualize the universe in 3D as a cylinder using Matplotlib.
    """
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')

    for time_slice in grid:
        theta = [2 * np.pi * node.ID / len(time_slice) for node in time_slice]
        z = [node.time_slice for node in time_slice]

        # Map spatial nodes onto a cylinder
        r = 1  # radius of the cylinder
        x = r * np.cos(theta)
        y = r * np.sin(theta)

        ax.scatter(x, y, z, label=f'Time Slice {time_slice[0].time_slice}', s=20, color='black')

        # Draw edges
        for i, node in enumerate(time_slice):
            if node.right_neighbour:
                ax.plot([x[i], x[(i + 1) % len(time_slice)]],
                        [y[i], y[(i + 1) % len(time_slice)]],
                        [z[i], z[(i + 1) % len(time_slice)]],
                        alpha=0.5, color='blue', linewidth=2)

                for future_node in node.future_nodes:
                    ax.plot([x[i], x[future_node.ID]],
                            [y[i], y[future_node.ID]],
                            [z[i], future_node.time_slice],
                            alpha=0.5, color='red', linewidth=2)

    ax.set_xlabel('X (Spatial Position)')
    ax.set_ylabel('Y (Spatial Position)')
    ax.set_zlabel('Time Slice')
    
    # Remove axes ticks x and y
    ax.set_xticks([])
    ax.set_yticks([])

    # Set time slice ticks in steps of 1
    ax.set_zticks(np.arange(0, len(grid), 1))

    plt.show()


# def createVTKDataObject(universe):
#     """Creates a VTK data object from the universe.

#     Args:
#         universe (Universe): The universe object.

#     Returns:
#         vtkPolyData: The VTK data object representing the universe.
#     """
#     import vtk
#     import math

#     pts = vtk.vtkPoints()
#     lines = vtk.vtkCellArray()
#     triangles = vtk.vtkCellArray()
#     colors = vtk.vtkUnsignedCharArray()
#     colors.SetNumberOfComponents(3)
#     colors.SetName("Colors")

#     # Get radius list
#     numPerTime = universe.getVolumePerTimeSlice()

#     # Get a node on t=0
#     pivot_node = None
#     for ID, node in universe.spatialNodes.items():
#         if node.t == 0:
#             pivot_node = node
#             break
    
#     # Create nodes
#     for i in range(universe.numTimeSlices):
#         numPts = numPerTime[i]
#         r = numPts / (2.0 * math.pi)
#         x0 = 0.0
#         y0 = 0.0
        
#         pivot_space = pivot_node
#         for j in range(numPts):
#             x = x0 + r * math.cos(2 * math.pi * j / numPts)
#             y = y0 + r * math.sin(2 * math.pi * j / numPts)
#             pts.InsertPoint(pivot_space.ID, (x, y, float(pivot_space.t)))
#             pivot_space = pivot_space.right_neighbour
        
#         if i < universe.numTimeSlices - 1:
#             pivot_node = pivot_node.future_nodes[0]
    
#     # Add edges
#     for ID, node in universe.spatialNodes.items():
#         if node.right_neighbour is not None:
#             line = vtk.vtkLine()
#             line.GetPointIds().SetId(0, node.ID) 
#             line.GetPointIds().SetId(1, node.right_neighbour.ID)
#             lines.InsertNextCell(line)
#             colors.InsertNextTuple3(0, 0, 0)
#         for e in node.future_nodes:
#             line = vtk.vtkLine()
#             line.GetPointIds().SetId(0, e.ID) 
#             line.GetPointIds().SetId(1, node.ID)
#             lines.InsertNextCell(line)
#             colors.InsertNextTuple3(0, 0, 0)
            
#     # Add triangles
#     for ID, node in universe.spatialNodes.items():
#         for e in node.past_nodes:
#             if e.right_neighbour is not None:
#                 if node in e.right_neighbour.future_nodes:
#                     triangle = vtk.vtkTriangle()
#                     triangle.GetPointIds().SetId(0, node.ID)
#                     triangle.GetPointIds().SetId(1, e.ID)
#                     triangle.GetPointIds().SetId(2, e.right_neighbour.ID)
#                     triangles.InsertNextCell(triangle)
#                     colors.InsertNextTuple3(255, 0, 0)

#         for e in node.future_nodes:
#             if e.right_neighbour is not None:
#                 if node in e.right_neighbour.past_nodes:
#                     triangle = vtk.vtkTriangle()
#                     triangle.GetPointIds().SetId(0, node.ID)
#                     triangle.GetPointIds().SetId(1, e.ID)
#                     triangle.GetPointIds().SetId(2, e.right_neighbour.ID)
#                     triangles.InsertNextCell(triangle)
#                     colors.InsertNextTuple3(0, 0, 255)
            
    
#     # Create a polydata to store everything in
#     linesPolyData = vtk.vtkPolyData()

#     linesPolyData.SetPoints(pts)
#     linesPolyData.SetLines(lines)
#     linesPolyData.SetPolys(triangles)
#     linesPolyData.GetCellData().SetScalars(colors)
#     linesPolyData.Modified()
    
#     return linesPolyData

# def createNewActor(universe):
#     """Creates a VTK actor from the universe.

#     Args:
#         universe (Universe): The universe object.

#     Returns:
#         vtkActor: The VTK actor representing the universe.
#     """
#     polyData = createVTKDataObject(universe)
    
#     polygonMapper = vtk.vtkPolyDataMapper()

#     polygonMapper.SetInputData(polyData)
#     polygonMapper.Update()

#     polygonActor = vtk.vtkActor()
#     polygonActor.SetMapper(polygonMapper)

#     return polygonActor

# def updateActor(render, polygonActor, universe):
#     """Updates the VTK actor representing the universe.

#     Args:
#         render (vtkRenderer): The VTK renderer.
#         polygonActor (vtkActor): The VTK actor representing the universe.
#         universe (Universe): The universe object.
#     """
#     render.RemoveActor(polygonActor)
#     polygonActor = createNewActor(universe)
#     render.AddActor(polygonActor)
