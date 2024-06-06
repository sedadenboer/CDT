import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Function to convert spherical coordinates to Cartesian coordinates
def spherical_to_cartesian(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z

# Create the graph
G = nx.Graph()
edges = [(0, 1), (0, 2), (0, 3), (1, 3), (1, 4), (1, 2), (2, 3), (2, 4), (3, 4)]
G.add_edges_from(edges)

# Define positions of points on the sphere using spherical coordinates
positions = {
    0: (1, 0, 0),
    1: (1, np.pi/2, 0),
    2: (1, np.pi/2, np.pi/2),
    3: (1, np.pi/2, np.pi),
    4: (1, np.pi, 0)
}

# Convert positions to Cartesian coordinates
cartesian_positions = {k: spherical_to_cartesian(*v) for k, v in positions.items()}

# Extract coordinates for plotting
x_nodes = [cartesian_positions[i][0] for i in G.nodes()]
y_nodes = [cartesian_positions[i][1] for i in G.nodes()]
z_nodes = [cartesian_positions[i][2] for i in G.nodes()]

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Draw nodes
ax.scatter(x_nodes, y_nodes, z_nodes, c='r', s=100)

# Draw edges
for edge in edges:
    x = [cartesian_positions[edge[0]][0], cartesian_positions[edge[1]][0]]
    y = [cartesian_positions[edge[0]][1], cartesian_positions[edge[1]][1]]
    z = [cartesian_positions[edge[0]][2], cartesian_positions[edge[1]][2]]
    ax.plot(x, y, z, color='b')

# Label nodes
for i in G.nodes():
    ax.text(cartesian_positions[i][0], cartesian_positions[i][1], cartesian_positions[i][2], str(i), color='black')

# Draw the sphere surface for better visualization
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 1 * np.outer(np.cos(u), np.sin(v))
y = 1 * np.outer(np.sin(u), np.sin(v))
z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='c', alpha=0.3, rstride=5, cstride=5)

# Set plot parameters
ax.set_box_aspect([1, 1, 1])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Spherical Topology Visualization')

plt.show()
