from pyvis.network import Network
import networkx as nx

# Create a graph using NetworkX
G = nx.cycle_graph(5)

# Convert to Pyvis for visualization
net = Network(notebook=True)
net.from_nx(G)
net.show("graph.html")
