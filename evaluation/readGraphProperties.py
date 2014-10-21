import packages.project as project
import packages.primitive as primitive
import packages.processing
import packages.relationGraph as relgraph
import packages.io
import argparse
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import isomorphism
import numpy as np


################################################################################
## UI Generation
def setupGraphUI(graph, title):
    fig, ax1 = plt.subplots()
    fig.canvas.set_window_title(title)
    
    graph.draw()


################################################################################
## Command line parsing
parser = argparse.ArgumentParser(description='Simply read a graph an output its statistics.')
parser.add_argument('primitives')
parser.add_argument('point_primitives')
parser.add_argument('cloud')

args = parser.parse_args()

projectdir = args.primitives

linesfile_it1  = args.primitives
assignfile_it1 = args.point_primitives
cloud    = packages.io.readPointCloudFromPly(args.cloud)

################################################################################
## Reading input files

lines_it1  = primitive.readPrimitivesFromFile(linesfile_it1)
assign_it1 = packages.io.readPointAssignementFromFiles(assignfile_it1)

################################################################################
## Build and display relation graphs
graph_it1 = relgraph.RelationGraph(lines_it1, assign_it1)

print "Number of points:      ", len(cloud)
print "Number of primitives : ",graph_it1.G.number_of_nodes()
print "Number of connections: ",graph_it1.G.number_of_edges()
print "Max nb of connections: ",graph_it1.G.number_of_nodes()*graph_it1.G.number_of_nodes()

max_conn = graph_it1.G.number_of_nodes()* ((graph_it1.G.number_of_nodes()-1)/2)
print "Coverage: ",float(graph_it1.G.number_of_edges())/float(max_conn)


plt.show()
