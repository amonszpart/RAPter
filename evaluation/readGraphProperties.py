import packages.project as project
import packages.primitive as primitive
import packages.processing
import packages.relationGraph as relgraph
import packages.io
import packages.colours as colours
import argparse
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import isomorphism
import numpy as np

try:
    from networkx import graphviz_layout
    layout=nx.graphviz_layout
except ImportError:
    print("PyGraphviz not found; drawing with spring layout; will be slow.")
    layout=nx.spring_layout
    
################################################################################
## UI Generation
def setupGraphUI(graph, primitives, title):
    fig, ax1 = plt.subplots()
    fig.canvas.set_window_title(title)
    
    lay = layout(graph.G)

    CmapObject = colours.Colours()
    cmap, gfilter = CmapObject.getDIDColourMap(primitives)
    
    print cmap
    print gfilter
    
    #nx.draw(graph.G, lay)
    nx.draw_networkx_edges(graph.G, lay)
    for did, colour in cmap.iteritems():    
        nx.draw_networkx_nodes(graph.G, lay, node_size=800, nodelist=gfilter[did], node_color=colour)
    nx.draw_networkx_labels(graph.G, lay)

angles = [0., 60., 90., 120., 180.]
tolerance = 0.1

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
graph_it1 = relgraph.RelationGraph(lines_it1, assign_it1, angles, tolerance)

print "Number of points:      ", len(cloud)
print "Number of primitives : ",graph_it1.G.number_of_nodes()
print "Number of connections: ",graph_it1.G.number_of_edges()
print "Max nb of connections: ",graph_it1.G.number_of_nodes()*graph_it1.G.number_of_nodes()

max_conn = graph_it1.G.number_of_nodes()* ((graph_it1.G.number_of_nodes()-1)/2)
print "Coverage: ",float(graph_it1.G.number_of_edges())/float(max_conn)

setupGraphUI(graph_it1, lines_it1, "Graph")


plt.show()
