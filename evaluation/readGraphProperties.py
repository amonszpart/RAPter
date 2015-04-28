import packages.project as project
import packages.primitive as primitive
import packages.processing as processing
import packages.relationGraph as relgraph
import packages.io
import packages.utils as utils
import packages.colours as colours
import packages.orderedSet as orderedSet
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
    
    nx.draw_networkx_edges(graph.G, lay)
    for did, colour in cmap.iteritems():    
        nx.draw_networkx_nodes(graph.G, lay, node_size=800, nodelist=gfilter[did], node_color=colour)
    nx.draw_networkx_labels(graph.G, lay)

tolerance = 0.1

################################################################################
## Command line parsing
parser = argparse.ArgumentParser(description='Simply read a graph an output its statistics.')
parser.add_argument('primitives')
parser.add_argument('point_primitives')
parser.add_argument('cloud')
parser.add_argument('--angles', nargs='*')
parser.add_argument('--iteration', default=" unknown")

args = parser.parse_args()

projectdir = args.primitives
angles     = utils.parseAngles(args.angles)
itId       = args.iteration

linesfile  = args.primitives
assignfile = args.point_primitives
cloud      = packages.io.readPointCloudFromPly(args.cloud)

################################################################################
## Reading input files

lines  = primitive.readPrimitivesFromFile(linesfile)
assign = packages.io.readPointAssignementFromFiles(assignfile)

#cleanlines = processing.removeUnassignedPrimitives(lines, assign)

################################################################################
## Build and display relation graphs
graph = relgraph.RelationGraph(lines, assign, angles, tolerance)

print "Number of points:      ", len(cloud)
print "Number of primitives : ",graph.G.number_of_nodes()
print "Number of metanodes:   ", graph.getNumberOfMetanodes()
print "Number of N to N rels: ", graph.getNumberOfNodeToNodeRelations()
#print "Number of connections: ",graph.G.number_of_edges()
#print "Max nb of connections: ",graph.G.number_of_nodes()*graph.G.number_of_nodes()

exit()

setupGraphUI(graph, lines, "Iteration "+itId)

plt.savefig("relationGraphs_it"+itId+".svg")
