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
parser = argparse.ArgumentParser(description='Compare ground truth noise distribution (continuous generator and generated samples) and the result of the optimisation.')
parser.add_argument('projectdir')

args = parser.parse_args()

projectdir = args.projectdir
if projectdir[-1] == '/':
    projectdir = projectdir[:-1]

projectname  = projectdir.split('/')[-1]
projectfile  = projectdir+'/gt/'+projectname+'.prj'
gtlinesfile  = projectdir+'/gt/primitives.csv'
gtassignfile = projectdir+'/gt/points_primitives.csv'
cloudfile    = projectdir+'/cloud.ply'

linesfile_it1  = projectdir+'/primitives_merged_it1.csv'
assignfile_it1 = projectdir+'/points_primitives_it1.csv'

print 'Processing project ', projectname

################################################################################
## Reading input files
project  = project.PyProject(projectfile)
cloud    = packages.io.readPointCloudFromPly(cloudfile)
gtlines  = primitive.readPrimitivesFromFile(gtlinesfile)
gtassign = packages.io.readPointAssignementFromFiles(gtassignfile)

lines_it1  = primitive.readPrimitivesFromFile(linesfile_it1)
assign_it1 = packages.io.readPointAssignementFromFiles(assignfile_it1)

gtlines    = packages.processing.removeUnassignedPrimitives(gtlines, gtassign)
lines_it1  = packages.processing.removeUnassignedPrimitives(lines_it1, assign_it1)
gtassign   = packages.processing.removeUnassignedPoint(gtlines, gtassign)
assign_it1 = packages.processing.removeUnassignedPoint(lines_it1, assign_it1)

################################################################################
## Build and display relation graphs
gtGraph   = relgraph.RelationGraph(gtlines, gtassign)
graph_it1 = relgraph.RelationGraph(lines_it1, assign_it1)

DiGM = isomorphism.GraphMatcher(gtGraph.G,graph_it1.G)
print DiGM.is_isomorphic()
print DiGM.mapping
print DiGM.subgraph_is_isomorphic()
print DiGM.mapping

print float(len(DiGM.mapping)) / float(len(gtlines))

setupGraphUI(gtGraph, 'Ground Truth')
setupGraphUI(graph_it1, 'Iteration 1')

plt.show()
