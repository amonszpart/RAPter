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
mappingfile  = projectdir+'/corresp.csv'

linesfile_it1  = projectdir+'/primitives_it0.bonmin.csv'
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

# associative arrays, mapping
# - the gt primitive to the estimated primitive
# - the gt uid to the estimated uid
primitiveCorres, primitiveCorresId = packages.io.readPrimitiveCorrespondancesFromFiles(mappingfile, gtlines, lines_it1)

#gtlines    = packages.processing.removeUnassignedPrimitives(gtlines, gtassign)
#lines_it1  = packages.processing.removeUnassignedPrimitives(lines_it1, assign_it1)
#gtassign   = packages.processing.removeUnassignedPoint(gtlines, gtassign)
#assign_it1 = packages.processing.removeUnassignedPoint(lines_it1, assign_it1)

################################################################################
## Build relation graphs
print "Processing GT relations...."
gtGraph   = relgraph.RelationGraph(gtlines, gtassign)
print "Processing estimated relations...."
graph_it1 = relgraph.RelationGraph(lines_it1, assign_it1)

#[e['matched']=0 for e in gtGraph.G.edges_iter()]
#[e['matched']=0 for e in graph_it1.G.edges_iter()]
#[e[-1]['matched']=0 for e in graph_it1.G.edges_iter(data=True)]
for e in graph_it1.G.edges_iter(data=True):
    e[-1]['matched']=0
for e in gtGraph.G.edges_iter(data=True):
    e[-1]['matched']=0

for p in gtlines:
    p_node = gtGraph.G.node[p.uid]
    if not primitiveCorres.has_key(p):
        print "Gt Primitive not matched (",p.uid,",",p.did,")"
    else:
        matched_p = primitiveCorres[p]
        matched_p_node = graph_it1.G.node[matched_p.uid]
        
        # iterate over all relations and check they have a counterpart in the estimated scene
        # cUid is the uid of the connected component
        # cUid can be used to access the connected component using
        # print cUid, primitiveCorresId[cUid]
        for idx, cUid in enumerate(gtGraph.G.edge[p.uid]):
            # now we are look for the connection starting from matched_p et going to primitiveCorresId[cUid]
            # matched_p.uid, primitiveCorresId[cUid]
            # 
            # if we find it, we increment the matched field of the edges, and move to the next one
            #print cUid, primitiveCorresId[cUid]
            for idx2, matched_cUid in enumerate(graph_it1.G.edge[matched_p.uid]):                
                #print "  ",matched_cUid, primitiveCorresId[cUid]
                if matched_cUid == primitiveCorresId[cUid]:
                    #print "  match found !"
                    gtGraph.G.edge[p.uid][cUid]['matched'] += 1
                    graph_it1.G.edge[matched_p.uid][matched_cUid]['matched'] += 1
                    break
        
def checkEdges(graph, doPrint):
    correct=0
    error = False
    for e in graph.G.edges_iter(data=True):
        count = e[-1]['matched']
        if (count == 2):correct+=1
        elif count == 0: 
            if doPrint: print "Missed edge detected..."
        else: error = True;
        
    return correct, error

gtcorrect, gtError = checkEdges(gtGraph, True)
correct_it1, error_it1 = checkEdges(graph_it1, False)

if gtError or error_it1: 
    print "Error occurred, invalid number of matches. ABORT"
    quit()
    
if gtcorrect != correct_it1:
    print "Error: non-symmetric detection"
    quit()

print "precision=",float(gtcorrect)/float(graph_it1.G.number_of_edges())
print "recall   =",float(gtcorrect)/float(gtGraph.G.number_of_edges())
