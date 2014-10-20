import packages.project as project
import packages.primitive as primitive
import packages.utils as utils
import packages.processing
import packages.io
import argparse
from matplotlib import pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm,kstest,skewtest,kurtosistest,normaltest
import os.path
import unitTests.relations as test_relations
import packages.relationGraph as relgraph

compareToGt = False

itmin=0
itmax=20



################################################################################
## Command line parsing
parser = argparse.ArgumentParser(description='Collect full run statistics.')
parser.add_argument('projectdir')

args = parser.parse_args()

projectdir = args.projectdir
if projectdir[-1] == '/':
    projectdir = projectdir[:-1]

projectname  = projectdir.split('/')[-1]
cloudfile    = projectdir+'/cloud.ply'
cloud        = packages.io.readPointCloudFromPly(cloudfile)

print 'Processing project ', projectname

projectfile  = projectdir+'/gt/'+projectname+'.prj'
gtlinesfile  = projectdir+'/gt/primitives.csv'
gtassignfile = projectdir+'/gt/points_primitives.csv'

compareToGt = os.path.isfile(projectfile) and os.path.isfile(gtlinesfile) and os.path.isfile(gtassignfile)

sigma_ref = 1.

if compareToGt:
    
    project  = project.PyProject(projectfile)
    gtlines  = primitive.readPrimitivesFromFile(gtlinesfile)
    gtassign = packages.io.readPointAssignementFromFiles(gtassignfile)
    
    gtlines  = packages.processing.removeUnassignedPrimitives(gtlines, gtassign)
    gtassign = packages.processing.removeUnassignedPoint(gtlines, gtassign)
    
    gtgraph  = relgraph.RelationGraph(gtlines, gtassign)
    
    ############################################################################
    ## Process noise
    sigma_ref = project.kernels[0].stdev()
    gtDistrib   = utils.distanceToPrimitives(cloud, gtassign,   gtlines)
    print "Ground truth sigma ", np.sqrt(np.var(gtDistrib)) / sigma_ref    
    print "# it_count sigma F-measure precision recall"
else:
    print "# it_count sigma"



for it in range(itmin, itmax):
    linesfile_it  = projectdir+'/primitives_it'+str(it)+'.bonmin.csv'
    assignfile_it = projectdir+'/points_primitives_it'+str(it-1)+'.csv'
    
    if it == 0:
        assignfile_it = projectdir+'/points_primitives.csv'
        
        
    if not os.path.isfile(linesfile_it):
        break    
    
    #print "Processing iteration",it
    #print "   ",linesfile_it
    #print "   ",assignfile_it


    ################################################################################
    ## Reading input files
    lines_it  = primitive.readPrimitivesFromFile(linesfile_it)
    assign_it = packages.io.readPointAssignementFromFiles(assignfile_it)
    
    #lines_it  = packages.processing.removeUnassignedPrimitives(lines_it, assign_it)
    #assign_it = packages.processing.removeUnassignedPoint(lines_it, assign_it)
    

    ################################################################################
    ## Process noise

    # Compute the distance between each point and its assigned primitive
    distrib_it = [x for x in utils.distanceToPrimitives(cloud, assign_it, lines_it)  if x != []]
    #print distrib_it    
    #distrib_it = distrib_it[np.logical_not(np.isnan(distrib_it))]
    sigma = np.sqrt(np.var(distrib_it)) / sigma_ref


    if compareToGt:
        mappingfile  = projectdir+'/primitives_corresp_it'+str(it)+'.csv'
        primitiveCorres, primitiveCorresId = packages.io.readPrimitiveCorrespondancesFromFiles(mappingfile, gtlines, lines_it)
        
        graph_it  = relgraph.RelationGraph(lines_it, assign_it)
        
        precision, recall = test_relations.process(gtlines, gtassign, lines_it, assign_it, primitiveCorres, primitiveCorresId, [], gtgraph, graph_it)
        F = 0.
        if precision != 0 and recall != 0:
            F = 2.*(precision * recall) / (precision + recall)
        print it, sigma, F, precision, recall
    else:
        print it, sigma

