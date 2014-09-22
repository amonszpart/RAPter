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

import unitTests.relations as test_relations


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
## Process
print test_relations.process(gtlines, gtassign, lines_it1, assign_it1, primitiveCorres, primitiveCorresId, True)


