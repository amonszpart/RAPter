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


compareToGt = False

itmin=0
itmax=3



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

sigma_ref = 1.

if compareToGt:
    projectfile  = projectdir+'/gt/'+projectname+'.prj'
    gtlinesfile  = projectdir+'/gt/primitives.csv'
    gtassignfile = projectdir+'/gt/points_primitives.csv'
    
    project  = project.PyProject(projectfile)
    gtlines  = primitive.readPrimitivesFromFile(gtlinesfile)
    gtassign = packages.io.readPointAssignementFromFiles(gtassignfile)
    
    gtlines  = packages.processing.removeUnassignedPrimitives(gtlines, gtassign)
    gtassign = packages.processing.removeUnassignedPoint(gtlines, gtassign)
    
    ############################################################################
    ## Process noise
    sigma_ref = project.kernels[0].stdev()
    gtDistrib   = utils.distanceToPrimitives(cloud, gtassign,   gtlines)
    print "Ground truth sigma ", np.sqrt(np.var(gtDistrib)) / sigma_ref    
    



for it in range(itmin, itmax):
    print "Processing iteration",it

    linesfile_it  = projectdir+'/primitives_merged_it'+str(it)+'.csv'
    assignfile_it = projectdir+'/points_primitives_it'+str(it)+'.csv'


    ################################################################################
    ## Reading input files
    lines_it  = primitive.readPrimitivesFromFile(linesfile_it)
    assign_it = packages.io.readPointAssignementFromFiles(assignfile_it)
    
    lines_it  = packages.processing.removeUnassignedPrimitives(lines_it, assign_it)
    assign_it = packages.processing.removeUnassignedPoint(lines_it, assign_it)


    ################################################################################
    ## Process noise

    # Compute the distance between each point and its assigned primitive
    distrib_it = utils.distanceToPrimitives(cloud, assign_it, lines_it)
    print "Estimated sigma ", np.sqrt(np.var(distrib_it)) / sigma_ref


