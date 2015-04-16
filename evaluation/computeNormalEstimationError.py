import argparse
from matplotlib import pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm,kstest,skewtest,kurtosistest,normaltest
import os.path
import unitTests.relations as test_relations
import packages.relationGraph as relgraph
import packages.color as color
import packages.polarPlot as polarplot
import math
import glob


def createGraph(anglesArrays, names):
    N = 20                                  # number of bins in the histogram
    mmin = 0.
    mmax = math.pi
    ind = np.arange(mmin, mmax, N)          # the x locations for the groups
                                            # we assume this to be uniform
    width = (mmax-mmin) / float(N*len(anglesArrays))   # width of each bar
    
    colors = color.getMediumPalette()
    
    
    fig, ax = plt.subplots()
    offset=0.
    for i, angle in enumerate(anglesArrays):
        n, bins = np.histogram(angle, N, (mmin, mmax), density=True)           
        ax.bar(offset+bins[:-1], n, width, color=colors[i])
        offset = offset+width
        
        print names[i]
        print "mean   = ", np.mean(angles)
        print "var    = ", np.var(angles)
        print "median = ", np.median(angles)

    plt.show()

## Load a mesh as GT (must have correct normals), the input point cloud,
## and the reconstruction output: primitives + assignment.
##
##   1. Compute the cloud to GT assignment
##      We now have two assigments per points: GT and estimated
##   2. For each point, compute a distance between GT and estimated normal 
##   3. Analyse the resulting distance distribution:
##        - variance and mean
##
##
##
################################################################################
## Command line parsing
parser = argparse.ArgumentParser(description='Compute reconstruction error in normal space.')
parser.add_argument('angles', help='Can be either a file listing angles or a directory containing multiple of these files.')

args = parser.parse_args()

anglesFile  = args.angles


if (os.path.isfile(anglesFile)) :
    angles = np.array( [np.float32(line.strip()) for line in open(anglesFile)] )
    createGraph( [angles],[anglesFile] )

else:
    filelist = glob.glob(anglesFile+"*.angles.csv")
    angles = []
    for f in filelist:
        angles.append( np.array( [np.float32(line.strip()) for line in open(f)] ) )

    createGraph( angles,filelist )


