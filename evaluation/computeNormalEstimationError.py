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
import packages.polarPlot as polarplot
import math

compareToGt = False

itmin=0
itmax=20


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
parser.add_argument('gtmesh')
parser.add_argument('cloud')
parser.add_argument('primitives')
parser.add_argument('assignment')

args = parser.parse_args()


