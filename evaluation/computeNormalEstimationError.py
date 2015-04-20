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

# Surrounds output with \textbf{} if minimum in row (or maximum in case of "precision" ).
def checkMin( value, values, unitStr = "", doMax = False ):
    # this is really stupid, I know...
    extrema = True;
    for v in values:
        if doMax:
            if values[v] > value:
                extrema = False;
                break;
        else:
            if values[v] < value:
                extrema = False;
                break;
    
    if not extrema:
        pattern =  "%%s & %%2.4f%s" % unitStr;
    else:
        pattern = "%%s & \\textbf{%%2.4f%s}" % unitStr;

    return pattern;

# Hack to get table column title (method name) from input angle filename
def guessShortNames( names ):
    shortNames = [];
    order = [];
    for name in names:
        if name.lower().find("schnabel") >= 0:
            shortNames.append( "RANSAC" );
            order.append(0);
        elif name.lower().find("globfit") >= 0 or name.lower().find("sub") >= 0:
            shortNames.append("Globfit");
            order.append(2);
        elif name.lower().find("pearl") >= 0:
            shortNames.append("\\textsc{Pearl}");
            order.append(1);
        elif name.lower().find("bonmin") >= 0:
            shortNames.append("Ours");
            order.append(3);
        else:
            shortNames.append(name);
            order.append(2);

    return shortNames,order;

def createGraph(anglesArrays, names):
    N = 20                                  # number of bins in the histogram
    mmin = 0.
    mmax = math.pi
    ind = np.arange(mmin, mmax, N)          # the x locations for the groups
                                            # we assume this to be uniform
    width = (mmax-mmin) / float(N*len(anglesArrays))   # width of each bar
    
    colors = color.getMediumPalette()
    
    # identify method names ("RANSAC, Pearl, ...") and their order in the table 
    shortNames,order = guessShortNames(names);
    # init storage dictionaries. key: shortName, value: statistic (i.e. mean)
    means       = {};
    variances   = {};
    medians     = {};
    counts      = {};
    precisions  = {};

    fig, ax = plt.subplots()
    offset=0.
    for i, angle in enumerate(anglesArrays):
        n, bins = np.histogram(angle, N, (mmin, mmax))
        ax.bar(offset+bins[:-1], n, width, color=colors[i])
        offset = offset+width
        
        # estimate
        mean      = np.mean(angle);
        variance  = np.var(angle);
        median    = np.median(angle);
        precision = len(angle[np.where( angle == 0. )]) / float(len(angle));

        # print
        print "mean   = ", mean, "rad,", mean * 180.0 / math.pi, "deg"
        print "var    = ", variance, "rad, ", variance * 180.0 / math.pi, "deg"
        print "median = ", median, "rad, ", median * 180.0 / math.pi, "deg"
        print "precision = ", precision * 100.0, "%";

        # record
        means[ shortNames[i] ] = mean;
        variances[ shortNames[i] ] = variance;
        medians[ shortNames[i] ]   = median;
        counts[ shortNames[i] ] = len(angle);
        precisions[ shortNames[i] ] = precision;

    # init latex strings
    strTitles       = "    &";
    strMeans        = "    & mean";
    strVars         = "    & variance";
    strMedians      = "    & median";
    strPrecisions   = "    & ==$0^{\\circ}$";
    strCounts       = "    & N ";

    # init ordered names
    keys = [None] * len(order);
    # determine final position
    for i,pos in enumerate(order):
        # hack for missing methods
        finalPos = sorted(order).index(pos);
        # store keys in order
        keys[finalPos] = shortNames[i];

    # print values in order:
    for key in keys:
        # print method name, just to be sure
        strTitles = "%s & %s" % (strTitles, key);
        
        # Check, if this method had the minimum score.
        # If yes, surround with \textbf{}. Appeend degree sign in either case.
        pattern = checkMin( means[key], means, "$^{\\circ}$" );
        strMeans = pattern % (strMeans, means[key] * 180.0 / math.pi );

        # Check, if this method had the minimum score.
        # If yes, surround with \textbf{}. Appeend degree sign in either case.
        pattern = checkMin( variances[key], variances, "$^{\\circ}$"  );
        strVars = pattern % (strVars,   variances[key] * 180.0 / math.pi );
        
        # Check, if this method had the minimum score.
        # If yes, surround with \textbf{}. Appeend degree sign in either case.
        pattern = checkMin( medians[key], medians, "$^{\\circ}$" );
        strMedians = pattern % (strMedians,  medians[key] * 180.0 / math.pi );

        # Check, if this method had the *maximum* score.
        # If yes, surround with \textbf{}. Appeend degree sign in either case.
        pattern = checkMin( precisions[key], precisions, "\\%%", True );
        strPrecisions = pattern % (strPrecisions, precisions[key] * 180.0 / math.pi );

        # Print number of samples (pairs of points compared)
        strCounts = "%s & %d" % (strCounts, counts[key] );

    # Print to console. 
    # TODO: print to file with scene name
    print ( "%s %s" % (strTitles , " \\\\ ") );
    print ( "%s %s" % (strMeans  , " \\\\ ") );
    print ( "%s %s" % (strVars   , " \\\\ ") );
    print ( "%s %s" % (strMedians, " \\\\ ") );
    print ( "%s %s" % (strPrecisions, " \\\\ ") );
    print ( "%s %s" % (strCounts, " \\\\ ") );

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


