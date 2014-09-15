import packages.project as project
import packages.primitive as primitive
import packages.io
import argparse
from matplotlib import pyplot as plt
import matplotlib.mlab as mlab
import numpy as np

################################################################################
## Visualization parameters
nbBins   = 200
rangeMin = -0.2
rangeMax =  0.2



################################################################################
## UI Generation
def setupUI(primArray, assignArray, title):
    # Generate UI
    fig, ax1 = plt.subplots()
    fig.canvas.set_window_title(title)
    ax1.set_ylabel('PDF')

    for k in project.kernels:
        pdf = k.getPDF(x)
        ax1.plot(x, pdf, color='#dbb971', linestyle='--')
        

    # Compute the distance between each point and its assigned primitive
    d = [ [primVar.distanceTo(cloud[a[0]]) for primVar in primArray if primVar.uid == a[1]][0] for a in assignArray]

    # Generate histogram and display it
    hist, bins = np.histogram(d, nbBins/5, (rangeMin, rangeMax))
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    ax2 = ax1.twinx()
    ax2.set_ylabel('Noise')
    ax2.bar(center, hist, align='center', color='#dbb971', width=width, edgecolor = "none")

    # Compute PDF of the gt histogram and display
    mean = np.mean(d)
    variance = np.var(d)
    sigma = np.sqrt(variance)
    mpdf = mlab.normpdf(x,mean,sigma)
    ax1.plot(x, mpdf, linewidth=3, color='#dbb971')
    ax1.fill_between(x, mpdf, 0, color='#eddebc' )




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

print 'Processing project ', projectname




################################################################################
## Reading input files
project  = project.PyProject(projectfile)
cloud    = packages.io.readPointCloudFromPly(cloudfile)
gtlines  = primitive.readPrimitivesFromFile(gtlinesfile)
gtassign = packages.io.readPointAssignementFromFiles(gtassignfile)



################################################################################
## Generate UI
x = np.linspace(rangeMin, rangeMax, nbBins)

setupUI(gtlines, gtassign, "Ground truth")
setupUI(gtlines, gtassign, "Estimation")

plt.show()
