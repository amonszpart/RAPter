import packages.project as project
import packages.primitive as primitive
import packages.io
import argparse
from matplotlib import pyplot as plt
import matplotlib.mlab as mlab
import numpy as np
from scipy.stats import norm,kstest,skewtest,kurtosistest,normaltest
from scipy.integrate import quad

################################################################################
## Visualization parameters
nbBins   = 200
rangeMin = -0.05
rangeMax =  0.05



################################################################################
## UI Generation
def setupHistogramUI(distribution, title, c1, c2, pdf=None):
    # Generate UI
    fig, ax1 = plt.subplots()
    fig.canvas.set_window_title(title)
    ax1.set_ylabel('PDF')

    if pdf != None:
        ax1.plot(x, pdf, color=c1, linestyle='--')
        
    # Generate histogram and display it
    hist, bins = np.histogram(distribution, nbBins/4, (rangeMin, rangeMax), normed=True)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2

    ax2 = ax1.twinx()
    ax2.set_ylabel('Noise')
    ax2.bar(center, hist, align='center', color=c1, width=width, edgecolor = "none")

    # Compute PDF of the gt histogram and display
    mean = np.mean(distribution)
    variance = np.var(distribution)
    sigma = np.sqrt(variance)
    mpdf = mlab.normpdf(x,mean,sigma)
    ax1.plot(x, mpdf, linewidth=3, color=c1)
    ax1.fill_between(x, mpdf, 0, color=c2 )
    

def setupComparisonUI(distrib1, title1, c1, 
                      distrib2, title2, c2,
                      kernel = None, c3 = None):
    fig, ax1 = plt.subplots()
    fig.canvas.set_window_title('CDFs')
    
    hist1, bins = np.histogram(distrib1, nbBins, (rangeMin, rangeMax), density=True)    
    hist2, bins = np.histogram(distrib2, nbBins, (rangeMin, rangeMax), density=True)
    
    ax1.plot(x, np.cumsum(hist1)/np.sum(hist1), color=c1)
    ax1.plot(x, np.cumsum(hist2)/np.sum(hist2), color=c2)
    
    if kernel != None:
      ax1.plot(x, kernel.cdf(x), color=c3, linestyle='--')
      
      print 'Analyzing ',title1,': ',kstest(distrib1, kernel.cdf, mode='asymp')
      print 'Analyzing ',title2,': ',kstest(distrib2, kernel.cdf, mode='asymp')
    
    
    



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



################################################################################
## Process data

# Compute the distance between each point and its assigned primitive
gtDistrib   = [ [primVar.distanceTo(cloud[a[0]]) for primVar in gtlines if primVar.uid == a[1]][0] for a in gtassign]
distrib_it1 = [ [primVar.distanceTo(cloud[a[0]]) for primVar in lines_it1 if primVar.uid == a[1]][0] for a in assign_it1]


################################################################################
## Generate UI
x = np.linspace(rangeMin, rangeMax, nbBins)

pdf = project.kernels[0]._pdf(x)

setupHistogramUI(gtDistrib, "Ground truth", '#dbb971', '#eddebc', pdf)
setupHistogramUI(distrib_it1, "Estimation", '#71b2db', '#b4d9f1', pdf)
setupComparisonUI(gtDistrib, "Ground truth", '#dbb971', 
                  distrib_it1, "Estimation", '#71b2db',
                  project.kernels[0], '#dbb971')
 
#print 'normal skewtest teststat = %6.3f pvalue = %6.4f' % skewtest(gtDistrib)
#print 'normal kurtosistest teststat = %6.3f pvalue = %6.4f' % kurtosistest(gtDistrib)
#print 'normal teststat = %6.3f pvalue = %6.4f' % normaltest(gtDistrib)
#skewtest(gtDistrib)                 
#print kstest(gtDistrib, project.kernels[0].cdf, mode='asymp')

plt.show()
