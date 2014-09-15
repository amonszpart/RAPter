import packages.project as project
import packages.primitive as primitive
import packages.io
import argparse
from matplotlib import pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Compare ground truth noise distribution (continuous generator and generated samples) and the result of the optimisation.')
parser.add_argument('projectdir')

args = parser.parse_args()

nbBins   = 200
rangeMin = -0.2
rangeMax =  0.2

projectdir = args.projectdir
if projectdir[-1] == '/':
    projectdir = projectdir[:-1]

projectname  = projectdir.split('/')[-1]
projectfile  = projectdir+'/gt/'+projectname+'.prj'
gtlinesfile  = projectdir+'/gt/primitives.csv'
gtassignfile = projectdir+'/gt/points_primitives.csv'
cloudfile    = projectdir+'/cloud.ply'

print 'Processing project ', projectname

project  = project.PyProject(projectfile)
cloud    = packages.io.readPointCloudFromPly(cloudfile)
gtlines  = primitive.readPrimitivesFromFile(gtlinesfile)
gtassign = packages.io.readPointAssignementFromFiles(gtassignfile)

# Generate UI
fig, ax1 = plt.subplots()
ax1.set_ylabel('Distribution', color='g')


# Generate noise kernel PDF curves
x = np.linspace(rangeMin, rangeMax, nbBins)

for k in project.kernels:
    pdf = k.getPDF(x)
    ax1.plot(x, pdf, color='g')
    ax1.fill_between(x, pdf, 0, color='g' )


# Compute the distance between each point and its assigned primitive
d = [ [x.distanceTo(cloud[a[0]]) for x in gtlines if x.uid == a[1]][0] for a in gtassign]

# Generate histogram and display it
hist, bins = np.histogram(d, nbBins/2, (rangeMin, rangeMax))
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

ax2 = ax1.twinx()
ax2.bar(center, hist, align='center', width=width)





#plt.plot(hist[1], hist[0])


plt.show()
