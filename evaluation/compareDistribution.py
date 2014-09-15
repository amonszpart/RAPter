import packages.project as project
import packages.primitive as primitive
import argparse
from matplotlib import pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Compare ground truth noise distribution (continuous generator and generated samples) and the result of the optimisation.')
parser.add_argument('projectdir')

args = parser.parse_args()

projectdir = args.projectdir
if projectdir[-1] == '/':
    projectdir = projectdir[:-1]

projectname = projectdir.split('/')[-1]
projectfile = projectdir+'/gt/'+projectname+'.prj'
gtlinesfile = projectdir+'/gt/primitives.csv'

print 'Processing project ', projectname

project = project.PyProject(projectfile)
gtlines = primitive.readPrimitivesFromFile(gtlinesfile)

x = np.linspace(-0.2, 0.2, 200)

for k in project.kernels:
    plt.plot(x, k.getPDF(x))


plt.show()
