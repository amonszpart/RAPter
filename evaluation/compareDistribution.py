import packages.project as project
import argparse
from matplotlib import pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description='Compare ground truth noise distribution (continuous generator and generated samples) and the result of the optimisation.')
parser.add_argument('projectfile', type=argparse.FileType('r'))

args = parser.parse_args()

project = project.PyProject(args.projectfile)

x = np.linspace(-0.2, 0.2, 200)

for k in project.kernels:
    plt.plot(x, k.getPDF(x))


plt.show()
