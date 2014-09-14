import packages.project as project
import argparse

parser = argparse.ArgumentParser(description='Compare ground truth noise distribution (continuous generator and generated samples) and the result of the optimisation.')
parser.add_argument('projectfile', type=argparse.FileType('r'))

args = parser.parse_args()

project = project.PyProject(args.projectfile)
