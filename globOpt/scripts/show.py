#!/usr/bin/python

from optparse import OptionParser
import os
# os.system()

parser = OptionParser()
parser.add_option("-s", "--scale", type="float", dest="scale",
                  help="scale (rho) parameter as point inlier distance for planes [0.001..0.05]")
parser.add_option("-p", "--primitives", type="string",
                  dest="primitivesPath", help="path to primitives file to display [primitives.csv]")
parser.add_option("-a", "--associations", type="string",
                  dest="assignmentsPath", help="path to points_primitives file containing point to primitive assignments [points_primitives.csv]")
parser.add_option("--pl", "--pop-limit", type="int", dest="popLimit", default=3,
                  help="Primitives having less points assigned, than this should not be rendered. [0..20]")
parser.add_option(
    "", "--title", type="string", dest="title", help="Window title")
parser.add_option("", "--angle-gens", type="string", dest="angleGens", default="90",
                  help="Which angles count as perfect angles. Only relevant, if relations are shown.")
parser.add_option("--2d", "--2D", action="store_true", dest="_2d", default=False, help="Interpret primitives as lines")
parser.add_option("--bg-colour", "", default=".9,.9,.9", dest="bgColour",
                  help="Background colour as r,g,b comma separated, no spaces")
parser.add_option("--show-scale", "", action="store_true", dest="showScale",
                  default=False, help="Show scale circle at origin")
parser.add_option("--paral-colours", "", type="float", dest="paralLimit", default=0.0001,
                  help="Show colours based on parallel normals. The threshold is this number in radians." )
parser.add_option("--dir-colours", "", action="store_true", dest="dirColours", default=False,
                  help="Colour based on cluster id and not parallelity" )

(options, args) = parser.parse_args()

if not options.scale:
    parser.error('-s, --scale is required')
if not options.primitivesPath:   # if filename is not given
    parser.error('-p, --primitives is required')
if not options.assignmentsPath:   # if filename is not given
    parser.error('-a, --assoc is required')
if not options.title:
    options.title = options.primitivesPath

defaultArgs = "--use-tags --no-clusters --no-pop --statuses -1,1 --no-rel --dir-colours"
if not options.showScale:
    defaultArgs += " --no-scale";
    
showOption = "--show" if options._2d else "--show3D"
print("default:", defaultArgs)
colourOption = "--dir-colours" if options.dirColours else "--paral-colours %f" % (options.paralLimit);

cmd = "../globOptVis %s --scale %f --pop-limit %d --title %s --angle-gens %s %s --bg-colour %s -p %s -a %s %s" % (
    showOption, options.scale, options.popLimit, options.title, options.angleGens, colourOption, options.bgColour, options.primitivesPath, options.assignmentsPath, defaultArgs)

print("[Calling] ", cmd)
print(os.system(cmd))
