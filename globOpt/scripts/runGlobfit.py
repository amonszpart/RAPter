#!/usr/bin/python

import os  # os.system()
from optparse import OptionParser


def call(cmd, dry=False, noExit=False):
    print("[Calling] %s" % (cmd) )
    if not dry:
        ret = os.system(cmd)
        if ret != 0:
            print("call returned error ", ret, ", aborting")
            if not noExit:
                exit
            else:
                print("call returned ", ret)


parser = OptionParser()
parser.add_option("-s", "--scale", type="float", dest="scale",
                  help="scale (rho) parameter as smallest feature size to preserve [0.001..0.05]")
parser.add_option("--matlab", "", type="string", dest="matlabExecutableFolder",
                  default="/home/bontius/matlab_symlinks/", help="Constains the executable \"matlab\"")
parser.add_option("--toGlobFit", "", type="string", dest="toGlobFitExecutable",
                  default="/home/bontius/workspace/globOpt/globOpt/build/Release/bin/toGlobFit", help="Path to executable \"toGlobFit\"")
parser.add_option("--globfit", "", type="string", dest="globfitExecutable",
                  default="/home/bontius/workspace/3rdparty/globfit/build/bin/globfit_release", help="Path to executable \"globfit\"")

parser.add_option("-p", "--primitives", type="string", dest="primitivesPath", default="segments.csv",
                  help="Primitives.csv to convert to globfit input [segments.csv]")
parser.add_option("-a", "--assoc", type="string", dest="assocPath", default="points_segments.csv",
                  help="Path to point to plane assignments [points_segments.csv]")
parser.add_option("-c", "--cloud", type="string", dest="cloudPath",
                  default="cloud.ply", help="Path to pointcloud ply. [cloud.ply]")
parser.add_option("-g", "--gf-output-path", type="string", dest="globfitOutputPath", default="segments_oa.globfit",
                  help="Path to globfit output to convert back to csv. [segments_oa.globfit, segments_pa.globfit]")
parser.add_option("-n", "--dry", action="store_true", default=False, help="Just print commands, don't run")

parser.add_option("--save-pa", "", type="string", dest="savePa", help="Don't do anything else, just copy the segments_pa.globfit file to segments_pa.csv" );


(options, args) = parser.parse_args()

if not options.primitivesPath:   # if filename is not given
    parser.error('-p, --primitives is required')
if not options.assocPath:   # if filename is not given
    parser.error('-a, --assoc is required')

if options.savePa:
    cmd = "%s --from %s --planes --prims %s --cloud %s -a %s --scale %f -o %s" % (options.toGlobFitExecutable, "segments_pa.globfit", options.primitivesPath, options.cloudPath, options.assocPath, 0.0, options.savePa)
    call(cmd, dry=options.dry)
    exit;

if not options.scale:   # if filename is not given
    parser.error('-s, --scale is required')

# Modify path
cmd = "export PATH=\"%s\":$PATH" % (options.matlabExecutableFolder)
# call(cmd);
print ("path: ", os.environ["PATH"])
os.environ["PATH"] = "%s:%s" % (options.matlabExecutableFolder, os.environ["PATH"])
print ("path: ", os.environ["PATH"])
exit

cmd = "%s --planes --prims %s --cloud %s -a %s --scale %f" % (
    options.toGlobFitExecutable, options.primitivesPath, options.cloudPath, options.assocPath, options.scale)
call(cmd, dry=options.dry)

cmd = "%s -i segments.globfit -v -o 3.0 -g 3.0 -a 0.1 -p %f -l %f -r %f" % (
    options.globfitExecutable, options.scale, options.scale, options.scale)
call(cmd, dry=options.dry)

cmd = "%s --from %s --planes --prims %s --cloud %s -a %s --scale %f" % (
    options.toGlobFitExecutable, options.globfitOutputPath, options.primitivesPath, options.cloudPath, options.assocPath, options.scale)
call(cmd, dry=options.dry)

# ../../globOptVis --show3D --scale $scale --pop-limit $poplimit -p primitives_it9.bonmin.csv -a points_primitives_it8.csv --title "GlobOpt" $visdefparam --paral-colours --no-rel &

# ../../globOptVis --show3D --scale $scale --pop-limit $poplimit -p primitives.globfit.csv -a points_primitives.globfit.csv --title "Globfit" $visdefparam --paral-colours --no-rel &

# ../../globOptVis --show3D --scale $scale --pop-limit $poplimit -p patches.csv -a points_primitives.csv --title "Input" $visdefparam --paral-colours --no-rel &
