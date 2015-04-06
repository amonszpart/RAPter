#!/usr/bin/python
from optparse import OptionParser

toGlobFit = "/home/bontius/workspace/globOpt/globOpt/build/Release/bin/toGlobFit";

def call(cmd, dry=True, noExit=False):

    print("%s" % (cmd))
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
                  help="Scale (rho) parameter as smallest feature size to preserve [0.001..0.05]")
parser.add_option("-r", "--subsample", type="float", dest="subsampleRatio", default=0.5,
                  help="Decides how much to keep from the original assignments [0.01..1.0]")
parser.add_option("--pl", "--popLimit", type="int", dest="popLimit", default=20,
                  help="Decides at minimum, how many points to keep for each plane. [ 20..100 ]");
parser.add_option("--prl", "--primLimit", type="int", dest="primLimit", default=0,
                  help="Decides how many primitives to keep. 0 means keep all. [ 0..n ]");
parser.add_option("-p", "--primitives", type="string", dest="primitivesPath", default="segments.csv",
                  help="Primitives.csv to convert to globfit input [segments.csv]")
parser.add_option("-a", "--assoc", type="string", dest="assocPath", default="points_segments.csv",
                  help="Path to point to plane assignments [points_segments.csv]")

(options, args) = parser.parse_args()
print("\n");

if not options.scale:
    print("Need scale -s, --scale")
    exit

print( "Setting --subsample to %f (keep 50 percent of the assignments)" % (options.subsampleRatio) )
print( "Setting --popLimit to %d, will keep at least %d points assigned to each plane" % (options.popLimit, options.popLimit) );
print( "Setting --assoc to %s" % (options.assocPath) );
print( "Setting --primitives to %s" % (options.primitivesPath) );

(options, args) = parser.parse_args()

cmd = "../runSegmentation.py -s %f --pl 4" % (options.scale)
call(cmd)

cmd = "%s --subsample-primitives %f --pop-limit %d --prim-limit %d --prims %s --cloud cloud.ply -a %s --scale %f" % (toGlobFit, options.subsampleRatio, options.popLimit, options.primLimit, options.primitivesPath, options.assocPath, options.scale)
call(cmd)

cmd = "../runGlobfit.py -s %f -p %s -a %s" % (options.scale, options.primitivesPath, options.assocPath)
call(cmd)

#../runGlobfit.py --save-pa segments_pa

#../show.py -p segments_pa.globfit.csv -a points_segments_pa.globfit.csv -s 0.05

print("\n");
