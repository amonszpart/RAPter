#!/usr/bin/python

from optparse import OptionParser
import os
# os.system()

parser = OptionParser()
parser.add_option("-s", "--scale", type="float", dest="scale",
                  help="scale (rho) parameter as smallest feature size to preserve [0.001..0.05]")
parser.add_option("--al", "--angle-limit", type="float", dest="angleLimit", default=0.1,
                  help="angleLimit (tau) in radians as mismatch between average normal of region and normal of new point to include [0.1..0.4]")
parser.add_option("--dl", "--distance-limit", type="float", dest="distanceLimitMultiplier", default=1.0,
                  help="Region growing merge distance is the result of scale*distLimitMult. Set it <1.0 for conservative region growing, and >1.0 for aggressive region growing.")
(options, args) = parser.parse_args()

print("options:", options)
if not options.scale:   # if filename is not given
    parser.error('-s, --scale is required');

cmd = "../glob_opt --segment3D --patch-pop-limit 0 --angle-limit %f --scale %f --dist-limit-mult %f --angle-gens 0" % (
    options.angleLimit, options.scale, options.distanceLimitMultiplier)
print("[Calling] ", cmd)
print(os.system(cmd))
os.rename("patches.csv", "segments.csv")
os.rename("points_primitives.csv", "points_segments.csv")
print("Output in segments.csv/points_segments.csv")
