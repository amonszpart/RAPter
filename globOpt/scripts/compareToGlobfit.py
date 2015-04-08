#!/usr/bin/python
import argparse
import os

toGlobFit = "/home/bontius/workspace/globOpt/globOpt/build/Release/bin/toGlobFit";
#toGlobFit = "/export/home/kandinsky/nmellado/workspace_globOpt/globOpt/globOpt/build/Release/bin/toGlobFit"

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


parser = argparse.ArgumentParser()
parser.add_argument("-s", "--scale", type=float, help="Scale (rho) parameter as smallest feature size to preserve [0.001..0.05]")
parser.add_argument("-r", "--subsample", type=float, default=0.5, help="Decides how much to keep from the original assignments [0.01..1.0]")
parser.add_argument("--popLimit", type=int, default=20, help="Decides at minimum, how many points to keep for each plane. [ 20..100 ]");
parser.add_argument("--primLimit", type=int, default=0, help="Decides how many primitives to keep. 0 means keep all. [ 0..n ]");
parser.add_argument("--primRand", type=float, default=0., help="Decides how many of the kept primitives should be randomly chosen. [ 0.0..1.0 ]");
parser.add_argument("-p", "--primitives", default="patches.csv", help="Primitives.csv to convert to globfit input [segments.csv]")
parser.add_argument("-a", "--assoc", default="points_patches.csv", help="Path to point to plane assignments [points_segments.csv]")
parser.add_argument("--angleThresh", type=float, default="10.0", help="Angle threshold given to globfit with -o and -g [0.1..10.0]")
parser.add_argument('--run', action="store_true", help="Call the scripts (disabled by default)")

args = parser.parse_args()
print("\n");

runIsNotDrySoGoAhead=args.run

if not args.scale:
    print("Need scale -s, --scale")
    exit

print( "Setting --subsample to %f (keep 50 percent of the assignments)" % (args.subsample) )
print( "Setting --popLimit to %d, will keep at least %d points assigned to each plane" % (args.popLimit, args.popLimit) );
print( "Setting --assoc to %s" % (args.assoc) );
print( "Setting --primitives to %s" % (args.primitives) );

#cmd = "../runSegmentation.py -s %f --pl 4" % (args.scale)
#call(cmd)

cmd = "%s --subsample-primitives %f --pop-limit %d --prim-limit %d --prim-random-ratio %f --prims %s --cloud cloud.ply -a %s --scale %f" % (toGlobFit, args.subsample, args.popLimit, args.primLimit, args.primRand,  args.primitives, args.assoc, args.scale)
call(cmd, dry=not runIsNotDrySoGoAhead)

root_prim_str  = args.primitives[:-4]
root_assoc_str = args.assoc[:-4]
sub_str        = ".sub_" + str(args.subsample) + "_" + str(args.primLimit) + ".csv"
subPrimitive   = root_prim_str  + sub_str
subAssoc       = root_assoc_str + sub_str

cmd = "../runGlobfit.py --angle-thresh %f -s %f -p %s -a %s" % (args.angleThresh, args.scale, subPrimitive, subAssoc)
call(cmd, dry=not runIsNotDrySoGoAhead)

cmd = "../show.py -s %f -p %s -a %s&" % (args.scale, subPrimitive, subAssoc)
call(cmd, dry=not runIsNotDrySoGoAhead)

outnames = [ "oa", "pa", "ae" ]
for n in outnames:
    #name = root_prim_str + "_" + n + ".globfit" 
    name = "segments_" + n + ".globfit" 
    if os.path.isfile(name): 
        out = root_prim_str  + ".primitives."+n
        print "# Converting " + name + " to " + out + ".globfit.csv"
        cmd = "%s --from %s --planes --prims %s --cloud cloud.ply -a %s --scale %f -o %s" % (toGlobFit, name, args.primitives, args.assoc, args.scale, out)
        call(cmd, dry=not runIsNotDrySoGoAhead)
        
        cmd = "../show.py -s %f -p %s -a %s &" % (args.scale, out + ".globfit.csv", "points_" + out + ".globfit.csv")
        call(cmd, dry=not runIsNotDrySoGoAhead)



#../runGlobfit.py --save-pa segments_pa

#../show.py -p segments_pa.globfit.csv -a points_segments_pa.globfit.csv -s 0.05

print("\n");





# ../compareToGlobfit.py -s 0.004 --primLimit 250

# ../toGlobFit --subsample-primitives 0.500000 --pop-limit 20 --prim-limit 300 --prims patches.csv --cloud cloud.ply -a points_primitives.csv --scale 0.0025
# ../runGlobfit.py -s 0.004 -p patches.sub_0.1_250.csv -a points_patches.sub_0.1_250.csv
# ../show.py -s 0.004 -p patches.sub_0.1_250.csv -a points_patches.sub_0.1_250.csv
# ../toGlobFit --from segments_ea.globfit --planes --prims patches.sub_0.1_250.csv --cloud cloud.ply -a points_patches.sub_0.1_250.csv --scale 0.004000 -o primitives_ea
# ../show.py -s 0.004 -p primitives_ea.globfit.csv -a points_primitives_ea.globfit.csv --save-poly
# Equivalent: ../globOptVis --show3D --scale 0.004000 --pop-limit 3 --title primitives_ea.globfit.csv --angle-gens 90 --paral-colours 0.000100 --bg-colour .9,.9,.9 -p primitives_ea.globfit.csv -a points_primitives_ea.globfit.csv --use-tags --no-clusters --no-pop --statuses -1,1 --no-rel --dir-colours --no-scale --save-poly --draw-mode 21
# python normal_distr.py /home/bontius/workspace/globOpt/data/scenes-paper/bu_nolaSigg/cloudRGBNormal_noUnass_noPrim.ply out.svg "toto"


# ../toGlobFit --subsample-primitives 0.500000 --pop-limit 20 --prim-limit 300 --prims patches.csv --cloud cloud.ply -a points_primitives.csv --scale 0.0025
# ../runGlobfit.py -s 0.0025 -p patches.sub_0.5_300.csv -a points_patches.sub_0.5_300.csv
# 
