ssh amonszpa@geometry.cs.ucl.ac.uk -L 2201:128.16.11.154:22
ssh globopt@localhost -p 2201

#python ../generatePolarPrimDirections.py segments.csv  --logscale --shownormals --useInactivePrims

# run schnabel
../ransac --schnabel3D --cloud cloud.ply -p patches.csv -a points_primitives.csv --scale 0.025 --minsup 600

#run pearl
../pearl --3D --scale 0.02 --pw 200 --cmp 250 -N 300

# show pearl
../globOptVis --show --scale 0.01 -p primitives.pearl.csv -a points_primitives.pearl.csv --pop-limit 3 --title "Ransac output" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel

# show ransac
../globOptVis --show --scale 0.01 -p primitives.ransac.csv -a points_primitives.ransac.csv --pop-limit 3 --title "Ransac output" --angle-gens 90 --use-tags --no-clusters --statuses -1,1 --no-pop --dir-colours --no-scale --bg-colour .9,.9,.9 --no-rel
