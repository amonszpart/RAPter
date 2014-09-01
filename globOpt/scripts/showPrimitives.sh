#!/bin/bash
../Release/bin/gurobi_opt --show --dir . --cloud cloud.ply --scale 0.05f --assoc points_primitives.txt --use-tags --no-clusters --prims $1 
