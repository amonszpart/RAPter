#!/bin/bash

~/workspace/pcd2ply/build$ for i in ls /home/bontius/workspace/globOpt/data/scenes/fall/fall15__noisy*.pcd; do `./pcd2ply $i $i.ply`; done;

