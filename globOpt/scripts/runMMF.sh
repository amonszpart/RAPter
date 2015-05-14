sudo modprobe nvidia-uvm
export PYTHONPATH=$PYTHONPATH:/home/bontius/workspace/3rdparty/opencv-trunk/install/lib/python2.7/dist-packages
python mmfExtract.py -i ./data/MIT_hallway_1_rgb.png --buffer --nMerge 39 --nIter 40 --plotInput --plotEval --cuda
