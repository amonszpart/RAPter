globOpt
=======

globOpt


##Dependencies
* 3rd Party:
** CoinBonmin (optimisation)
** libfbi (fast box intersection)
* SmartGeometryGroup:
** QCQPcpp
** geometry-tools
All dependencies are assumed to be available in a workspace folder and a `3rdparty` subdirectory.

GlobOpt also requires OpenCV and PointCloudLibrary, and we recommend to use respectively the last stable release and the development branch for these packages (tested with Ubuntu 14.04 and Debian Testing).

###CoinBonmin
First, download Bonmin in `${WORKSPACE}/3rdparty`:
```
svn co https://projects.coin-or.org/svn/Bonmin/stable/1.5 CoinBonmin-stable
```
Then:
* install the following packages (Debian/Ubuntu): `liblapack-dev libblas-dev fortran-compiler`.
* install 3rd party solver using the script `path/to/Bonmin/ThirdParty/Mumps/get.Mumps`.
* Compile and install:
```
mkdir build
cd build
../configure -C
make
make -install
```

###libfbi
Download in `${WORKSPACE}/3rdparty` from [Github](https://github.com/mkirchner/libfbi.git).
No installation required.


###QCQPcpp
Download QCQPcpp in `${WORKSPACE}` from [Github](https://github.com/amonszpart/QCQPcpp), then build:
```
mkdir build
cd build
cmake ..
make
```

You can run cmake with `-DUSE_MOSEK` if you need it.

###geometry-tools
Download geometry-tools in `${WORKSPACE}` from [Github](https://github.com/smartgeometry-ucl/geometry-tools), then build:
```
mkdir build
cd build
cmake ..
make
```

##Compilation
###GlobOpt
Once all the dependencies are satisfied, you can build GlobOpt. By default, dependencies are expected to be in `/home:${USER}/workspace/` (you may need to edit CMakelist.txt to change that behaviour). 
```
mkdir build
cd build
cmake ..
make
```

###InputGen
Just compile using cmake like usual.
