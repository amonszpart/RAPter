globOpt
=======

globOpt


##Dependencies

###CoinBonmin
First, download Bonmin:
```
svn co https://projects.coin-or.org/svn/Bonmin/stable/1.5 CoinBonmin-stable
```
Then:
* install the following packages (Debian/Ubuntu): 'liblapack-dev libblas-dev fortran-compiler'.
* install 3rd party solver using the script 'path/to/Bonmin/ThirdParty/Mumps/get.Mumps'.
* Compile and install:
```
mkdir build
cd build
../configure -C
make
make -install
```


###QCQPcpp
Download and install from [Github](https://github.com/amonszpart/QCQPcpp)
```
mkdir build
cd build
cmake ..
make
```

You can run cmake with '-DUSE_MOSEK' if you need it.


##Compilation
###GlobOpt
You need to have all the dependencies in the folder '~/workspace/3rdparty/'
```
mkdir build
cd build
cmake ..
make
```

###InputGen