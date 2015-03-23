grep -rlw libs/thread -e "TIME_UTC" | xargs sed -i 's/TIME_UTC/TIME_UTC_/g'
./bootstrap.sh --prefix=/home/bontius/workspace/3rdparty/boost_1_49_0/install cxxflags=-fPIC
./b2 --build-dir=/home/bontius/workspace/3rdparty/boost_1_49_0/build --with-thread install
./b2 link=static -j 4 --build-dir=/home/bontius/workspace/3rdparty/boost_1_49_0/build runtime-link=static --with-thread --with-system --with-program_options --with-filesystem -d+2 -a install

sudo apt-get install csh