#cd ../../..
#git pull git://github.com/conanxu1/hello.git
#cd shu*/dde*/*++



rm -rf build_by_cmake
mkdir build_by_cmake
cd build_by_cmake

echo $1 
cmake .. $1 
make -j4
echo "======>>>>>"
chmod 777 ddetest
./ddetest
