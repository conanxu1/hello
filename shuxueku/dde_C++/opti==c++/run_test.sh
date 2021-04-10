cd ../../..
git pull
<<<<<<< HEAD


=======
>>>>>>> 7c0cc035a2e82d50ab769b570c848bfbbaa25144
cd shu*/dde*/*++



rm -rf build_by_cmake
mkdir build_by_cmake
cd build_by_cmake

echo $1 
cmake .. $1 
make -j4
echo "======>>>>>"
chmod 777 ddetest
./ddetest
