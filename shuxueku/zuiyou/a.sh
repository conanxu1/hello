cd /root

#git clone git://github.com/xianyi/OpenBLAS
#apt install aria2 p7zip-full

#aria2c -c -x 16 -s 16 https://github.com/xianyi/OpenBLAS/archive/v0.3.3.zip

#7z x v0.3.3.zip


#apt install gfortran

cd OpenBLAS-0.3.3
#sudo make FC=gfortran

#make PREFIX=/root/Openblas  install

ln -s /root/Openblas/lib/libopenblas.so  /usr/lib/libblas.so.3

ln -s /root/Openblas/lib/liblapack.so.3 /usr/lib/liblapack.so.3

ln -s /root/Openblas/lib/libopenblas.so.0 /usr/lib/libopenblas.so.0 


#git submodule update –init –recursive 
#make -j 8 -f Makefile.install install
