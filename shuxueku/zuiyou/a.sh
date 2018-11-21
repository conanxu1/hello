
cd /root

#git clone git://github.com/xianyi/OpenBLAS
apt install -y aria2 p7zip-full make gcc cmake


rm OpenBLAS-0.3.3*

aria2c -c -x 16 -s 16 https://github.com/xianyi/OpenBLAS/archive/v0.3.3.zip

7z x OpenBLAS-0.3.3.zip

cd 
apt install gfortran

cd OpenBLAS-0.3.3
sudo make FC=gfortran -j4

make PREFIX=/root/Openblas  install

rm /usr/lib/libblas.so.3
rm /usr/lib/liblapack.so.3
rm /usr/lib/libopenblas.so.0



ln -s /root/Openblas/lib/libopenblas.so  /usr/lib/libblas.so.3

ln -s /root/Openblas/lib/liblapack.so.3 /usr/lib/liblapack.so.3

ln -s /root/Openblas/lib/libopenblas.so.0 /usr/lib/libopenblas.so.0


#git submodule update –init –recursive
#make -j 8 -f Makefile.install install
