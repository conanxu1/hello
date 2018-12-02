


# aria2c -c -s 16 -x 16 https://github.com/downloads/chokkan/liblbfgs/liblbfgs-1.10.tar.gz 
# 7z x liblbfgs-1.10.tar.gz 
7z x liblbfgs-1.10.tar
cd liblbfgs-1.10



./configure

sudo make

sudo make install
