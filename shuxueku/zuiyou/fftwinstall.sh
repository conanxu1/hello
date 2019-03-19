./configure --prefix=/usr    \
            --enable-shared  \
            --enable-threads \
           --enable-sse2    \
           --enable-avx     &&
sudo make
sudo make install

sudo make clean &&

./configure --prefix=/usr    \
           --enable-shared  \
          --enable-threads \
          --enable-sse2    \
          --enable-avx     \
          --enable-float   &&
sudo make

sudo make install

sudo make clean &&
./configure --prefix=/usr    \
           --enable-shared  \
            --enable-threads \
          --enable-long-double &&
sudo make


sudo make install
