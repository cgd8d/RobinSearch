# Install GMP from source.

wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar -xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --help
CC=clang-18 CFLAGS="-march=native -O2" AR=llvm-ar-18 NM=llvm-nm-18 RANLIB=llvm-ranlib-18 ./configure --disable-shared --disable-fft
make -j --output-sync
make -j --output-sync check
#cd tune
#make speed
ls -a
