# Install GMP from source.

wget https://ftp.gnu.org/gnu/gmp/gmp-6.2.1.tar.xz
tar -xf gmp-6.2.1.tar.xz
cd gmp-6.2.1
./configure --help
CC=clang-12 CFLAGS="-fuse-ld=gold -flto=thin -march=native -O3" AR=llvm-ar-12 NM=llvm-nm-12 RANLIB=llvm-ranlib-12 ./configure --disable-shared --disable-fft
make -j --output-sync
make -j --output-sync check
#cd tune
#make speed
ls -a
