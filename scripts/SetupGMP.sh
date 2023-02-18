# Install GMP from source.

wget https://ftp.gnu.org/gnu/gmp/gmp-6.2.1.tar.xz
tar -xf gmp-6.2.1.tar.xz
cd gmp-6.2.1
./configure --help
CC=clang-15 CFLAGS="-fuse-ld=gold -flto=thin -march=native -O3" AR=llvm-ar-15 NM=llvm-nm-15 RANLIB=llvm-ranlib-15 ./configure --disable-shared --disable-fft
make -j --output-sync
make -j --output-sync check
#cd tune
#make speed
ls -a
