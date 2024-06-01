# Install GMP from source.

wget https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.xz
tar -xf gmp-6.3.0.tar.xz
cd gmp-6.3.0
./configure --help
CC=clang-$MY_LLVM_VER CFLAGS="-march=native -O2" AR=llvm-ar-$MY_LLVM_VER NM=llvm-nm-$MY_LLVM_VER RANLIB=llvm-ranlib-$MY_LLVM_VER ./configure --disable-shared --disable-fft
make -j --output-sync
make -j --output-sync check
#cd tune
#make speed
ls -a
