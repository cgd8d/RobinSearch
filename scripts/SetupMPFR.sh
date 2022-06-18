# Install MPFR from source.

git clone --depth 1 --branch 4.1 https://gitlab.inria.fr/mpfr/mpfr
cd mpfr
autoreconf -i
./configure --help
CC=clang-12 CFLAGS="-fuse-ld=gold -flto=thin -march=native -O3" AR=llvm-ar-12 NM=llvm-nm-12 RANLIB=llvm-ranlib-12 ./configure --with-gmp-build=${PWD%/*}/gmp-6.2.1 --disable-shared --disable-decimal-float --disable-float128
#echo "make tune"
#cd tune
#make tune
#cd ..
echo make
make
echo make check
make check
echo ls -a
ls -a
echo ls -a src
cd src
ls -a
cd .libs
echo ls -a src/.libs
ls -a
