# Install MPFR from source.
# Parallel build doesn't seem to speed it up (6/18/2022),
# but left here anyway.

git clone --depth 1 --branch 4.2 https://gitlab.inria.fr/mpfr/mpfr
cd mpfr
autoreconf -i
./configure --help
CC=clang-15 CFLAGS="-fuse-ld=gold -flto=thin -march=native -O3" AR=llvm-ar-15 NM=llvm-nm-15 RANLIB=llvm-ranlib-15 ./configure --with-gmp-build=${PWD%/*}/gmp-6.3.0 --disable-shared --disable-decimal-float --disable-float128
#echo "make tune"
#cd tune
#make tune
#cd ..
echo make
make -j --output-sync
echo make check
make -j --output-sync check
echo ls -a
ls -a
echo ls -a src
cd src
ls -a
cd .libs
echo ls -a src/.libs
ls -a
