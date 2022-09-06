# Install MPFR from source.
# Parallel build doesn't seem to speed it up (6/18/2022),
# but left here anyway.

git clone --depth 1 --branch 4.1 https://gitlab.inria.fr/mpfr/mpfr
cd mpfr
autoreconf -i
./configure --help
CC=clang-14 CFLAGS="-fuse-ld=gold -flto=thin -march=native -O3" AR=llvm-ar-14 NM=llvm-nm-14 RANLIB=llvm-ranlib-14 ./configure --with-gmp-build=${PWD%/*}/gmp-6.2.1 --disable-shared --disable-decimal-float --disable-float128
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
