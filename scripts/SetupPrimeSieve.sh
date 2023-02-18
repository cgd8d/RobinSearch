# Compile libprimesieve from source.

git clone --depth 1 --branch v11.0 https://github.com/kimwalisch/primesieve
cp -f next_prime2.c primesieve/test/next_prime2.c
cd primesieve
sed -i 's/private/public/g' include/primesieve/iterator.hpp
CC="clang-15 -march=native" CXX="clang++-15 -march=native" cmake -DBUILD_TESTS=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=TRUE -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .
# CC="clang-15 -fuse-ld=gold -march=native" CXX="clang++-15 -fuse-ld=gold -march=native" cmake -DBUILD_TESTS=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=TRUE -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .
make -j VERBOSE=1
#ctest
echo "starting next_prime2"
./test/next_prime2

echo "Run ls"
ls
echo "Run ls CMakeFiles"
ls
