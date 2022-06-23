# Compile libprimesieve from source.

git clone --depth 1 https://github.com/kimwalisch/primesieve
cd primesieve
CC="clang-12 -fuse-ld=gold -march=native" CXX="clang++-12 -fuse-ld=gold -march=native" cmake -DBUILD_TESTS=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=TRUE -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .
make -j
ctest
echo "Run ls"
ls
echo "Run ls CMakeFiles"
ls
