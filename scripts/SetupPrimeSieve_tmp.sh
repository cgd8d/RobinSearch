# Compile libprimesieve from source.

git clone --depth 1 --branch 1bc574635fc755c8e3f0a0e887de0d4459743eb1 https://github.com/kimwalisch/primesieve
cd primesieve
CC="clang-12 -march=native" CXX="clang++-12 -march=native" cmake -DBUILD_TESTS=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .
make -j
ctest
echo "Run ls"
ls
echo "Run ls CMakeFiles"
ls
