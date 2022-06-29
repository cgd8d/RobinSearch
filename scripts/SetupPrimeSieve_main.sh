# Compile libprimesieve from source.

git clone --depth 1 https://github.com/kimwalisch/primesieve
cd primesieve
CC="clang-12 -march=native -Wall -Wextra" CXX="clang++-12 -march=native -Wall -Wextra" cmake -DBUILD_TESTS=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .
make -j
ctest
echo "Run ls"
ls
echo "Run ls CMakeFiles"
ls
