# Compile libprimesieve from source.
# llvm 15 LTO seems to break something, but I
# doubt there's much to gain from LTO for this
# user case anyway so just turn off LTO.
# gcc LTO does seem to work, if I needed it later.

git clone --depth 1 --branch v11.0 https://github.com/kimwalisch/primesieve
cd primesieve
sed -i 's/private/public/g' include/primesieve/iterator.hpp
CC="clang-15 -march=native" CXX="clang++-15 -march=native" cmake -DBUILD_TESTS=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .
make -j VERBOSE=1
ctest

echo "Run ls"
ls
echo "Run ls CMakeFiles"
ls
