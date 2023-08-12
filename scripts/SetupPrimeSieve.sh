# Compile libprimesieve from source.
# llvm 15 LTO seems to break something, but I
# doubt there's much to gain from LTO for this
# user case anyway so just turn off LTO.
# gcc LTO does seem to work, if I needed it later.

git clone --depth 1 --branch v11.1 https://github.com/kimwalisch/primesieve
cd primesieve
sed -i 's/private/public/g' include/primesieve/iterator.hpp
#sed -i 's/BUCKET_BYTES = 8 << 10/BUCKET_BYTES = 8 << 11/' include/primesieve/config.hpp
sed -i 's/FACTOR_ERATSMALL = 0.2/FACTOR_ERATSMALL = 0.15/' include/primesieve/config.hpp
cat include/primesieve/config.hpp
CC="clang-15 -march=native" CXX="clang++-15 -march=native" cmake -DBUILD_TESTS=ON -DBUILD_SHARED_LIBS=OFF -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON .
make -j VERBOSE=1
ctest

echo "Run ls"
ls
echo "Run ls CMakeFiles"
ls
