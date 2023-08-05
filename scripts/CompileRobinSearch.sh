# Compile RobinSearch including profile-guided
# optimization.

clang++-15 --version
/usr/bin/ld -v
lscpu
# Get macros only.
clang++-15 -march=native -E - -dM </dev/null # print macros
# Compile and run instrumented version for PGO.
clang++-15 -fprofile-instr-generate -std=c++20 -fuse-ld=gold -g -Wall -O3 -flto=thin -march=native -Wl,-plugin-opt=whole-program-visibility -fopenmp $MY_DEPENDENCY_FLAGS Search.cpp -lmpfr -lgmp -lboost_serialization -lprimesieve -o Search
echo "Starting instrumented run."
time ./Search 38
llvm-profdata-15 merge -output=code.profdata default.profraw
# Now compile the real version.
clang++-15 -fprofile-instr-use=code.profdata -std=c++20 -fuse-ld=gold -g -Wall -O3 -flto=thin -march=native -Wl,-plugin-opt=whole-program-visibility -fopenmp $MY_DEPENDENCY_FLAGS Search.cpp -lmpfr -lgmp -lboost_serialization -lprimesieve -o Search
#clang++ -std=c++20 -DNDEBUG -Wall -O3 -flto -march=native -lmpfr -lprimesieve Search.cpp -o Search_ndebug
