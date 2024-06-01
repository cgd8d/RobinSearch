# Compile RobinSearch including profile-guided
# optimization.

clang++-$MY_LLVM_VER --version
/usr/bin/ld -v
lscpu
# Get macros only.
clang++-$MY_LLVM_VER -march=native -E - -dM </dev/null # print macros
# Compile and run instrumented version for PGO.
clang++-$MY_LLVM_VER -fprofile-instr-generate -std=c++20 -fuse-ld=gold -g -Wall -O3 -flto=thin -march=native -Wl,-plugin-opt=whole-program-visibility -fopenmp $MY_DEPENDENCY_FLAGS Search.cpp -lmpfr -lgmp -lboost_serialization -lprimesieve -o Search
echo "Starting instrumented run."
PROF_JOB_STOP_TIME=`date --date="+2 mins" +%s`
time ./Search 60 $PROF_JOB_STOP_TIME checkpoints/Checkpoint_42hrs_20231211.txt
llvm-profdata-$MY_LLVM_VER merge -output=code.profdata default.profraw
# Now compile the real version.
clang++-$MY_LLVM_VER -fprofile-instr-use=code.profdata -std=c++20 -fuse-ld=gold -g -Wall -O3 -flto=thin -march=native -Wl,-plugin-opt=whole-program-visibility -fopenmp $MY_DEPENDENCY_FLAGS Search.cpp -lmpfr -lgmp -lboost_serialization -lprimesieve -o Search
#clang++ -std=c++20 -DNDEBUG -Wall -O3 -flto -march=native -lmpfr -lprimesieve Search.cpp -o Search_ndebug
