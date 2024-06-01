sudo apt-get update
sudo apt-get install gnuplot llvm-$MY_LLVM_VER clang-$MY_LLVM_VER autoconf-archive libomp-$MY_LLVM_VER-dev libboost-serialization-dev hwloc
perf config llvm.clang-path=clang-$MY_LLVM_VER
chmod a+x addr2line
chmod a+x objdump
PATH=`pwd`:$PATH which addr2line
PATH=`pwd`:$PATH which objdump

echo "lstopo-no-graphics says:"
lstopo-no-graphics -p

echo "Setting up libprimesieve..."
time bash scripts/SetupPrimeSieve.sh

echo "Setting up GMP..."
time bash scripts/SetupGMP.sh

echo "Setting up MPFR..."
time bash scripts/SetupMPFR.sh

echo "Compiling RobinSearch..."
time bash scripts/CompileRobinSearch.sh
