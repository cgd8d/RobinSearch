echo "JOB_START_TIME=`date +%s`" >> $GITHUB_ENV
sudo apt-get update
sudo apt-get install gnuplot llvm-15 clang-15 autoconf-archive libomp-15-dev libboost-serialization-dev
perf config llvm.clang-path=clang-15
chmod a+x addr2line
chmod a+x objdump
PATH=`pwd`:$PATH which addr2line
PATH=`pwd`:$PATH which objdump

bash scripts/SetupPrimeSieve.sh
bash scripts/SetupGMP.sh
bash scripts/SetupMPFR.sh
bash scripts/CompileRobinSearch.sh
