

pid=0

term_handler() {
    if [ $pid -ne 0 ]; then
        kill -SIGTERM "$pid"
        wait "$pid"
    fi
    exit 0;
}
trap 'echo "received SIGINT"; kill ${!}; term_handler' SIGINT
trap 'echo "received SIGTERM"; kill ${!}; term_handler' SIGTERM
trap 'echo "received SIGHUP"; kill ${!}; term_handler' SIGHUP


sudo perf record --call-graph dwarf -F 10 -o perf.data ./Search 58 &
pid="$!"
while true
do
    tail -f /dev/null & wait ${!}
done
