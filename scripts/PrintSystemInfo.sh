echo "Script to print system info and exit."

# sudo apt-get install hwloc

echo "lstopo-no-graphics says:"
lstopo-no-graphics -p
echo
echo "lscpu says:"
lscpu --extended --output-all
