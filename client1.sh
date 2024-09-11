#!/bin/bash

# Compile the C++ program
g++ -o client pollclient.cpp

run() {
    local commands=$1
    ./client localhost << EOL
$commands
EOL
}

# Add the graph
run "ADD_GRAPH 1 5 0 1 2 0 3 6 1 3 8 1 4 5 1 2 3 2 4 7 3 4 9\n"

# Sleep for 1 second to ensure the graph is added before asking for algorithms
sleep 1

# Solve MST using pipeline algorithm
"SOLVE_MST_PIPELINE 1 PRIM\n"

# Sleep for 1 second to ensure the previous command has been processed
sleep 1

# Solve MST using leader-follower pattern with different algorithms
run "SOLVE_MST_LF 1 PRIM"
sleep 1
run "SOLVE_MST_LF 1 KRUSKAL"
sleep 1
run "SOLVE_MST_LF 1 BORUVKA"
sleep 1
run "SOLVE_MST_LF 1 TARJAN"


# Print message to indicate the script has finished
echo "Commands executed."
