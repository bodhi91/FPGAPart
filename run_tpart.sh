 #!/bin/bash

# accept the following arguments from command line
OR_EXEC=""
hypergraph_file=$1
num_parts=$2
balance_constraint=$3
seed=$4
timing_file=$5
pre_packed_file=$6

# Create a run_triton_part.tcl file containing the following
cat > run_triton_part.tcl <<EOF
#!/usr/bin/env tclsh

# Execute the triton_part_hypergraph with passed parameters
triton_part_fpga -hypergraph_file $hypergraph_file -num_parts $num_parts -balance_constraint $balance_constraint -timing_file $timing_file -pre_packed_file $pre_packed_file -seed $seed
exit
EOF

$OR_EXEC run_triton_part.tcl | tee run_triton_part.log
