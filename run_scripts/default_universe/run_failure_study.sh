#!/bin/bash

### Set variables ###
num_cuts_values=(3)  # Add your num_cuts values here
percent_wires_cut_values=(80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96)  # Add values for percent wires cut
delay_increase=1
placer_cost_constant=1
allow_fanin_transfer=on
allow_fanout_transfer=on
allow_bidir_interposer_wires=on
allow_chanx_conn=on
allow_additional_fanin=off
allow_additional_fanout=off
pct_interp_to_drive=0
pct_interp_to_be_driven_by=0
timing_driven=on
routing_failure_predictor=aggressive
astar_fac=1.2
constant_type=0
graph_model=star
graph_edge_weight=1/f
design_min_chan_width=0

logical_paths="logical.top_n_paths"
clb_paths="clb.top_n_paths"

ub_factor=5
seed=0

# Configuration
PARENT_DIR=""
RUN_DIR=$(pwd)
VTR_ROOT="${PARENT_DIR}/vtr-verilog-to-routing"
VTR_ARCH="${VTR_ROOT}/vtr_flow/arch/timing/k6_frac_N10_mem32K_40nm.xml"
PART_RUNNER="${VTR_ROOT}/vpr/vpr_jvds.py"
VTR_EXE="${VTR_ROOT}/vpr/vpr"

VTR_BLIF_LOCATION="${PARENT_DIR}/vtr-verilog-to-routing/vtr_flow/benchmarks/vtr_benchmarks_blif"
KOIOS_BLIF_LOCATION=""
touch "${RUN_DIR}/design_vtr_run_cmds"

vtr_designs=("bgm" "mcml" "boundtop" "diffeq1" "diffeq2") 
koios_designs=("eltwise_layer" "conv_layer_hls" "softmax" "reduction_layer" "robot_rl")
vtr_cwidth=(170 198 94 122 114) 
koios_width=(200 94 154 192 172)
total_seeds=20

# For each design, seed, num_cuts, and percent_wires_cut, write out the VTR run command to an output file
for percent_wires_cut in "${percent_wires_cut_values[@]}"; do
    for num_cuts in "${num_cuts_values[@]}"; do
        num_dies=$((num_cuts + 1))  # Calculate the number of dies based on num_cuts

        COMMON_OPTS="--timing_analysis $timing_driven \
                     --timing_driven_clustering $timing_driven \
                     --routing_failure_predictor $routing_failure_predictor \
                     --astar_fac $astar_fac \
                     --num_cuts $num_cuts \
                     --delay_increase $delay_increase \
                     --percent_wires_cut $percent_wires_cut \
                     --placer_cost_constant $placer_cost_constant \
                     --constant_type $constant_type \
                     --allow_bidir_interposer_wires $allow_bidir_interposer_wires \
                     --allow_chanx_conn $allow_chanx_conn \
                     --allow_additional_fanin $allow_additional_fanin \
                     --allow_additional_fanout $allow_additional_fanout \
                     --allow_fanin_transfer $allow_fanin_transfer \
                     --allow_fanout_transfer $allow_fanout_transfer \
                     --pct_interp_to_drive $pct_interp_to_drive"

        for i in "${!vtr_designs[@]}"; do
            design=${vtr_designs[$i]}
            current_channel_width=${vtr_cwidth[$i]}

            for seed in $(seq 0 $((total_seeds - 1))); do
                # Define the output directory structure
                output_dir="${RUN_DIR}/percent_wires_cut_${percent_wires_cut}/num_dies_${num_dies}/${design}/seed_${seed}"

                echo "Writing out VTR script for $design with seed $seed, num_cuts $num_cuts, and percent_wires_cut $percent_wires_cut"

                # Generate command and output it to the script file
                echo "mkdir -p ${output_dir}; cd ${output_dir}; rm -rf *; \
                                ${VTR_EXE} $VTR_ARCH $VTR_BLIF_LOCATION/${design}.blif \
                                --route_chan_width $current_channel_width $COMMON_OPTS --seed $seed --nodisp; \
                                mv vpr_stdout.log ${design}_cuts${num_cuts}_pwcut${percent_wires_cut}_seed${seed}.vpr_default_stdout.log;" >> "${RUN_DIR}/design_vtr_run_cmds"

                # Command to remove all files except the log file
                echo "find ${output_dir}/ -type f ! -name '*.log' -delete; cd ..;" >> "${RUN_DIR}/design_vtr_run_cmds"
            done
        done

        for i in "${!koios_designs[@]}"; do
            design=${koios_designs[$i]}
            current_channel_width=${koios_width[$i]}

            for seed in $(seq 0 $((total_seeds - 1))); do
                # Define the output directory structure
                output_dir="${RUN_DIR}/percent_wires_cut_${percent_wires_cut}/num_dies_${num_dies}/${design}/seed_${seed}"

                echo "Writing out VTR script for $design with seed $seed, num_cuts $num_cuts, and percent_wires_cut $percent_wires_cut"

                # Generate command and output it to the script file
                echo "mkdir -p ${output_dir}; cd ${output_dir}; rm -rf *; \
                                ${VTR_EXE} $VTR_ARCH $KOIOS_BLIF_LOCATION/${design}.abc.blif \
                                --route_chan_width $current_channel_width $COMMON_OPTS --seed $seed --nodisp; \
                                mv vpr_stdout.log ${design}_cuts${num_cuts}_pwcut${percent_wires_cut}_seed${seed}.vpr_default_stdout.log;" >> "${RUN_DIR}/design_vtr_run_cmds"

                # Command to remove all files except the log file
                echo "find ${output_dir}/ -type f ! -name '*.log' -delete; cd ..;" >> "${RUN_DIR}/design_vtr_run_cmds"
            done
        done
    done
done
