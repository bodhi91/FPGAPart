# A Partitioning-Based CAD Flow for Interposer-Based Multi-Die FPGAs

This repository exists for the sole purpose of a double-blind review for FCCM'25. 

`FPGAPart` adopts the multilevel partitioning paradigm inspired by `TritonPart`. This flow consists of three main phases:

## 1. Multilevel Clustering
A sequence of progressively coarser hypergraphs is constructed. At each level of coarsening, clusters of vertices are identified and merged into a single vertex, representing the cluster in the coarser hypergraph. The input to this flow is a primitive netlist. To incorporate a hint of heterogeneity that can guide the partitioning process, we run VTR 7.0's pre-packing before partitioning. Our experiments demonstrate that incorporating pre-packing information during partitioning improves the quality of results compared to no pre-packing before partitioning.

We extract the top |P| timing-critical paths using VTR’s STA tool. From the P paths, we identify frequently occurring vertex sets (patterns) using pattern-mining techniques. These vertex sets serve as additional clustering guidance. Since these sets include vertices that frequently appear in P, clustering them together can reduce the snaking factor. We use a parallelized FP-Growth algorithm for this. After completion of FP-Growth, the rest of the multilevel hierarchy is constructed by running TritonPart’s multilevel timing-driven clustering.

## 2. Initial Partitioning
After clustering, an initial partitioning solution is computed for the coarsest hypergraph. The reduced size of this hypergraph enables the use of various partitioning methods:
- **Random**
- **VILE**
- **ILP**

We use CPLEX as our ILP solver and leverage the `UserCallback` API to dynamically add cutting planes.

## 3. Multilevel Refinement
Starting with a feasible solution on the coarsest hypergraph, the flow performs uncoarsening and move-based refinement to improve the partitioning solution. This is carried out level by level, gradually refining the solution as the hypergraph is uncoarsened. Two variants of refinement are used:
- **Direct K-way FM**
- **K-way Pairwise FM (K-way P-FM)**

## FPGAPart Overview
The overall FPGAPart framework is illustrated below:
![FPGAPart Framework](https://github.com/bodhi91/FPGAPart/blob/main/fpga_part.png)

## Flow Overview
The overall FPGAPart-based flow is illustrated below:
![Overall Flow](https://github.com/bodhi91/FPGAPart/blob/main/vtr_flow_fccm.png)

## Running Instructions 
Follow these steps to set up and run FPGAPart:

```bash
# Step 1: Download OpenROAD
# Step 2: Clone the OpenROAD repository
$ git clone --recursive https://github.com/The-OpenROAD-Project/OpenROAD.git

# Step 3: Copy contents of FPGAPart to OpenROAD/src/par/
# (Ensure you adjust paths appropriately)

# Step 4: Create a build directory and navigate into it
$ mkdir build && cd build

# Step 5: Run cmake
$ cmake ..

# Step 6: Build the project
$ make -j

# Step 7: Clone VTR 7.0 interposer branch
$ git clone --branch interposer https://github.com/verilog-to-routing/vtr-verilog-to-routing.git

# Step 8: Copy all contents of vtr_7_0 to vtr-verilog-to-routing/vpr/SRC/

# Step 9: Run any of the scripts in the `run_scripts` directory
# (Update path variables in the scripts before running)
```

## Example Run
Below is an example output from running the FPGAPart flow:

```
OpenROAD v2.0-16943-gf2c7a7f 
Features included (+) or not (-): +Charts +GPU +GUI +Python
This program is licensed under the BSD-3 license. See the LICENSE file for details.
Components of this program may be licensed under more restrictive licenses which must be honored.
[INFO] Hypergraph file : logical_hypergraph_tpart.txt
[INFO] Timing file : logical.top_n_paths
[STATUS] Read 1000 timing paths from file
[INFO] Total prepacked clusters : 193
[DEBUG PAR-multilevel_partitioning] Starting multilevel partitioning.
[DEBUG PAR-multilevel_partitioning] no base balance is specified. Use default value.
[DEBUG PAR-multilevel_partitioning] no scale factor is specified. Use default value.
[DEBUG PAR-multilevel_partitioning] No vertex weighting is specified. Use default value of 1.
[DEBUG PAR-multilevel_partitioning] No placement weighting is specified. Use default value of 1.

Multilevel Partitioning Parameters:

hyperedge weight factor : [ 1.000000  ]
vertex weight factor : [ 1.000000 1.000000 1.000000 1.000000 1.000000  ]
placement weight factor : [ 1.000000 1.000000  ]

net_timing_factor : 1.0
path_timing_factor : 1.0
path_snaking_factor : 1.0
timing_exp_factor : 2.0

coarsen order : RANDOM
thr_coarsen_hyperedge_size_skip : 200
thr_coarsen_vertices : 10
thr_coarsen_hyperedges : 50
coarsening_ratio : 1.6
max_coarsen_iters : 30
adj_diff_ratio : 0.0001
min_num_vertices_each_part : 4

num_initial_solutions : 50
num_best_initial_solutions : 10

refine_iters : 10
max_moves (FM or greedy refinement) : 60
early_stop_ratio : 0.5
total_corking_passes : 25
v_cycle_flag : true
max_num_vcycle : 1
num_coarsen_solutions : 3
num_vertices_threshold_ilp : 50

[INFO] Running FP-Growth with min_support = 10
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 472, num_hyperedges = 678, num_timing_paths = 354
[DEBUG PAR-coarsening] Level 1 :: num_vertices = 295, num_hyperedges = 497, num_timing_paths = 290
[DEBUG PAR-coarsening] Level 2 :: num_vertices = 184, num_hyperedges = 314, num_timing_paths = 199
[DEBUG PAR-coarsening] Level 3 :: num_vertices = 115, num_hyperedges = 179, num_timing_paths = 114
[DEBUG PAR-coarsening] Level 4 :: num_vertices = 90, num_hyperedges = 110, num_timing_paths = 47
[DEBUG PAR-coarsening] Level 5 :: num_vertices = 87, num_hyperedges = 103, num_timing_paths = 43
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.014511137 seconds
...
[DEBUG PAR-multilevel_partitioning] The runtime of multilevel partitioner : 9.185599043 seconds
```
---
