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
min_num_vertcies_each_part : 4

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
[DEBUG PAR-refinement] Set the max_move to 87
[DEBUG PAR-refinement] Reset the max_move to 60
[DEBUG PAR-refinement] Reset the refiner_iters to 10
[DEBUG PAR-refinement] Level 1 :: num_vertices = 90, num_hyperedges = 110, cutcost = 204.89392, best_solution_id = 0
[DEBUG PAR-refinement] Level 2 :: num_vertices = 115, num_hyperedges = 179, cutcost = 238.13416, best_solution_id = 0
[DEBUG PAR-refinement] Level 3 :: num_vertices = 184, num_hyperedges = 314, cutcost = 283.64417, best_solution_id = 0
[DEBUG PAR-refinement] Level 4 :: num_vertices = 295, num_hyperedges = 497, cutcost = 284.44174, best_solution_id = 6
[DEBUG PAR-refinement] Level 5 :: num_vertices = 472, num_hyperedges = 678, cutcost = 282.493, best_solution_id = 8
[DEBUG PAR-refinement] Level 6 :: num_vertices = 748, num_hyperedges = 745, cutcost = 298.20822, best_solution_id = 5
[DEBUG PAR-refinement] Level 7 :: num_vertices = 941, num_hyperedges = 938, cutcost = 313.38925, best_solution_id = 5
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 123, num_hyperedges = 133, num_timing_paths = 39
[DEBUG PAR-coarsening] Level 1 :: num_vertices = 100, num_hyperedges = 107, num_timing_paths = 31
[DEBUG PAR-coarsening] Level 2 :: num_vertices = 98, num_hyperedges = 101, num_timing_paths = 23
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.0007470840000000001 seconds
[DEBUG PAR-refinement] Level 1 :: num_vertices = 100, num_hyperedges = 107, cutcost = 187.23839, best_solution_id = 0
[DEBUG PAR-refinement] Level 2 :: num_vertices = 123, num_hyperedges = 133, cutcost = 192.69112, best_solution_id = 0
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 472, num_hyperedges = 678, num_timing_paths = 354
[DEBUG PAR-coarsening] Level 1 :: num_vertices = 295, num_hyperedges = 497, num_timing_paths = 281
[DEBUG PAR-coarsening] Level 2 :: num_vertices = 184, num_hyperedges = 315, num_timing_paths = 191
[DEBUG PAR-coarsening] Level 3 :: num_vertices = 115, num_hyperedges = 160, num_timing_paths = 80
[DEBUG PAR-coarsening] Level 4 :: num_vertices = 93, num_hyperedges = 115, num_timing_paths = 55
[DEBUG PAR-coarsening] Level 5 :: num_vertices = 92, num_hyperedges = 110, num_timing_paths = 54
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.0046548200000000005 seconds
[DEBUG PAR-refinement] Set the max_move to 92
[DEBUG PAR-refinement] Reset the max_move to 60
[DEBUG PAR-refinement] Reset the refiner_iters to 10
[DEBUG PAR-refinement] Level 1 :: num_vertices = 93, num_hyperedges = 115, cutcost = 217.83842, best_solution_id = 0
[DEBUG PAR-refinement] Level 2 :: num_vertices = 115, num_hyperedges = 160, cutcost = 232.31691, best_solution_id = 7
[DEBUG PAR-refinement] Level 3 :: num_vertices = 184, num_hyperedges = 315, cutcost = 266.8295, best_solution_id = 7
[DEBUG PAR-refinement] Level 4 :: num_vertices = 295, num_hyperedges = 497, cutcost = 271.69446, best_solution_id = 6
[DEBUG PAR-refinement] Level 5 :: num_vertices = 472, num_hyperedges = 678, cutcost = 246.47977, best_solution_id = 1
[DEBUG PAR-refinement] Level 6 :: num_vertices = 748, num_hyperedges = 745, cutcost = 268.61902, best_solution_id = 1
[DEBUG PAR-refinement] Level 7 :: num_vertices = 941, num_hyperedges = 938, cutcost = 277.38443, best_solution_id = 6
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 123, num_hyperedges = 131, num_timing_paths = 34
[DEBUG PAR-coarsening] Level 1 :: num_vertices = 109, num_hyperedges = 114, num_timing_paths = 26
[DEBUG PAR-coarsening] Level 2 :: num_vertices = 106, num_hyperedges = 108, num_timing_paths = 24
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.00090764 seconds
[DEBUG PAR-refinement] Level 1 :: num_vertices = 109, num_hyperedges = 114, cutcost = 177.3, best_solution_id = 0
[DEBUG PAR-refinement] Level 2 :: num_vertices = 123, num_hyperedges = 131, cutcost = 183.57533, best_solution_id = 0
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 472, num_hyperedges = 678, num_timing_paths = 354
[DEBUG PAR-coarsening] Level 1 :: num_vertices = 295, num_hyperedges = 490, num_timing_paths = 265
[DEBUG PAR-coarsening] Level 2 :: num_vertices = 184, num_hyperedges = 310, num_timing_paths = 182
[DEBUG PAR-coarsening] Level 3 :: num_vertices = 115, num_hyperedges = 162, num_timing_paths = 91
[DEBUG PAR-coarsening] Level 4 :: num_vertices = 90, num_hyperedges = 116, num_timing_paths = 54
[DEBUG PAR-coarsening] Level 5 :: num_vertices = 88, num_hyperedges = 111, num_timing_paths = 52
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.004830589000000001 seconds
[DEBUG PAR-refinement] Set the max_move to 88
[DEBUG PAR-refinement] Reset the max_move to 60
[DEBUG PAR-refinement] Reset the refiner_iters to 10
[DEBUG PAR-refinement] Level 1 :: num_vertices = 90, num_hyperedges = 116, cutcost = 222.22754, best_solution_id = 3
[DEBUG PAR-refinement] Level 2 :: num_vertices = 115, num_hyperedges = 162, cutcost = 244.06097, best_solution_id = 3
[DEBUG PAR-refinement] Level 3 :: num_vertices = 184, num_hyperedges = 310, cutcost = 278.04462, best_solution_id = 3
[DEBUG PAR-refinement] Level 4 :: num_vertices = 295, num_hyperedges = 490, cutcost = 275.58133, best_solution_id = 3
[DEBUG PAR-refinement] Level 5 :: num_vertices = 472, num_hyperedges = 678, cutcost = 275.27628, best_solution_id = 0
[DEBUG PAR-refinement] Level 6 :: num_vertices = 748, num_hyperedges = 745, cutcost = 300.84888, best_solution_id = 0
[DEBUG PAR-refinement] Level 7 :: num_vertices = 941, num_hyperedges = 938, cutcost = 306.21698, best_solution_id = 3
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 129, num_hyperedges = 131, num_timing_paths = 19
[DEBUG PAR-coarsening] Level 1 :: num_vertices = 126, num_hyperedges = 128, num_timing_paths = 19
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.000583417 seconds
[DEBUG PAR-refinement] Level 1 :: num_vertices = 129, num_hyperedges = 131, cutcost = 186.6535, best_solution_id = 0
[DEBUG PAR-multilevel_partitioning] Finish Candidate Solutions Generation
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 124, num_hyperedges = 128, num_timing_paths = 25
[DEBUG PAR-coarsening] Level 1 :: num_vertices = 117, num_hyperedges = 118, num_timing_paths = 19
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.00042914300000000005 seconds
[DEBUG PAR-refinement] Level 1 :: num_vertices = 124, num_hyperedges = 128, cutcost = 181.95738, best_solution_id = 0
[DEBUG PAR-multilevel_partitioning] Finish Cut-Overlay Clustering and Optimal Partitioning
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 471, num_hyperedges = 596, num_timing_paths = 323
[DEBUG PAR-coarsening] Level 1 :: num_vertices = 294, num_hyperedges = 444, num_timing_paths = 266
[DEBUG PAR-coarsening] Level 2 :: num_vertices = 183, num_hyperedges = 316, num_timing_paths = 202
[DEBUG PAR-coarsening] Level 3 :: num_vertices = 114, num_hyperedges = 180, num_timing_paths = 121
[DEBUG PAR-coarsening] Level 4 :: num_vertices = 86, num_hyperedges = 107, num_timing_paths = 65
[DEBUG PAR-coarsening] Level 5 :: num_vertices = 84, num_hyperedges = 97, num_timing_paths = 56
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.009134893 seconds
[DEBUG PAR-refinement] Level 1 :: num_vertices = 86, num_hyperedges = 107, cutcost = 173.71864, best_solution_id = 0
[DEBUG PAR-refinement] Level 2 :: num_vertices = 114, num_hyperedges = 180, cutcost = 201.713, best_solution_id = 0
[DEBUG PAR-refinement] Level 3 :: num_vertices = 183, num_hyperedges = 316, cutcost = 225.77992, best_solution_id = 0
[DEBUG PAR-refinement] Level 4 :: num_vertices = 294, num_hyperedges = 444, cutcost = 248.44463, best_solution_id = 0
[DEBUG PAR-refinement] Level 5 :: num_vertices = 471, num_hyperedges = 596, cutcost = 257.15015, best_solution_id = 0
[DEBUG PAR-refinement] Level 6 :: num_vertices = 748, num_hyperedges = 745, cutcost = 268.05017, best_solution_id = 0
[DEBUG PAR-refinement] Level 7 :: num_vertices = 941, num_hyperedges = 938, cutcost = 271.45096, best_solution_id = 0
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 93, num_hyperedges = 92, num_timing_paths = 7
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.000199365 seconds
[DEBUG PAR-multilevel_partitioning] Finish Vcycle Refinement
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 471, num_hyperedges = 596, num_timing_paths = 323
[DEBUG PAR-coarsening] Level 1 :: num_vertices = 294, num_hyperedges = 444, num_timing_paths = 266
[DEBUG PAR-coarsening] Level 2 :: num_vertices = 183, num_hyperedges = 316, num_timing_paths = 202
[DEBUG PAR-coarsening] Level 3 :: num_vertices = 114, num_hyperedges = 180, num_timing_paths = 121
[DEBUG PAR-coarsening] Level 4 :: num_vertices = 86, num_hyperedges = 107, num_timing_paths = 65
[DEBUG PAR-coarsening] Level 5 :: num_vertices = 84, num_hyperedges = 97, num_timing_paths = 56
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.009045712000000001 seconds
[DEBUG PAR-refinement] Level 1 :: num_vertices = 86, num_hyperedges = 107, cutcost = 173.71864, best_solution_id = 0
[DEBUG PAR-refinement] Level 2 :: num_vertices = 114, num_hyperedges = 180, cutcost = 201.713, best_solution_id = 0
[DEBUG PAR-refinement] Level 3 :: num_vertices = 183, num_hyperedges = 316, cutcost = 225.77992, best_solution_id = 0
[DEBUG PAR-refinement] Level 4 :: num_vertices = 294, num_hyperedges = 444, cutcost = 248.44463, best_solution_id = 0
[DEBUG PAR-refinement] Level 5 :: num_vertices = 471, num_hyperedges = 596, cutcost = 257.15015, best_solution_id = 0
[DEBUG PAR-refinement] Level 6 :: num_vertices = 748, num_hyperedges = 745, cutcost = 268.05017, best_solution_id = 0
[DEBUG PAR-refinement] Level 7 :: num_vertices = 941, num_hyperedges = 938, cutcost = 271.45096, best_solution_id = 0
[DEBUG PAR-coarsening] Running FC Multilevel Coarsening...
[DEBUG PAR-coarsening] Level 0 :: num_vertices = 86, num_hyperedges = 85, num_timing_paths = 7
[DEBUG PAR-coarsening] Hierarchical coarsening time 0.000179941 seconds

Print Statistics for Partition

Cutcost : 1140.1669
Vertex balance of block_0 : 0.40123  ( 65.00000 )    1.00000  ( 96.00000 )    1.00000  ( 193.00000 )    1.00000  ( 5.00000 )    0.97320  ( 472.00000 )    
Vertex balance of block_1 : 0.30864  ( 50.00000 )    0.00000  ( 0.00000 )    0.00000  ( 0.00000 )    0.00000  ( 0.00000 )    0.00000  ( 0.00000 )    
Vertex balance of block_2 : 0.29012  ( 47.00000 )    0.00000  ( 0.00000 )    0.00000  ( 0.00000 )    0.00000  ( 0.00000 )    0.02680  ( 13.00000 )    

Constraints and Cut Evaluation

Satisfy the balance constraint : false
Satisfy the group constraint : true
Satisfy the fixed vertices constraint : true
Display Timing Path Cuts Statistics
	Total number of timing paths = 1000
	Total number of timing-critical paths = 0
	Total number of timing-noncritical paths = 1000
	The worst number of cuts on timing-critical paths = 0
	The average number of cuts on timing-critical paths = 0.0
	Total number of timing-noncritical to timing critical paths = 24
	The worst number of cuts on timing-non2critical paths = 2
	The average number of cuts on timing-non2critical paths = 1.5416666
[DEBUG PAR-multilevel_partitioning] The runtime of multilevel partitioner : 9.185599043 seconds
```
---
