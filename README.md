# A Partitioning-Based CAD Flow for Interposer-Based Multi-Die FPGAs

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

---
