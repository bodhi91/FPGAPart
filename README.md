# A Partitioning-Based CAD Flow for Interposer-Based Multi-Die FPGAs

`FPGAPart` adopts the multilevel partitioning paradigm inspired by `TritonPart`. This flow consists of three main phases:

## 1. Multilevel Clustering
A sequence of progressively coarser hypergraphs is constructed. At each level of coarsening, clusters of vertices are identified and merged into a single vertex, representing the cluster in the coarser hypergraph.

## 2. Initial Partitioning
After clustering, an initial partitioning solution is computed for the coarsest hypergraph. The reduced size of this hypergraph enables the use of various partitioning methods:
- **Random**
- **VILE**
- **ILP**

## 3. Multilevel Refinement
Starting with a feasible solution on the coarsest hypergraph, the flow performs uncoarsening and move-based refinement to improve the partitioning solution. This is carried out level by level, gradually refining the solution as the hypergraph is uncoarsened. Two variants of refinement are used:
- **Direct K-way FM**
- **K-way Pairwise FM (K-way P-FM)**

## Flow Overview
The overall flow is illustrated below:
![FPGAPart Framework](https://github.com/bodhi91/FPGAPart/blob/main/fpga_part.png)

---

Feel free to explore the repository for more details on implementation and examples.



