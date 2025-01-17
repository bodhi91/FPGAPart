#pragma once

#include <functional>
#include <map>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// This file contains the classes for Coarsening Phase
// Coarsening Phase is to generate a sequence of coarser hypergraph
// We define Coarsening as an operator class.
// It will accept a HGraphPtr (std::shared_ptr<Hypergraph>) as input
// and return a sequence of coarser hypergraphs

#include "Evaluator.h"
#include "Hypergraph.h"
#include "utl/Logger.h"

namespace par {

// a sequence of coarser hypergraphs
using CoarseGraphPtrs = std::vector<HGraphPtr>;

// nickname for shared pointer of Coarsening
class Coarsener;
using CoarseningPtr = std::shared_ptr<Coarsener>;

// the type for vertex ordering
enum class CoarsenOrder
{
  RANDOM,
  DEGREE,
  SIZE,
  DEFAULT
};

// function : convert CoarsenOrder to string
std::string ToString(CoarsenOrder order);

// coarsening class
// during coarsening, all the coarser hypergraph will not have vertex type
// Because the timing information will become messy in the coarser hypergraph
// And we will never perform slack propagation on a coarser hypergraph

class Coarsener
{
 public:
  Coarsener(
      int num_parts,
      int thr_coarsen_hyperedge_size_skip,  // ignore large hyperedge
      int thr_coarsen_vertices,    // the number of vertices of coarsest
                                   // hypergraph
      int thr_coarsen_hyperedges,  // the number of vertices of coarsest
                                   // hypergraph
      float coarsening_ratio,  // coarsening ratio of two adjacent hypergraphs
      int max_coarsen_iters,   // the number of iterations
      float
          adj_diff_ratio,  // the minimum difference of two adjacent hypergraphs
      const std::vector<float>&
          thr_cluster_weight,  // the weight of largest cluster in a hypergraph
      int random_seed,
      CoarsenOrder vertex_order_choice,  // vertex order
      EvaluatorPtr evaluator,            // evaluator to calculate score
      utl::Logger* logger);

  // the function of coarsen a hypergraph
  // The main function pf Coarsener class
  // The input is a hypergraph
  // The output is a sequence of coarser hypergraphs
  // Notice that the input hypergraph is not const,
  // because the hgraphs returned can be edited
  // The timing cost of hgraph will be initialized if it has been not.
  CoarseGraphPtrs LazyFirstChoice(const HGraphPtr& hgraph) const;

  HGraphPtr ContractHypergraph(
      const HGraphPtr& hgraph,
      const std::vector<int>& vertex_cluster_id_vec) const;

  // create a coarser hypergraph based on specified grouping information
  // for each vertex.
  // each vertex has been map its group
  // (1) remove single-vertex hyperedge
  // (2) remove lager hyperedge
  // (3) detect parallel hyperedges
  // (4) handle group information
  // (5) group fixed vertices based on each block
  // group vertices based on group_attr and hgraph->fixed_attr_
  HGraphPtr GroupVertices(
      const HGraphPtr& hgraph,
      const std::vector<std::vector<int>>& group_attr) const;

  // Set the size of hyperedge to skip
  void SetThrCoarsenHyperedgeSizeSkip(int thr_coarsen_hyperedge_size_skip)
  {
    thr_coarsen_hyperedge_size_skip_ = thr_coarsen_hyperedge_size_skip;
    // if the size of a hyperedge is larger than
    // thr_coarsen_hyperedge_size_skip_, then we ignore this
    // hyperedge during coarsening
  }

  void IncreaseRandomSeed() { random_seed_++; }

 private:
  // private functions (utilities)

  // Single-level Coarsening
  // The input is a hypergraph
  // The output is a coarser hypergraph
  // The input hypergraph will be updated
  HGraphPtr Aggregate(const HGraphPtr& hgraph) const;

  // find the vertex matching scheme
  // the inputs are the hgraph and the attributes of clusters
  // the lazy update means that we do not change the hgraph itself,
  // but during the matching process, we do dynamically update
  // placement_attr_c. vertex_weights_c, fixed_attr_c and community_attr_c
  void VertexMatching(
      const HGraphPtr& hgraph,
      std::vector<int>&
          vertex_cluster_id_vec,  // map current vertex_id to cluster_id
      // the remaining arguments are related to clusters
      Matrix<float>& vertex_weights_c,
      std::vector<int>& community_attr_c,
      std::vector<int>& fixed_attr_c,
      Matrix<float>& placement_attr_c) const;

  // order the vertices based on user-specified parameters
  void OrderVertices(const HGraphPtr& hgraph, std::vector<int>& vertices) const;

  // Similar to the VertexMatching,
  // handle group information
  // group fixed vertices based on each block
  // group vertices based on group_attr and hgraph->fixed_attr_
  void ClusterBasedGroupInfo(
      const HGraphPtr& hgraph,
      const std::vector<std::vector<int>>&
          group_attr,  // Please pass by value here because we need to update
                       // the group_attr
      std::vector<int>&
          vertex_cluster_id_vec,  // map current vertex_id to cluster_id
      // the remaining arguments are related to clusters
      Matrix<float>& vertex_weights_c,
      std::vector<int>& community_attr_c,
      std::vector<int>& fixed_attr_c,
      Matrix<float>& placement_attr_c) const;

  // create the contracted hypergraph based on the vertex matching in
  // vertex_cluster_id_vec
  HGraphPtr Contraction(
      const HGraphPtr& hgraph,
      const std::vector<int>&
          vertex_cluster_id_vec,  // map current vertex_id to cluster_id
      // the remaining arguments are related to clusters
      const Matrix<float>& vertex_weights_c,
      const std::vector<int>& community_attr_c,
      const std::vector<int>& fixed_attr_c,
      const Matrix<float>& placement_attr_c) const;

  const int num_parts_ = 2;
  // coarsening related parameters (stop conditions)

  // If the size of a hyperedge is larger than
  // thr_coarsen_hyperedge_size_skip_, then we ignore this hyperedge
  // during coarsening
  int thr_coarsen_hyperedge_size_skip_ = 50;

  // If the size of a hyperedge is larger than
  // thr_coarsen_hyperedge_size_skip_, then we ignore this hyperedge
  // during coarsening
  const int thr_coarsen_vertices_ = 200;

  // The minimum threshold of number of hyperedges in the coarsest
  // hypergraph
  const int thr_coarsen_hyperedges_ = 50;

  // The ratio of number of vertices of adjacent coarse hypergraphs
  const float coarsening_ratio_ = 1.5;

  // Maxinum number of coarsening iterations
  const int max_coarsen_iters_ = 20;

  // The ratio of number of vertices of adjacent coarse hypergraphs if
  // the ratio is less than adj_diff_ratio_, then stop coarsening
  const float adj_diff_ratio_ = 0.01;

  std::vector<float> thr_cluster_weight_;  // the maximum weight of a cluster
  int random_seed_ = 0;
  CoarsenOrder vertex_order_choice_ = CoarsenOrder::RANDOM;
  EvaluatorPtr evaluator_ = nullptr;
  utl::Logger* logger_ = nullptr;
};

// infrastructure for FPGrowth algorithm

using item_t = int;
using transaction_t = std::vector<item_t>;
using dataset_t = std::vector<transaction_t>;
using item_count_map_t = std::unordered_map<item_t, int>;

class FPGrowth
{
 public:
  FPGrowth(int min_support, int num_threads);
  void fit(const dataset_t& transactions);
  void print_frequent_patterns() const;
  std::vector<std::vector<int>> cluster(std::vector<int> comm_attr = {}) const;

 private:
  struct VectorHash
  {
    template <typename T>
    std::size_t operator()(const std::vector<T>& vec) const
    {
      std::size_t hash = 0;
      std::hash<T> hasher;
      for (const auto& elem : vec) {
        hash ^= hasher(elem) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
      }
      return hash;
    }
  };

  struct Node
  {
    item_t item;
    int count;
    std::map<item_t, std::shared_ptr<Node>> children;

    Node(const item_t& item, int count) : item(item), count(count) {}
  };

  item_count_map_t build_frequent_itemsets(const dataset_t& transactions);
  void construct_fp_tree(const dataset_t& transactions,
                         const item_count_map_t& frequent_items);
  void parallel_mine(const std::shared_ptr<Node>& node,
                     const std::vector<item_t>& prefix,
                     int thread_index);
  void mine_fp_tree(
      const std::shared_ptr<Node>& node,
      const std::vector<item_t>& prefix,
      std::unordered_map<std::vector<item_t>, int, VectorHash>& result);

  std::shared_ptr<Node> root_;
  int min_support_;
  int num_threads_;
  std::mutex tree_mutex_;
  std::vector<std::thread> workers_;
  std::vector<std::unordered_map<std::vector<item_t>, int, VectorHash>>
      patterns_;
  std::unordered_set<std::vector<item_t>, VectorHash> global_pattern_set_;
};

using FPGrower = std::shared_ptr<FPGrowth>;

}  // namespace par
