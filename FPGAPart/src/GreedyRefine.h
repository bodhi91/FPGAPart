#pragma once

#include <set>

#include "Refiner.h"

namespace par {

class GreedyRefine;
using GreedyRefinerPtr = std::shared_ptr<GreedyRefine>;

// ------------------------------------------------------------------------------
// K-way hyperedge greedy refinement
// Basically during hyperedge greedy refinement, we try to move the straddle
// hyperedge into some block, to minimize the cost.
// Moving the entire hyperedge can help to escape the local minimum caused
// by moving vertex one by one
// ------------------------------------------------------------------------------
class GreedyRefine : public Refiner
{
 public:
  using Refiner::Refiner;

 private:
  // In each pass, we only move the boundary vertices
  // here we pass block_balance and net_degrees as reference
  // because we only move a few vertices during each pass
  // i.e., block_balance and net_degs will not change too much
  // so we precompute the block_balance and net_degs
  // the return value is the gain improvement
  float Pass(const HGraphPtr& hgraph,
             const Matrix<float>& upper_block_balance,
             const Matrix<float>& lower_block_balance,
             Matrix<float>& block_balance,        // the current block balance
             Matrix<int>& net_degs,               // the current net degree
             std::vector<float>& cur_paths_cost,  // the current path cost
             Partitions& solution,
             std::vector<bool>& visited_vertices_flag) override;
};

}  // namespace par