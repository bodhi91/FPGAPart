// High-level description
// This file includes the basic utility functions for operations
///////////////////////////////////////////////////////////////////////////////
#pragma once
#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef LOAD_CPLEX
// for ILP solver in CPLEX
#include "ilcplex/cplex.h"
#include "ilcplex/ilocplex.h"
#endif

namespace par {

// Matrix is a two-dimensional vectors
template <typename T>
using Matrix = std::vector<std::vector<T>>;

struct Rect
{
  // all the values are in db unit
  int lx = 0;
  int ly = 0;
  int ux = 0;
  int uy = 0;

  Rect(int lx_, int ly_, int ux_, int uy_) : lx(lx_), ly(ly_), ux(ux_), uy(uy_)
  {
  }

  // check if the Rect is valid
  bool IsValid() const { return ux > lx && uy > ly; }

  // reset the fence
  void Reset()
  {
    lx = 0;
    ly = 0;
    ux = 0;
    uy = 0;
  }
};

// Define the type for vertices
enum VertexType
{
  COMB_STD_CELL,  // combinational standard cell
  SEQ_STD_CELL,   // sequential standard cell
  MACRO,          // hard macros
  PORT            // IO ports
};

std::string GetVectorString(const std::vector<float>& vec);

// Split a string based on deliminator : empty space and ","
std::vector<std::string> SplitLine(const std::string& line);

// Add right vector to left vector
void Accumulate(std::vector<float>& a, const std::vector<float>& b);

// weighted sum
std::vector<float> WeightedSum(const std::vector<float>& a,
                               float a_factor,
                               const std::vector<float>& b,
                               float b_factor);

// divide the vector
std::vector<float> DivideFactor(const std::vector<float>& a, float factor);

// divide the vectors element by element
std::vector<float> DivideVectorElebyEle(const std::vector<float>& emb,
                                        const std::vector<float>& factor);

// multiplty the vector
std::vector<float> MultiplyFactor(const std::vector<float>& a, float factor);

// operation for two vectors +, -, *,  ==, <
std::vector<float> operator+(const std::vector<float>& a,
                             const std::vector<float>& b);

std::vector<float> operator*(const std::vector<float>& a, float factor);

std::vector<float> operator-(const std::vector<float>& a,
                             const std::vector<float>& b);

std::vector<float> operator*(const std::vector<float>& a,
                             const std::vector<float>& b);

bool operator<(const std::vector<float>& a, const std::vector<float>& b);

bool operator<=(const Matrix<float>& a, const Matrix<float>& b);

bool operator==(const std::vector<float>& a, const std::vector<float>& b);

// Basic functions for a vector
std::vector<float> abs(const std::vector<float>& a);

float norm2(const std::vector<float>& a);

float norm2(const std::vector<float>& a, const std::vector<float>& factor);

// ILP-based Partitioning Instance
// Call ILP Solver to partition the design
bool ILPPartitionInst(
    int num_parts,
    int vertex_weight_dimension,
    std::vector<int>& solution,
    const std::map<int, int>& fixed_vertices,     // vertex_id, block_id
    const Matrix<int>& hyperedges,                // hyperedges
    const std::vector<float>& hyperedge_weights,  // one-dimensional
    const Matrix<float>& vertex_weights,          // two-dimensional
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance);

// Call CPLEX to solve the ILP Based Partitioning
#ifdef LOAD_CPLEX

class HgrLazyConstraintCallback : public IloCplex::LazyConstraintCallbackI
{
 public:
  HgrLazyConstraintCallback(IloEnv env,
                            IloArray<IloNumVarArray> x_vars,
                            int num_parts,
                            int num_vertices,
                            int vertex_dimensions,
                            Matrix<float> vertex_weights,
                            Matrix<float> upper_block_balance,
                            Matrix<float> lower_block_balance)
      : IloCplex::LazyConstraintCallbackI(env),
        x_vars_(x_vars),
        num_parts_(num_parts),
        num_vertices_(num_vertices),
        vertex_dimensions_(vertex_dimensions),
        vertex_weights_(vertex_weights),
        upper_block_balance_(upper_block_balance),
        lower_block_balance_(lower_block_balance)
  {
  }

  IloCplex::CallbackI* duplicateCallback() const override
  {
    return new (getEnv()) HgrLazyConstraintCallback(*this);
  }

 protected:
  void main() override
  {
    for (int i = 0; i < vertex_dimensions_; ++i) {
      for (int block_id = 0; block_id < num_parts_; ++block_id) {
        IloExpr balance_expr(getEnv());
        for (int v = 0; v < num_vertices_; ++v) {
          balance_expr += vertex_weights_[v][i] * x_vars_[block_id][v];
        }
        if (getValue(balance_expr) < lower_block_balance_[block_id][i]
            || getValue(balance_expr) > upper_block_balance_[block_id][i]) {
          IloConstraint lower_bound
              = balance_expr >= lower_block_balance_[block_id][i];
          IloConstraint upper_bound
              = balance_expr <= upper_block_balance_[block_id][i];
          add(lower_bound);
          add(upper_bound);
        }
        balance_expr.end();
      }
    }
  }

 private:
  IloArray<IloNumVarArray> x_vars_;
  int vertex_dimensions_;
  int num_parts_;
  int num_vertices_;
  Matrix<float> vertex_weights_;
  Matrix<float> upper_block_balance_;
  Matrix<float> lower_block_balance_;
};

ILOSTLBEGIN
static int depth_counter_ = 0;

class DepthMIPCallback : public IloCplex::MIPCallbackI
{
 public:
  DepthMIPCallback(IloEnv env) : IloCplex::MIPCallbackI(env) {}

  void main() override
  {
    try {
      int nodeDepth = getNnodes() - getNremainingNodes();
      std::cout << "Approximate node depth: " << nodeDepth << std::endl;
    } catch (IloException& e) {
      std::cerr << "Error: " << e.getMessage() << std::endl;
    }
  }
};

class NeighborhoodInfluenceCuttingPlanes : public IloCplex::UserCutCallbackI
{
 public:
  NeighborhoodInfluenceCuttingPlanes(
      IloEnv env,
      IloArray<IloNumVarArray> x_vars,
      IloArray<IloNumVarArray> y_vars,
      const Matrix<int>& hyperedges,
      const Matrix<float>& vertex_weights,
      const std::vector<float>& hyperedge_weights,
      const Matrix<float>& upper_block_balance,
      const Matrix<float>& lower_block_balance,
      float max_vertex_weight,
      float max_hyperedge_weight,
      int num_parts)
      : IloCplex::UserCutCallbackI(env),
        x_vars_(x_vars),
        y_vars_(y_vars),
        hyperedges_(hyperedges),
        vertex_weights_(vertex_weights),
        hyperedge_weights_(hyperedge_weights),
        upper_block_balance_(upper_block_balance),
        lower_block_balance_(lower_block_balance),
        max_vertex_weight_(max_vertex_weight),
        max_hyperedge_weight_(max_hyperedge_weight),
        num_parts_(num_parts)
  {
  }

  IloCplex::CallbackI* duplicateCallback() const override
  {
    return new (getEnv()) NeighborhoodInfluenceCuttingPlanes(*this);
  }

 protected:
  void main() override
  {
    int num_nodes = getNnodes();
    if (num_nodes < 50) {
      return;
    }

    // Parameters for determining indeterminate assignments
    float tolerance = 0.1;            // Threshold for fractional values
    float ambiguity_threshold = 0.5;  // Minimum gap to avoid ambiguity
    // Iterate over vertices
    for (int vertex = 0; vertex < x_vars_[0].getSize(); ++vertex) {
      // Find the maximum and second-maximum values in x_vars_[vertex]
      float max_val = -1.0;
      float second_max_val = -1.0;
      int max_partition = -1;
      std::vector<float> influence_scores;
      for (int partition = 0; partition < x_vars_.getSize(); ++partition) {
        double value = getValue(x_vars_[partition][vertex]);
        if (value > max_val) {
          second_max_val = max_val;
          max_val = value;
          max_partition = partition;
        } else if (value > second_max_val) {
          second_max_val = value;
        }
      }
      // Calculate gap
      float gap = max_val - second_max_val;
      // Determine if the vertex is indeterminate
      if (max_val < 1.0 - tolerance || gap < ambiguity_threshold) {
        for (int partition = 0; partition < x_vars_.getSize(); ++partition) {
          influence_scores.push_back(calculateInfluence(vertex, partition));
        }
        auto max_influence_it = std::max_element(influence_scores.begin(),
                                                 influence_scores.end());
        int target_partition
            = std::distance(influence_scores.begin(), max_influence_it);
        // make this target partition's probability > 50 %
        // make every other partition's probability < 50 %
        for (int partition = 0; partition < x_vars_.getSize(); ++partition) {
          IloExpr cut_expr(getEnv());
          if (partition == target_partition) {
            cut_expr += x_vars_[partition][vertex] - 0.7;
          } else {
            cut_expr += 0.7 - x_vars_[partition][vertex];
          }
          add(cut_expr >= 0.0);
          cut_expr.end();
        }
        //cut_expr = x_vars_[target_partition][vertex] - 0.5;
        //add(cut_expr >= 0.0);
        //cut_expr.end();
      }
    }
  }

 private:
  int findTargetPartition(const std::vector<double>& influence_scores)
  {
    // Partition influence scores to accumulate total influence per partition
    std::vector<double> partition_influences(num_parts_, 0.0);

    // Loop over each variable and add its influence to the target partition
    for (int var_index = 0; var_index < influence_scores.size(); ++var_index) {
      float vertex_wt = std::accumulate(vertex_weights_[var_index].begin(),
                                        vertex_weights_[var_index].end(),
                                        0.0);
      for (int partition_id = 0; partition_id < num_parts_; ++partition_id) {
        double influence = influence_scores[var_index];
        partition_influences[partition_id] += influence * vertex_wt;
      }
    }

    int target_partition
        = std::distance(partition_influences.begin(),
                        std::max_element(partition_influences.begin(),
                                         partition_influences.end()));

    return target_partition;
  }

  float calculateInfluence(int vertex, int partition)
  {
    float pull = 0.0;
    float max_possible_pull = 0.0;  // Upper bound for normalization

    for (int e : hyperedges_.at(vertex)) {
      int size_he = hyperedges_.at(e).size();
      if (size_he > 50) {
        continue;  // Skip large hyperedges
      }

      // Calculate the cut probability for this hyperedge
      float cut_probability = 0.0;
      for (int p = 0; p < num_parts_; ++p) {
        float partition_sum = 0.0;
        for (int neighbor : hyperedges_.at(e)) {
          partition_sum += getValue(x_vars_[p][neighbor]);
        }
        if (partition_sum > 0 && partition_sum < 1.0) {
          cut_probability += partition_sum * (1.0 - partition_sum);
        }
      }

      // Adjust the hyperedge weight by its cut probability
      float hyperedge_weight = hyperedge_weights_.at(e) / max_hyperedge_weight_;
      hyperedge_weight *= cut_probability;

      // Iterate over neighbors to calculate pull
      for (int neighbor : hyperedges_.at(e)) {
        if (neighbor != vertex) {
          float fractional_assignment = getValue(x_vars_[partition][neighbor]);

          // Normalize the vertex weight of the neighbor
          std::vector<float> vertex_weight_nbr = vertex_weights_.at(neighbor);
          float vertex_weight_nbr_sum = std::accumulate(
              vertex_weight_nbr.begin(), vertex_weight_nbr.end(), 0.0);
          float normalized_vertex_weight
              = vertex_weight_nbr_sum / max_vertex_weight_;

          // Calculate pull as the product of fractional assignment, normalized
          // weight, and adjusted hyperedge weight
          pull += fractional_assignment * normalized_vertex_weight
                  * hyperedge_weight;

          // Update the maximum possible pull for normalization
          max_possible_pull += hyperedge_weight;
        }
      }
    }

    // Ensure pull is between 0 and 1 by dividing by the maximum possible pull
    return (max_possible_pull > 0) ? (pull / max_possible_pull) : 0.0;
  }

  IloArray<IloNumVarArray> x_vars_;
  IloArray<IloNumVarArray> y_vars_;
  const Matrix<int>& hyperedges_;
  const Matrix<float>& vertex_weights_;
  const std::vector<float>& hyperedge_weights_;
  const Matrix<float>& upper_block_balance_;
  const Matrix<float>& lower_block_balance_;
  const float max_vertex_weight_;
  const float max_hyperedge_weight_;
  int num_parts_;
  float threshold_ = 0.7;
};

class NeighborhoodInfluenceBranching : public IloCplex::BranchCallbackI
{
 public:
  NeighborhoodInfluenceBranching(IloEnv env,
                                 IloArray<IloNumVarArray> x_vars,
                                 IloArray<IloNumVarArray> y_vars,
                                 const Matrix<int>& hyperedges,
                                 const Matrix<float>& vertex_weights,
                                 const std::vector<float>& hyperedge_weights,
                                 const Matrix<float>& upper_block_balance,
                                 const Matrix<float>& lower_block_balance,
                                 float max_vertex_weight,
                                 float max_hyperedge_weight)
      : IloCplex::BranchCallbackI(env),
        x_vars_(x_vars),
        y_vars_(y_vars),
        hyperedges_(hyperedges),
        vertex_weights_(vertex_weights),
        hyperedge_weights_(hyperedge_weights),
        upper_block_balance_(upper_block_balance),
        lower_block_balance_(lower_block_balance),
        max_vertex_weight_(max_vertex_weight),
        max_hyperedge_weight_(max_hyperedge_weight)
  {
  }
  IloCplex::CallbackI* duplicateCallback() const override
  {
    return (new (getEnv()) NeighborhoodInfluenceBranching(*this));
  }

 protected:
  void main() override
  {
    // get the depth of current tree

    // NodeId node_id = getNodeId();
    // std::cout << "Get the node ID: " << node_id << std::endl;
    //  get the depth of the node
    // NodeData* node = getNodeData();

    if (depth_counter_++ > 10) {
      return;
    }

    std::cout << "Calling Branching Callback at Depth: " << depth_counter_
              << std::endl;

    // Calculate influence scores on a subset of fractional variables
    std::unordered_map<IloIntVar, double, IloIntVarHash, IloIntVarEqual>
        influence_scores;
    for (int i = 0; i < x_vars_.getSize(); ++i) {
      for (int j = 0; j < x_vars_[i].getSize(); ++j) {
        double val = getValue(x_vars_[i][j]);
        if (std::abs(val - std::round(val)) > 1e-6) {
          // Calculate influence for only a subset to reduce computation
          if (influence_scores.size() < 100) {  // Limiting to top 100
                                                // fractional variables
            double influence = calculateInfluence(j, i);
            influence_scores[x_vars_[i][j]] = influence;
          }
        }
      }
    }

    // Find the variable with the maximum influence score
    if (!influence_scores.empty()) {
      auto max_influence_var
          = std::max_element(influence_scores.begin(),
                             influence_scores.end(),
                             [](const std::pair<IloIntVar, double>& a,
                                const std::pair<IloIntVar, double>& b) {
                               return a.second < b.second;
                             })
                ->first;

      std::cout << "Branching on variable with ID: "
                << max_influence_var.getId()
                << " with influence: " << influence_scores[max_influence_var]
                << std::endl;

      IloNum objestimate
          = getObjValue();  // Use current objective value as base estimate

      std::cout << "Current Objective Value: " << objestimate << std::endl;

      /* do not make a branch
       */

      // Branch on the variable with the maximum influence
      makeBranch(max_influence_var,
                 getValue(max_influence_var),
                 IloCplex::BranchUp,
                 objestimate);
      makeBranch(max_influence_var,
                 getValue(max_influence_var),
                 IloCplex::BranchDown,
                 objestimate);
    }
  }

 private:
  // Define a custom hash and equality comparator for IloIntVar
  struct IloIntVarHash
  {
    std::size_t operator()(const IloIntVar& var) const
    {
      return std::hash<int>()(
          var.getId());  // Assuming getId() uniquely identifies the variable
    }
  };

  struct IloIntVarEqual
  {
    bool operator()(const IloIntVar& var1, const IloIntVar& var2) const
    {
      return var1.getId()
             == var2.getId();  // Check for equality based on unique ID
    }
  };

  float calculateInfluence(int vertex,
                           int partition)  // Calculate influence of a vertex
  {
    float pull = 0.0;
    for (int e : hyperedges_.at(vertex)) {
      int size_he = hyperedges_.at(e).size();
      if (size_he > 50) {
        continue;
      }
      float hyperedge_weight = hyperedge_weights_.at(e);
      for (int neighbor : hyperedges_.at(e)) {
        if (neighbor != vertex) {
          float fractional_assignment = getValue(x_vars_[partition][neighbor]);
          // Normalize the vertex weight of the neighbor
          std::vector<float> vertex_weight_nbr = vertex_weights_.at(neighbor);
          float vertex_weight_nbr_sum = std::accumulate(
              vertex_weight_nbr.begin(), vertex_weight_nbr.end(), 0.0);
          float normalized_vertex_weight
              = vertex_weight_nbr_sum / max_vertex_weight_;
          // Calculate pull as the product of fractional assignment, normalized
          // weight, and hyperedge weight
          pull += fractional_assignment * normalized_vertex_weight
                  * hyperedge_weight;
        }
      }
    }
    return pull;
  }

 private:
  IloArray<IloNumVarArray> x_vars_;
  IloArray<IloNumVarArray> y_vars_;
  const Matrix<int>& hyperedges_;
  const Matrix<float>& vertex_weights_;
  const std::vector<float>& hyperedge_weights_;
  const Matrix<float>& upper_block_balance_;
  const Matrix<float>& lower_block_balance_;
  const float max_vertex_weight_;
  const float max_hyperedge_weight_;
};

inline IloCplex::Callback NbrInfluencecallback(
    IloEnv env,
    IloArray<IloNumVarArray> x_vars,
    IloArray<IloNumVarArray> y_vars,
    const Matrix<int>& hyperedges,
    const Matrix<float>& vertex_weights,
    const std::vector<float>& hyperedge_weights,
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance,
    float max_vertex_weight,
    float max_hyperedge_weight)
{
  return IloCplex::Callback(
      new (env) NeighborhoodInfluenceBranching(env,
                                               x_vars,
                                               y_vars,
                                               hyperedges,
                                               vertex_weights,
                                               hyperedge_weights,
                                               upper_block_balance,
                                               lower_block_balance,
                                               max_vertex_weight,
                                               max_hyperedge_weight));
}

bool OptimalPartCplexv2(
    int num_parts,
    int vertex_weight_dimension,
    std::vector<int>& solution,
    const std::map<int, int>& fixed_vertices,     // vertex_id, block_id
    const Matrix<int>& hyperedges,                // hyperedges
    const std::vector<float>& hyperedge_weights,  // one-dimensional
    const Matrix<float>& vertex_weights,          // two-dimensional
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance);

bool OptimalPartCplex(
    int num_parts,
    int vertex_weight_dimension,
    std::vector<int>& solution,
    const std::map<int, int>& fixed_vertices,     // vertex_id, block_id
    const Matrix<int>& hyperedges,                // hyperedges
    const std::vector<float>& hyperedge_weights,  // one-dimensional
    const Matrix<float>& vertex_weights,          // two-dimensional
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance);
#endif

}  // namespace par