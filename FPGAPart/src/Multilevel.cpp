
#include "Multilevel.h"

#include <functional>
#include <queue>
#include <random>
#include <thread>

#include "Evaluator.h"
#include "Hypergraph.h"
#include "Partitioner.h"
#include "utl/Logger.h"

namespace par {

using utl::PAR;

MultilevelPartitioner::MultilevelPartitioner(
    const int num_parts,
    const bool v_cycle_flag,
    const int num_initial_solutions,
    const int num_best_initial_solutions,
    const int num_vertices_threshold_ilp,
    const int max_num_vcycle,
    const int num_coarsen_solutions,
    const int seed,
    CoarseningPtr coarsener,
    FPGrower fpgrower,
    PartitioningPtr partitioner,
    KWayFMRefinerPtr k_way_fm_refiner,
    KWayPMRefinerPtr k_way_pm_refiner,
    GreedyRefinerPtr greedy_refiner,
    IlpRefinerPtr ilp_refiner,
    EvaluatorPtr evaluator,
    utl::Logger* logger)
    : num_parts_(num_parts),
      num_vertices_threshold_ilp_(num_vertices_threshold_ilp),
      num_initial_random_solutions_(num_initial_solutions),
      num_best_initial_solutions_(num_best_initial_solutions),
      max_num_vcycle_(max_num_vcycle),
      num_coarsen_solutions_(num_coarsen_solutions),
      seed_(seed),
      v_cycle_flag_(v_cycle_flag)
{
  coarsener_ = std::move(coarsener);
  fpgrower_ = std::move(fpgrower);
  partitioner_ = std::move(partitioner);
  k_way_fm_refiner_ = std::move(k_way_fm_refiner);
  k_way_pm_refiner_ = std::move(k_way_pm_refiner);
  greedy_refiner_ = std::move(greedy_refiner);
  ilp_refiner_ = std::move(ilp_refiner);
  evaluator_ = std::move(evaluator);
  logger_ = logger;
}

// Main function
// here the hgraph should not be const
// Because our slack-rebudgeting algorithm will change hgraph
std::vector<int> MultilevelPartitioner::Partition(
    const HGraphPtr& hgraph,
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance)
{
  // Main implementation
  // Step 1: run initial partitioning with different random seed
  // Step 2: run cut-overlay clustering to enhance the solution
  // Step 3: Guided v-cycle

  // Step 1: run initial partitioning with different random seed
  // In experiments, we observe that the benefits of increasing number of
  // vcycles is very limited. However, the quality of solutions will change a
  // lot with different random seed
  Matrix<int> top_solutions;
  float best_cost = std::numeric_limits<float>::max();
  int best_solution_id = -1;
  for (int id = 0; id < num_coarsen_solutions_; id++) {
    coarsener_->IncreaseRandomSeed();
    SetFPGrowth(true);
    //SetPrePack(true);
    top_solutions.push_back(
        SingleLevelPartition(hgraph, upper_block_balance, lower_block_balance));
    const float cost
        = evaluator_->CutEvaluator(hgraph, top_solutions.back(), false).cost;
    if (cost <= best_cost) {
      best_cost = cost;
      best_solution_id = id;
    }
  }

  debugPrint(logger_,
             PAR,
             "multilevel_partitioning",
             1,
             "Finish Candidate Solutions Generation");

  // Step 2: run cut-overlay clustering to enhance the solution
  SetFPGrowth(false);
  //SetPrePack(false);
  std::vector<int> best_solution = CutOverlayILPPart(hgraph,
                                                     upper_block_balance,
                                                     lower_block_balance,
                                                     top_solutions,
                                                     best_solution_id);

  debugPrint(logger_,
             PAR,
             "multilevel_partitioning",
             1,
             "Finish Cut-Overlay Clustering and Optimal Partitioning");

  // Step 3: Guided v-cycle. Note that hgraph has been updated.
  // The best_solution will be refined.
  // The initial value of best solution will be used to guide the coarsening
  // process and use as the initial solution
  if (v_cycle_flag_ == true) {
    VcycleRefinement(
        hgraph, upper_block_balance, lower_block_balance, best_solution);
  }

  debugPrint(
      logger_, PAR, "multilevel_partitioning", 1, "Finish Vcycle Refinement");

  return best_solution;
}

// Private functions (Utilities)

// Run single-level partitioning
std::vector<int> MultilevelPartitioner::SingleLevelPartition(
    const HGraphPtr& hgraph,
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance)
{
  // Step 1: run coarsening
  // Step 2: run initial partitioning
  // Step 3: run refinement
  // Step 4: cut-overlay clustering and ILP-based partitioning

  CoarseGraphPtrs hierarchy;
  HGraphPtr clustered_hgraph = hgraph;

  // Step 1: Run initial coarsening if pre_pack_ is enabled
  if (pre_pack_) {
    hierarchy.push_back(hgraph);
    clustered_hgraph
        = coarsener_->ContractHypergraph(hgraph, pre_packed_solution_);
    //hierarchy.push_back(clustered_hgraph);
  }

  // Step 2: Run FP-Growth-based clustering if enabled
  if (fp_growth_) {
    std::vector<std::vector<int>> transactions;
    if (pre_pack_) {
      // if pre_packing happened then we need to re-cluster the transactions
      for (auto& transaction : transactions_) {
        std::set<int> new_transaction;
        for (auto& vertex : transaction) {
          new_transaction.insert(pre_packed_solution_[vertex]);
        }
        if (new_transaction.size() > 1) {
          transactions.push_back(std::vector<int>(new_transaction.begin(),
                                                new_transaction.end()));
        }
      }
      if (clustered_hgraph->HasTiming() && clustered_hgraph->GetTimingPathCostSize() == 0
        && !clustered_hgraph->HasHyperedgeTimingCost()) {
        evaluator_->InitializeTiming(clustered_hgraph);
      } 
    } else {
      transactions = transactions_;
    }

    fpgrower_->fit(transactions);

    // Step 2.1: Perform FP-Growth clustering
    std::vector<std::vector<int>> clusters = fpgrower_->cluster();
    std::vector<int> vertex_cluster_id_vec(clustered_hgraph->GetNumVertices(),
                                           -1);
    int cluster_id = 0;  // the id of cluster

    for (const auto& group : clusters) {
      for (const auto& v : group) {
        vertex_cluster_id_vec[v] = cluster_id;
      }
      cluster_id++;
    }

    // assign each other vertex to a single-vertex cluster
    for (int v = 0; v < clustered_hgraph->GetNumVertices(); v++) {
      if (vertex_cluster_id_vec[v] == -1) {
        vertex_cluster_id_vec[v] = cluster_id++;
      }
    }

    // Step 2.2: Contract hypergraph based on FP-Growth clustering
    auto clustered_hgraph_fp = coarsener_->ContractHypergraph(
        clustered_hgraph, vertex_cluster_id_vec);

    // Step 2.3: Perform additional lazy-first-choice coarsening
    hierarchy.push_back(clustered_hgraph);
    CoarseGraphPtrs fc_hierarchy
        = coarsener_->LazyFirstChoice(clustered_hgraph_fp);

    // Step 2.4: Add all levels of hierarchy to the final hierarchy
    for (auto& hgraph : fc_hierarchy) {
      hierarchy.push_back(hgraph);
    }
  } else {
    // Step 3: Default lazy-first-choice coarsening
    hierarchy = coarsener_->LazyFirstChoice(hgraph);
  }

  // CoarseGraphPtrs hierarchy = coarsener_->LazyFirstChoice(hgraph);
  // Step 2: run initial partitioning
  HGraphPtr coarsest_hgraph = hierarchy.back();

  // pick top num_best_initial_solutions_ solutions from
  // num_initial_random_solutions_ solutions
  // here we reserve GetNumVertices() for each top_solution
  Matrix<int> top_solutions;
  for (int i = 0; i < num_best_initial_solutions_; i++) {
    std::vector<int> solution{};
    top_solutions.push_back(solution);
    top_solutions.back().reserve(hgraph->GetNumVertices());  // reserve the size
  }
  int best_solution_id = -1;
  InitialPartition(coarsest_hgraph,
                   upper_block_balance,
                   lower_block_balance,
                   top_solutions,
                   best_solution_id);

  // Step 3: run refinement
  // Here we need to do rebugetting on the best solution
  RefinePartition(hierarchy,
                  upper_block_balance,
                  lower_block_balance,
                  top_solutions,
                  best_solution_id);

  // Step 4: cut-overlay clustering and ILP-based partitioning
  // Perform cut-overlay clustering and ILP-based partitioning
  // The ILP-based partitioning uses top_solutions[best_solution_id] as a hint,
  // such that the runtime can be signficantly reduced
  SetFPGrowth(
      false);  // we need to set it to false so that the ILP does not FP-growth
  //SetPrePack(false);
  return CutOverlayILPPart(hgraph,
                           upper_block_balance,
                           lower_block_balance,
                           top_solutions,
                           best_solution_id);
}

// Use the initial solution as the community feature
// Call Vcycle refinement
void MultilevelPartitioner::VcycleRefinement(
    const HGraphPtr& hgraph,
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance,
    std::vector<int>& best_solution)
{
  Matrix<int> candidate_solutions;
  candidate_solutions.push_back(best_solution);
  for (int num_cycles = 0; num_cycles < max_num_vcycle_; num_cycles++) {
    // use the initial solution as the community feature
    hgraph->SetCommunity(best_solution);
    SetFPGrowth(true);
    //SetPrePack(true);
    best_solution = SingleCycleRefinement(
        hgraph, upper_block_balance, lower_block_balance);
    candidate_solutions.push_back(best_solution);
    const float cost
        = evaluator_->CutEvaluator(hgraph, best_solution, false).cost;
    debugPrint(logger_,
               PAR,
               "v_cycle_refinement",
               1,
               "num_cycles = {}, cutcost = {}",
               num_cycles,
               cost);
  }

  SetFPGrowth(false);
  //SetPrePack(false);
  // Perform Cut-overlay clustering and ILP-based partitioning
  best_solution
      = CutOverlayILPPart(hgraph,
                          upper_block_balance,
                          lower_block_balance,
                          candidate_solutions,
                          static_cast<int>(candidate_solutions.size()) - 1);
}

// Single Vcycle Refinement
std::vector<int> MultilevelPartitioner::SingleCycleRefinement(
    const HGraphPtr& hgraph,
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance) const
{
  // Step 1: run coarsening
  // Step 2: run refinement

  // Step 1: run coarsening
  CoarseGraphPtrs hierarchy;
  HGraphPtr clustered_hgraph = hgraph;

  // Step 1: Run initial coarsening if pre_pack_ is enabled
  std::vector<int> community;
  if (hgraph->HasCommunity()) {
    for (int i = 0; i < hgraph->GetNumVertices(); i++) {
      community.push_back(hgraph->GetCommunity(i));
    }
  }

  // use pre-packced solution filtered through community info.
  if (pre_pack_) {
    // re-do clustering based on the community
    std::vector<int> comm_cluster_id(hgraph->GetNumVertices(), -1);
    std::unordered_map<int, std::vector<int>> clusters;

    for (int i = 0; i < hgraph->GetNumVertices(); i++) {
      int cluster_id = pre_packed_solution_[i];
      clusters[cluster_id].push_back(i);
    }

    // now redo: clustering
    // 1. for each cluster see if we have vertices belonging to different
    // communities
    // 2. if yes, split the cluster into separate clusters each containing
    // vertices from the same community

    int cluster_id = 0;
    for (const auto& cluster : clusters) {
      std::unordered_map<int, std::vector<int>> comm_vertices;
      for (const auto& v : cluster.second) {
        int comm = community[v];
        comm_vertices[comm].push_back(v);
      }

      for (const auto& comm_vertex : comm_vertices) {
        for (const auto& v : comm_vertex.second) {
          comm_cluster_id[v] = cluster_id;
        }
        cluster_id++;
      }
    }

    clustered_hgraph = coarsener_->ContractHypergraph(hgraph, comm_cluster_id);
    hierarchy.push_back(hgraph);
  }

  if (fp_growth_ == true) {
    // cluster transactions_ based on pre_packed clusters
    std::vector<std::vector<int>> transactions;
    if (pre_pack_) {
      // if pre_packing happened then we need to re-cluster the transactions
      for (auto& transaction : transactions_) {
        std::set<int> new_transaction;
        for (auto& vertex : transaction) {
          new_transaction.insert(pre_packed_solution_[vertex]);
        }
        transactions.push_back(std::vector<int>(new_transaction.begin(),
                                                new_transaction.end()));
      }

      if (clustered_hgraph->HasTiming() && clustered_hgraph->GetTimingPathCostSize() == 0
        && !clustered_hgraph->HasHyperedgeTimingCost()) {
        evaluator_->InitializeTiming(clustered_hgraph);
      } 
    } else {
      transactions = transactions_;
    }

    fpgrower_->fit(transactions);

    std::vector<std::vector<int>> clusters = fpgrower_->cluster();

    std::vector<int> vertex_cluster_id_vec(clustered_hgraph->GetNumVertices(),
                                           -1);
    int cluster_id = 0;  // the id of cluster

    for (const auto& group : clusters) {
      for (const auto& v : group) {
        vertex_cluster_id_vec[v] = cluster_id;
      }
      cluster_id++;
    }

    // assign each other vertex to a single-vertex cluster
    for (int v = 0; v < clustered_hgraph->GetNumVertices(); v++) {
      if (vertex_cluster_id_vec[v] == -1) {
        vertex_cluster_id_vec[v] = cluster_id++;
      }
    }

    auto clustered_hgraph_fp = coarsener_->ContractHypergraph(
        clustered_hgraph, vertex_cluster_id_vec);
    hierarchy.push_back(clustered_hgraph);
    CoarseGraphPtrs fc_hierarchy
        = coarsener_->LazyFirstChoice(clustered_hgraph_fp);
    hierarchy.insert(hierarchy.end(), fc_hierarchy.begin(), fc_hierarchy.end());
  } else {
    hierarchy = coarsener_->LazyFirstChoice(hgraph);
  }
  // CoarseGraphPtrs hierarchy = coarsener_->LazyFirstChoice(hgraph);

  // Step 2: run initial refinement
  HGraphPtr coarsest_hgraph = hierarchy.back();
  Matrix<int> top_solutions(1);
  coarsest_hgraph->CopyCommunity(top_solutions[0]);
  int best_solution_id = 0;  // only one solution
  if (coarsest_hgraph->GetNumVertices() <= num_vertices_threshold_ilp_) {
    partitioner_->Partition(coarsest_hgraph,
                            upper_block_balance,
                            lower_block_balance,
                            top_solutions[best_solution_id],
                            PartitionType::INIT_DIRECT_ILP);
  }
  // Here we need to do rebugetting on the best solution
  RefinePartition(hierarchy,
                  upper_block_balance,
                  lower_block_balance,
                  top_solutions,
                  best_solution_id);

  return top_solutions[best_solution_id];
}

// Generate initial partitioning
// Include random partitioning, Vile partitioning and ILP partitioning
void MultilevelPartitioner::InitialPartition(
    const HGraphPtr& hgraph,
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance,
    Matrix<int>& top_initial_solutions,
    int& best_solution_id) const
{
  debugPrint(logger_,
             PAR,
             "initial_partitioning",
             1,
             "Running Initial Partitioning...");
  std::mt19937 gen;
  gen.seed(seed_);
  std::uniform_real_distribution<> dist(0.0, 1.0);
  std::vector<float> initial_solutions_cost;
  std::vector<bool>
      initial_solutions_flag;  // if the solutions statisfy balance constraint
  Matrix<int> initial_solutions;
  if (hgraph->GetNumVertices() <= num_vertices_threshold_ilp_) {
    // random partitioning + Vile + ILP
    initial_solutions.resize(num_initial_random_solutions_ * 2 + 2);
  } else {
    // random partitioning + Vile
    initial_solutions.resize(num_initial_random_solutions_ * 2 + 1);
  }
  // We need k_way_fm_refiner to generate a balanced partitioning
  k_way_fm_refiner_->SetMaxMove(hgraph->GetNumVertices());
  // generate random seed
  for (int i = 0; i < num_initial_random_solutions_; ++i) {
    const int seed = std::numeric_limits<int>::max() * dist(gen);
    auto& solution = initial_solutions[i];
    // call random partitioning
    partitioner_->SetRandomSeed(seed);
    partitioner_->Partition(hgraph,
                            upper_block_balance,
                            lower_block_balance,
                            solution,
                            PartitionType::INIT_RANDOM);
    // call FM refiner to improve the solution
    k_way_fm_refiner_->Refine(
        hgraph, upper_block_balance, lower_block_balance, solution);
    const auto token = evaluator_->CutEvaluator(hgraph, solution, false);
    initial_solutions_cost.push_back(token.cost);
    // Here we only check the upper bound to make sure more possible solutions
    initial_solutions_flag.push_back(token.block_balance
                                     <= upper_block_balance);
    debugPrint(logger_,
               PAR,
               "initial_partitioning",
               1,
               "{} :: Random part cutcost = {}, balance_flag = {}",
               i,
               initial_solutions_cost.back(),
               (bool) initial_solutions_flag.back());
  }
  // generate random vile solution
  for (int i = 0; i < num_initial_random_solutions_; ++i) {
    const int seed = std::numeric_limits<int>::max() * dist(gen);
    auto& solution = initial_solutions[i + num_initial_random_solutions_];
    // call random partitioning
    partitioner_->SetRandomSeed(seed);
    partitioner_->Partition(hgraph,
                            upper_block_balance,
                            lower_block_balance,
                            solution,
                            PartitionType::INIT_RANDOM_VILE);
    // call FM refiner to improve the solution
    k_way_fm_refiner_->Refine(
        hgraph, upper_block_balance, lower_block_balance, solution);
    const auto token = evaluator_->CutEvaluator(hgraph, solution, false);
    initial_solutions_cost.push_back(token.cost);
    // Here we only check the upper bound to make sure more possible solutions
    initial_solutions_flag.push_back(token.block_balance
                                     <= upper_block_balance);
    debugPrint(logger_,
               PAR,
               "initial_partitioning",
               1,
               "{} :: Random VILE part cutcost = {}, balance_flag = {}",
               i,
               initial_solutions_cost.back(),
               (bool) initial_solutions_flag.back());
  }

  // Vile partitioning. Vile partitioning needs refiner to generated a balanced
  // partitioning
  auto& vile_solution = initial_solutions[num_initial_random_solutions_ * 2L];
  partitioner_->Partition(hgraph,
                          upper_block_balance,
                          lower_block_balance,
                          vile_solution,
                          PartitionType::INIT_VILE);
  // We need k_way_fm_refiner to generate a balanced partitioning
  k_way_fm_refiner_->Refine(
      hgraph, upper_block_balance, lower_block_balance, vile_solution);
  k_way_fm_refiner_->RestoreDefaultParameters();
  const auto vile_token
      = evaluator_->CutEvaluator(hgraph, vile_solution, false);
  initial_solutions_cost.push_back(vile_token.cost);
  initial_solutions_flag.push_back(vile_token.block_balance
                                   <= upper_block_balance);
  debugPrint(logger_,
             PAR,
             "initial_partitioning",
             1,
             "VILE part cutcost = {}, balance_flag = {}",
             initial_solutions_cost.back(),
             (bool) initial_solutions_flag.back());
  // ILP partitioning
  if (hgraph->GetNumVertices() <= num_vertices_threshold_ilp_) {
    auto& ilp_solution = initial_solutions.back();
    // Use previous best solution as a starting point
    int ilp_solution_id = 0;
    float ilp_solution_cost = std::numeric_limits<float>::max();
    for (auto id = 0; id < initial_solutions_cost.size(); id++) {
      if (initial_solutions_flag[id] == true
          && initial_solutions_cost[id] < ilp_solution_cost) {
        ilp_solution_id = id;
        ilp_solution_cost = initial_solutions_cost[id];
      }
    }
    ilp_solution = initial_solutions[ilp_solution_id];
    partitioner_->Partition(hgraph,
                            upper_block_balance,
                            lower_block_balance,
                            ilp_solution,
                            PartitionType::INIT_DIRECT_ILP);
    const auto ilp_token
        = evaluator_->CutEvaluator(hgraph, ilp_solution, false);
    initial_solutions_cost.push_back(ilp_token.cost);
    initial_solutions_flag.push_back(ilp_token.block_balance
                                     <= upper_block_balance);
    debugPrint(logger_,
               PAR,
               "initial_partitioning",
               1,
               "ILP part cutcost = {}, balance_flag = {}",
               initial_solutions_cost.back(),
               (bool) initial_solutions_flag.back());
  }

  // sort the solutions based on cost
  std::vector<int> solution_ids(initial_solutions_cost.size(), 0);
  std::iota(solution_ids.begin(), solution_ids.end(), 0);
  // define compare function
  auto lambda_sort_criteria = [&](int& x, int& y) -> bool {
    return initial_solutions_cost[x] < initial_solutions_cost[y];
  };
  std::sort(solution_ids.begin(), solution_ids.end(), lambda_sort_criteria);
  // pick the top num_best_initial_solutions_ solutions
  // while satisfying the balance constraint
  int num_chosen_best_init_solution = 0;
  float best_initial_cost = 0.0;
  std::vector<bool> visited_solution_flag(solution_ids.size(), false);
  for (auto id : solution_ids) {
    if (initial_solutions_flag[id] == true) {
      top_initial_solutions[num_chosen_best_init_solution]
          = initial_solutions[id];
      visited_solution_flag[id] = true;
      num_chosen_best_init_solution++;
      if (num_chosen_best_init_solution == 1) {
        best_initial_cost = initial_solutions_cost[id];
      }
      if (num_chosen_best_init_solution >= num_best_initial_solutions_) {
        break;
      }
    }
  }

  if (num_chosen_best_init_solution < num_best_initial_solutions_) {
    for (auto id : solution_ids) {
      if (visited_solution_flag[id] == false) {
        top_initial_solutions[num_chosen_best_init_solution]
            = initial_solutions[id];
        visited_solution_flag[id] = true;
        num_chosen_best_init_solution++;
        if (num_chosen_best_init_solution >= num_best_initial_solutions_) {
          break;
        }
      }
    }
  }

  // remove invalid solution
  while (top_initial_solutions.size() > num_chosen_best_init_solution) {
    top_initial_solutions.pop_back();
  }

  debugPrint(logger_,
             PAR,
             "initial_partitioning",
             1,
             "Number of chosen best initial solutions = {}",
             num_chosen_best_init_solution);
  best_solution_id = 0;  // the first one is the best one
  debugPrint(logger_,
             PAR,
             "initial_partitioning",
             1,
             "Best initial cutcost {}",
             best_initial_cost);
}

// Refine the solutions in top_solutions in parallel with multi-threading
// the top_solutions and best_solution_id will be updated during this process
void MultilevelPartitioner::RefinePartition(
    CoarseGraphPtrs hierarchy,
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance,
    Matrix<int>& top_solutions,
    int& best_solution_id) const
{
  if (hierarchy.size() <= 1) {
    return;  // no need to refine.
  }

  int num_level = 0;
  // rebudget based on the best solution
  auto hgraph_iter = hierarchy.rbegin();
  while (hgraph_iter != hierarchy.rend()) {
    HGraphPtr coarse_hgraph = *hgraph_iter;
    hgraph_iter++;
    if (hgraph_iter == hierarchy.rend()) {
      return;
    }
    HGraphPtr hgraph = *hgraph_iter;

    // convert the solution in coarse_hgraph to the solution of hgraph
    for (auto& top_solution : top_solutions) {
      std::vector<int> refined_solution;
      refined_solution.resize(hgraph->GetNumVertices());
      for (int cluster_id = 0; cluster_id < coarse_hgraph->GetNumVertices();
           cluster_id++) {
        const int part_id = top_solution[cluster_id];
        for (const auto& v : coarse_hgraph->GetVertexCAttr(cluster_id)) {
          refined_solution[v] = part_id;
        }
      }
      top_solution = refined_solution;
    }

    // Parallel refine all the solutions
    std::vector<std::thread> threads;
    threads.reserve(top_solutions.size());
    for (auto& top_solution : top_solutions) {
      threads.emplace_back(&par::MultilevelPartitioner::CallRefiner,
                           this,
                           hgraph,
                           std::ref(upper_block_balance),
                           std::ref(lower_block_balance),
                           std::ref(top_solution));
    }
    for (auto& th : threads) {
      th.join();
    }
    threads.clear();

    // update the best_solution_id
    float best_cost = std::numeric_limits<float>::max();
    for (auto i = 0; i < top_solutions.size(); i++) {
      const float cost
          = evaluator_->CutEvaluator(hgraph, top_solutions[i], false).cost;
      if (best_cost > cost) {
        best_cost = cost;
        best_solution_id = i;
      }
    }
    debugPrint(logger_,
               PAR,
               "refinement",
               1,
               "Level {} :: num_vertices = {}, num_hyperedges = {},"
               " cutcost = {}, best_solution_id = {}",
               ++num_level,
               hgraph->GetNumVertices(),
               hgraph->GetNumHyperedges(),
               best_cost,
               best_solution_id);
  }
}

// Refine function
// k_way_pm_refinement,
// k_way_fm_refinement and greedy refinement
void MultilevelPartitioner::CallRefiner(
    const HGraphPtr& hgraph,
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance,
    std::vector<int>& solution) const
{
  if (num_parts_ > 1) {  // Pair-wise FM only used for multi-way partitioning
    k_way_pm_refiner_->Refine(
        hgraph, upper_block_balance, lower_block_balance, solution);
  }
  k_way_fm_refiner_->Refine(
      hgraph, upper_block_balance, lower_block_balance, solution);
  greedy_refiner_->Refine(
      hgraph, upper_block_balance, lower_block_balance, solution);
}

// Perform cut-overlay clustering and ILP-based partitioning
// The ILP-based partitioning uses top_solutions[best_solution_id] as a hint,
// such that the runtime can be signficantly reduced
std::vector<int> MultilevelPartitioner::CutOverlayILPPart(
    const HGraphPtr& hgraph,
    const Matrix<float>& upper_block_balance,
    const Matrix<float>& lower_block_balance,
    const Matrix<int>& top_solutions,
    int best_solution_id) const
{
  std::vector<int> optimal_solution = top_solutions[best_solution_id];
  std::vector<int> vertex_cluster_vec(hgraph->GetNumVertices(), -1);
  // check if the hyperedge is cut by solutions
  std::vector<bool> hyperedge_mask(hgraph->GetNumHyperedges(), false);
  for (const auto& solution : top_solutions) {
    for (int e = 0; e < hgraph->GetNumHyperedges(); e++) {
      if (hyperedge_mask[e] == true) {
        continue;  // This hyperedge has been cut
      }

      const auto range = hgraph->Vertices(e);
      const int block_id = solution[*range.begin()];
      for (const int vertex :
           boost::make_iterator_range(range.begin() + 1, range.end())) {
        if (solution[vertex] != block_id) {
          hyperedge_mask[e] = true;
          break;  // end this hyperedge
        }
      }
    }
  }

  // pre-order BFS to traverse the hypergraph
  auto lambda_detect_connected_components = [&](int v, int cluster_id) -> void {
    std::queue<int> wavefront;
    wavefront.push(v);
    while (wavefront.empty() == false) {
      const int u = wavefront.front();
      wavefront.pop();
      for (const int e : hgraph->Edges(u)) {
        if (hyperedge_mask[e] == true) {
          continue;  // this hyperedge has been cut
        }
        for (const int v_nbr : hgraph->Vertices(e)) {
          if (vertex_cluster_vec[v_nbr] == -1) {
            vertex_cluster_vec[v_nbr] = cluster_id;
            wavefront.push(v_nbr);
          }
        }
      }
    }
  };

  // detect the connected components and mask each connected component as a
  // cluster
  int cluster_id = -1;
  // map the initial optimal solution to the solution of clustered hgraph
  std::vector<int> init_solution;
  for (int v = 0; v < hgraph->GetNumVertices(); v++) {
    if (vertex_cluster_vec[v] == -1) {
      vertex_cluster_vec[v] = ++cluster_id;
      init_solution.push_back(top_solutions[best_solution_id][v]);
      lambda_detect_connected_components(v, cluster_id);
    }
  }

  const int num_clusters = cluster_id + 1;
  std::vector<std::vector<int>> cluster_attr;
  cluster_attr.reserve(num_clusters);
  for (int id = 0; id < num_clusters; id++) {
    std::vector<int> group_cluster{};
    cluster_attr.push_back(group_cluster);
  }
  for (int v = 0; v < hgraph->GetNumVertices(); v++) {
    cluster_attr[vertex_cluster_vec[v]].push_back(v);
  }

  // Call ILP-based partitioning
  HGraphPtr clustered_hgraph = coarsener_->GroupVertices(hgraph, cluster_attr);
  debugPrint(logger_,
             PAR,
             "cut_overlay_clustering",
             1,
             "num_vertices = {}, num_hyperedges = {}",
             clustered_hgraph->GetNumVertices(),
             clustered_hgraph->GetNumHyperedges());

  if (num_clusters <= num_vertices_threshold_ilp_) {
    partitioner_->Partition(clustered_hgraph,
                            upper_block_balance,
                            lower_block_balance,
                            init_solution,
                            PartitionType::INIT_DIRECT_ILP);
  } else {
    clustered_hgraph->SetCommunity(init_solution);
    init_solution = SingleCycleRefinement(
        clustered_hgraph, upper_block_balance, lower_block_balance);
  }

  // map the solution back to the original hypergraph
  for (int c_id = 0; c_id < clustered_hgraph->GetNumVertices(); c_id++) {
    const int block_id = init_solution[c_id];
    for (const auto& v : clustered_hgraph->GetVertexCAttr(c_id)) {
      optimal_solution[v] = block_id;
    }
  }

  debugPrint(logger_,
             PAR,
             "cut_overlay_clustering",
             1,
             "Statistics of cut-overlay solution:");
  evaluator_->CutEvaluator(hgraph, optimal_solution, false);
  return optimal_solution;
}

}  // namespace par
