// Headers from belief propagation algorithm
#include <dai/alldai.h>
#include <dai/jtree.h>
#include <dai/bp.h>
#include <dai/decmap.h>

#include <fstream>
#include <iostream>
#include <string>

class Node_GT {
public:
  int id, gt;
  Node_GT(int node_id, int node_gt){
    id = node_id;
    gt = node_gt;
  }
};

class Node_Pair_GT {
public:
  int id_1, id_2;
  int gt_a, gt_b;
  Node_Pair_GT(int node_id_1, int node_id_2, int node_gt_a, int node_gt_b){
    id_1 = node_id_1; id_2 = node_id_2;
    gt_a = node_gt_a; gt_b = node_gt_b;
  }
};

void read_node_genotypes(std::string input_file, std::vector<Node_GT>& node_gts){
  int node_id, genotype;
  std::ifstream input(input_file);
  while (input >> node_id >> genotype)
    node_gts.push_back(Node_GT(node_id, genotype));
  input.close();
}

void read_pair_genotypes(std::string input_file, std::vector<Node_Pair_GT>& pair_gts){
  int node_id_1, node_id_2;
  int node_gt_a, node_gt_b;
  std::ifstream input(input_file);
  while (input >> node_id_1 >> node_id_2 >> node_gt_a >> node_gt_b)
    pair_gts.push_back(Node_Pair_GT(node_id_1, node_id_2, node_gt_a, node_gt_b));
  input.close();
}

void read_node_ids(std::string input_file, std::vector<int>& node_ids){
  int node_id;
  std::ifstream input(input_file);
  while (input >> node_id)
    node_ids.push_back(node_id);
  input.close();
}

void fix_genotypes(std::vector< std::pair<int,int> >& gt_info, dai::FactorGraph& fg){
  for (int i = 0; i < gt_info.size(); i++)
    fg.clamp(gt_info[i].first, gt_info[i].second, true);
}

int main(int argc, char *argv[]) {
  dai::FactorGraph fg;
  std::cerr << "Reading factor graph" << std::endl;
  fg.ReadFromFile(argv[1]);

  size_t maxiter    = 2000;
  dai::Real tol     = 1e-9;
  size_t verb       = 3;
  dai::Real damping = 0.01;

  dai::PropertySet opts;
  opts.set("maxiter", maxiter);   // Maximum number of iterations
  opts.set("tol", tol);           // Tolerance for convergence
  opts.set("verbose", verb);      // Verbosity (amount of output generated)
  opts.set("damping", damping);   // Degree to which messages should be damped 

  std::vector<Node_Pair_GT> pair_gts;
  read_pair_genotypes(argv[2], pair_gts);

  std::vector<int> node_ids;
  read_node_ids(argv[3], node_ids);
  
  dai::BP bp(fg, opts("updates", std::string("SEQRND"))("logdomain",true));
  bp.init();
  bp.run();


  for (unsigned int i = 0; i < pair_gts.size(); ++i){
    std::cout << pair_gts[i].id_1 << "\t" << pair_gts[i].id_2 << "\t"
	      << pair_gts[i].gt_a << "\t" << pair_gts[i].gt_b << "\t" 
	      << bp.belief(bp.var(pair_gts[i].id_1))[pair_gts[i].gt_a] << "\t" 
	      << bp.belief(bp.var(pair_gts[i].id_1))[pair_gts[i].gt_b] << std::endl;
  }
  
  for (unsigned int i = 0; i < node_ids.size(); ++i){
    std::cout << node_ids[i] << bp.belief(bp.var(node_ids[i])) << std::endl;
  }

  
  return 0;
}
