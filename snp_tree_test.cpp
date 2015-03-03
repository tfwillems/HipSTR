#include <iostream>
#include <thread>
#include <chrono>
#include <random>
#include <time.h>
#include <assert.h>

#include "snp_tree.h"

uint32_t randKey(uint32_t floor, uint32_t ceiling) {
  uint32_t range = ceiling - floor;
  return floor + range * ((double) rand() / (double) (RAND_MAX + 1.0));
}

std::pair<uint32_t, uint32_t> randomInterval(uint32_t maxStart, uint32_t maxLength, uint32_t maxStop){
  uint32_t start = randKey(0, maxStart);
  uint32_t stop  = std::min(randKey(start, start + maxLength), maxStop);
  return std::pair<uint32_t, uint32_t>(start, stop);
}

SNP randomSNP(uint32_t maxCoord){
  uint32_t coord = randKey(0, maxCoord);
  return SNP(coord, 'A', 'C');
}

int main() {  
  // a simple sanity check
  std::vector<SNP> sanity_snps;
  sanity_snps.push_back(SNP(30, 'A', 'C'));
  sanity_snps.push_back(SNP(50, 'A', 'C'));
  SNPTree sanity_tree(sanity_snps);
  std::vector<SNP> sanity_results;
  sanity_tree.findContained(35, 60, sanity_results);
  assert(sanity_results.size() == 1);
  sanity_results.clear();
  sanity_tree.findContained(15, 32, sanity_results);
  assert(sanity_results.size() == 1);
  sanity_results.clear();
  sanity_tree.findContained(15, 60, sanity_results);
  assert(sanity_results.size() == 2);
  sanity_results.clear();
  sanity_tree.findContained(15, 20, sanity_results);
  assert(sanity_results.size() == 0);

  srand((unsigned)time(NULL));
  std::vector<SNP> snps;
  std::vector< std::pair<uint32_t, uint32_t> > queries;
    
  // generate a test set of target intervals
  for (int i = 0; i < 100000; ++i)
    snps.push_back(randomSNP(10000000));

  // and queries
  for (int i = 0; i < 10000; ++i)
    queries.push_back(randomInterval(10000000, 1000, 10000000 + 1));

  typedef std::chrono::high_resolution_clock Clock;
  typedef std::chrono::milliseconds milliseconds;

  // using brute-force search
  std::vector<size_t> bruteforcecounts;
  Clock::time_point t0 = Clock::now();
  for (auto query_iter = queries.begin(); query_iter != queries.end(); ++query_iter) {
    std::vector<SNP> results;
    for (auto snp_iter = snps.begin(); snp_iter != snps.end(); ++snp_iter) {
      if (snp_iter->pos() >= query_iter->first && snp_iter->pos() <= query_iter->second)
	results.push_back(*snp_iter);
    }
    bruteforcecounts.push_back(results.size());
  }

  Clock::time_point t1 = Clock::now();
  milliseconds ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
  std::cout << "brute force:\t" << ms.count() << "ms" << std::endl;

  // using the SNP tree
  SNPTree tree = SNPTree(snps);
  std::vector<size_t> treecounts;
  t0 = Clock::now();
  for (auto query_iter = queries.begin(); query_iter != queries.end(); ++query_iter) {
    std::vector<SNP> results;
    tree.findContained(query_iter->first, query_iter->second, results);
    treecounts.push_back(results.size());
  }
  t1 = Clock::now();
  ms = std::chrono::duration_cast<milliseconds>(t1 - t0);
  std::cout << "interval tree:\t" << ms.count() << "ms" << std::endl;
  
  // check that the same number of results are returned
  auto bfc_iter = bruteforcecounts.begin();
  for (auto tree_iter = treecounts.begin(); tree_iter != treecounts.end(); ++tree_iter, ++bfc_iter)
    assert(*bfc_iter == *tree_iter);

  return 0;
}

