#pragma once
#include <string>
#include <vector>

#include "pkd/knowledge_base.hpp"
#include "pkd/repurpose.hpp"

namespace pkd {

struct RankingMetrics {
  double auroc = 0, mrr = 0, recall1 = 0, recall3 = 0, recall5 = 0;
  int n_pos = 0, n_pairs = 0;
};

// Rank every disease for every drug and compare with known indications.
// Ranks are "filtered": other known indications of the same drug are not
// counted against a positive. Ties are resolved pessimistically.
RankingMetrics evaluate_ranking(const std::vector<DrugRecord>& drugs, const KnowledgeBase& kb,
                                const ScoringWeights& w);

struct ValidationOptions {
  int permutations = 1000;
  unsigned seed = 42;
  int top_k = 3;
  std::string pdf_dir;    // folder with monograph PDFs for retrospective rediscovery
  std::string truth_csv;  // pdf,later_indications,note
};

// Prints a full validation report to stdout. Returns 0.
int run_validation(const KnowledgeBase& kb, const ValidationOptions& opt);

}  // namespace pkd
