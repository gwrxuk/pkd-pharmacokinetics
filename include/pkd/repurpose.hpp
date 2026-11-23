#pragma once
#include <string>
#include <vector>

#include "pkd/knowledge_base.hpp"

namespace pkd {

struct ScoringWeights {
  double pathway = 0.3;      // weight of pathway (Jaccard) similarity
  double conflict = 1.0;     // penalty for acting on a target in the wrong direction
  bool use_exposure = true;  // scale by PK target coverage
};

struct Hypothesis {
  std::string disease;
  double score = 0;
  double target_score = 0;  // noisy-OR over matched disease targets (direction-consistent)
  double pathway_score = 0; // Jaccard(drug pathways, disease pathways)
  double coverage = 0;      // Css,avg / (Css,avg + Ceff); NaN if unknown
  double conflict = 0;      // noisy-OR over targets hit in the opposite direction
  bool known = false;       // already an indication of this drug
  bool safety_flag = false; // conflict > 0 -> potential contraindication
  std::vector<std::string> evidence;
};

// Steady-state average concentration under the labelled regimen, mg/L (NaN if unknown).
double css_avg(const DrugRecord& d);
// Fraction of "effect" reached at steady state: Css/(Css+Ceff). NaN if Ceff unknown.
double exposure_coverage(const DrugRecord& d);

// Score every disease in the knowledge base for this drug, best first.
std::vector<Hypothesis> rank_diseases(const DrugRecord& d, const KnowledgeBase& kb, const ScoringWeights& w = {});

}  // namespace pkd
