#pragma once
#include <string>
#include <utility>
#include <vector>

#include "pkd/util.hpp"

namespace pkd {

enum class Direction { Down, Up, Any };  // inhibit/antagonize, activate/agonize, unknown/modulate
Direction parse_direction(const std::string& action);
const char* to_string(Direction d);

struct TargetAction {
  std::string symbol;
  Direction dir = Direction::Any;
  std::string action;  // free text, e.g. "inhibitor"
};

// A drug, either curated (data/drugs.csv) or extracted from a monograph PDF.
// Units: dose mg, time h, CL L/h, V L, concentrations mg/L, AUC mg*h/L.
struct DrugRecord {
  std::string name;
  double mw = nan_v();
  double dose_mg = nan_v();
  std::string route = "oral";
  double tau_h = nan_v();
  double F = nan_v();
  double ka = nan_v();
  double CL = nan_v();
  double V = nan_v();
  double Q = 0, V2 = 0;  // optional peripheral compartment
  double fu = nan_v();
  double ceff = nan_v();  // effective (therapeutic) plasma concentration, mg/L
  double half_life_h = nan_v();
  double obs_cmax = nan_v(), obs_tmax = nan_v(), obs_auc = nan_v();
  std::vector<TargetAction> targets;
  std::vector<std::string> indications;  // disease keys
  std::vector<std::string> notes;        // how values were derived
  std::string source;
};

struct DiseaseTarget {
  std::string symbol;
  Direction dir = Direction::Down;  // what a treatment should do to this target
  double weight = 0.5;              // 0..1 strength of target-disease evidence
};

struct Disease {
  std::string key;
  std::string area;
  std::vector<std::string> synonyms;
  std::vector<DiseaseTarget> targets;
  std::vector<std::string> pathways;
};

struct TargetInfo {
  std::string symbol;
  std::vector<std::string> synonyms;
  std::vector<std::string> pathways;
};

struct Contraindication {
  std::string drug, disease, reason;
};

class KnowledgeBase {
 public:
  static KnowledgeBase load(const std::string& dir);

  std::vector<DrugRecord> drugs;
  std::vector<Disease> diseases;
  std::vector<TargetInfo> targets;
  std::vector<Contraindication> contraindications;

  const TargetInfo* target(const std::string& symbol) const;
  const Disease* disease(const std::string& key) const;
  const DrugRecord* drug(const std::string& name) const;  // case-insensitive
};

}  // namespace pkd
