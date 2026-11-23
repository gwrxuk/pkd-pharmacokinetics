#include "pkd/knowledge_base.hpp"

#include <stdexcept>

namespace pkd {

Direction parse_direction(const std::string& action) {
  std::string s = lower(action);
  if (s == "down" || s.find("inhib") != std::string::npos || s.find("antagon") != std::string::npos ||
      s.find("block") != std::string::npos || s.find("suppress") != std::string::npos ||
      s.find("inverse") != std::string::npos)
    return Direction::Down;
  if (s == "up" || s.find("agonist") != std::string::npos || s.find("activ") != std::string::npos ||
      s.find("stimul") != std::string::npos || s.find("opener") != std::string::npos ||
      s.find("enhanc") != std::string::npos || s.find("potentiat") != std::string::npos)
    return Direction::Up;
  return Direction::Any;
}

const char* to_string(Direction d) {
  switch (d) {
    case Direction::Down: return "down";
    case Direction::Up: return "up";
    default: return "any";
  }
}

namespace {
std::vector<TargetAction> parse_targets(const std::string& field) {
  std::vector<TargetAction> out;
  for (auto& item : split_nonempty(field, ';')) {
    auto parts = split(item, ':');
    TargetAction t;
    t.symbol = trim(parts[0]);
    t.action = parts.size() > 1 ? trim(parts[1]) : "";
    t.dir = parse_direction(t.action);
    out.push_back(t);
  }
  return out;
}
}  // namespace

KnowledgeBase KnowledgeBase::load(const std::string& dir) {
  KnowledgeBase kb;

  auto tt = read_csv(dir + "/targets.csv");
  for (size_t r = 0; r < tt.rows.size(); ++r) {
    TargetInfo t;
    t.symbol = tt.get(r, "symbol");
    t.synonyms = split_nonempty(tt.get(r, "synonyms"), '|');
    t.pathways = split_nonempty(tt.get(r, "pathways"), ';');
    kb.targets.push_back(t);
  }

  auto dt = read_csv(dir + "/diseases.csv");
  for (size_t r = 0; r < dt.rows.size(); ++r) {
    Disease d;
    d.key = dt.get(r, "key");
    d.area = dt.get(r, "area");
    d.synonyms = split_nonempty(dt.get(r, "synonyms"), '|');
    for (auto& item : split_nonempty(dt.get(r, "targets"), ';')) {
      auto p = split(item, ':');
      if (p.size() < 3) throw std::runtime_error("bad disease target '" + item + "' in " + d.key);
      DiseaseTarget x;
      x.symbol = trim(p[0]);
      x.dir = parse_direction(trim(p[1]));
      x.weight = to_double(p[2], 0.5);
      d.targets.push_back(x);
    }
    d.pathways = split_nonempty(dt.get(r, "pathways"), ';');
    kb.diseases.push_back(d);
  }

  auto rt = read_csv(dir + "/drugs.csv");
  for (size_t r = 0; r < rt.rows.size(); ++r) {
    DrugRecord d;
    d.name = rt.get(r, "name");
    d.mw = to_double(rt.get(r, "mw"), nan_v());
    d.dose_mg = to_double(rt.get(r, "dose_mg"), nan_v());
    d.route = rt.get(r, "route");
    d.tau_h = to_double(rt.get(r, "tau_h"), 24);
    d.F = to_double(rt.get(r, "F"), nan_v());
    d.ka = to_double(rt.get(r, "ka_h"), nan_v());
    d.CL = to_double(rt.get(r, "CL_L_h"), nan_v());
    d.V = to_double(rt.get(r, "V_L"), nan_v());
    d.fu = to_double(rt.get(r, "fu"), nan_v());
    d.ceff = to_double(rt.get(r, "ceff_mg_L"), nan_v());
    d.targets = parse_targets(rt.get(r, "targets"));
    d.indications = split_nonempty(rt.get(r, "indications"), ';');
    d.obs_cmax = to_double(rt.get(r, "obs_cmax_mg_L"), nan_v());
    d.obs_tmax = to_double(rt.get(r, "obs_tmax_h"), nan_v());
    d.obs_auc = to_double(rt.get(r, "obs_auc_mg_h_L"), nan_v());
    if (d.CL > 0 && d.V > 0) d.half_life_h = kLn2 * d.V / d.CL;
    d.source = "data/drugs.csv";
    kb.drugs.push_back(d);
  }

  try {
    auto ct = read_csv(dir + "/contraindications.csv");
    for (size_t r = 0; r < ct.rows.size(); ++r)
      kb.contraindications.push_back({ct.get(r, "drug"), ct.get(r, "disease"), ct.get(r, "reason")});
  } catch (const std::exception&) {
    // optional file
  }

  // sanity: every disease target should be known
  for (auto& d : kb.diseases)
    for (auto& t : d.targets)
      if (!kb.target(t.symbol)) throw std::runtime_error("disease " + d.key + " uses unknown target " + t.symbol);
  return kb;
}

const TargetInfo* KnowledgeBase::target(const std::string& symbol) const {
  for (auto& t : targets)
    if (t.symbol == symbol) return &t;
  return nullptr;
}

const Disease* KnowledgeBase::disease(const std::string& key) const {
  for (auto& d : diseases)
    if (d.key == key) return &d;
  return nullptr;
}

const DrugRecord* KnowledgeBase::drug(const std::string& name) const {
  std::string n = lower(name);
  for (auto& d : drugs)
    if (lower(d.name) == n) return &d;
  return nullptr;
}

}  // namespace pkd
