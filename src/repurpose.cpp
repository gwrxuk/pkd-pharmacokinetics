#include "pkd/repurpose.hpp"

#include <algorithm>
#include <cmath>
#include <set>

namespace pkd {

double css_avg(const DrugRecord& d) {
  if (!(d.CL > 0) || !(d.dose_mg > 0) || !(d.tau_h > 0)) return nan_v();
  double F = std::isnan(d.F) ? 1.0 : d.F;
  return F * d.dose_mg / (d.CL * d.tau_h);
}

double exposure_coverage(const DrugRecord& d) {
  double css = css_avg(d);
  if (std::isnan(css) || !(d.ceff > 0)) return nan_v();
  return css / (css + d.ceff);
}

std::vector<Hypothesis> rank_diseases(const DrugRecord& d, const KnowledgeBase& kb, const ScoringWeights& w) {
  std::set<std::string> drug_paths;
  for (auto& t : d.targets)
    if (auto ti = kb.target(t.symbol)) drug_paths.insert(ti->pathways.begin(), ti->pathways.end());

  double cov = exposure_coverage(d);
  // Unknown coverage -> neutral factor 0.75 (same as coverage 0.5).
  double cov_factor = !w.use_exposure ? 1.0 : std::isnan(cov) ? 0.75 : 0.5 + 0.5 * cov;

  std::vector<Hypothesis> out;
  for (auto& dis : kb.diseases) {
    Hypothesis h;
    h.disease = dis.key;
    h.coverage = cov;
    double keep_agree = 1, keep_conf = 1;
    for (auto& need : dis.targets)
      for (auto& t : d.targets) {
        if (t.symbol != need.symbol) continue;
        if (t.dir == Direction::Any || need.dir == Direction::Any) {
          keep_agree *= 1 - 0.8 * need.weight;
          h.evidence.push_back(t.symbol + " (direction unknown) ~ disease needs " + to_string(need.dir) +
                               " (w=" + fmt(need.weight, 2) + ")");
        } else if (t.dir == need.dir) {
          keep_agree *= 1 - need.weight;
          h.evidence.push_back(t.symbol + " " + to_string(t.dir) + " matches disease need (w=" +
                               fmt(need.weight, 2) + ")");
        } else {
          keep_conf *= 1 - need.weight;
          h.evidence.push_back("CONFLICT: " + t.symbol + " " + to_string(t.dir) + " but disease needs " +
                               to_string(need.dir) + " (w=" + fmt(need.weight, 2) + ")");
        }
      }
    h.target_score = 1 - keep_agree;
    h.conflict = 1 - keep_conf;

    std::set<std::string> dp(dis.pathways.begin(), dis.pathways.end());
    size_t inter = 0;
    for (auto& p : drug_paths) inter += dp.count(p);
    size_t uni = drug_paths.size() + dp.size() - inter;
    h.pathway_score = uni ? (double)inter / uni : 0;
    if (inter) {
      std::vector<std::string> shared;
      for (auto& p : drug_paths)
        if (dp.count(p)) shared.push_back(p);
      h.evidence.push_back("shared pathways: " + join(shared, ", "));
    }

    h.score = (h.target_score + w.pathway * h.pathway_score) * cov_factor - w.conflict * h.conflict;
    h.known = std::find(d.indications.begin(), d.indications.end(), dis.key) != d.indications.end();
    h.safety_flag = h.conflict > 0;
    out.push_back(h);
  }
  std::stable_sort(out.begin(), out.end(), [](const Hypothesis& a, const Hypothesis& b) { return a.score > b.score; });
  return out;
}

}  // namespace pkd
