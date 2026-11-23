#include "pkd/validate.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>

#include "pkd/monograph.hpp"
#include "pkd/pk.hpp"

namespace pkd {

RankingMetrics evaluate_ranking(const std::vector<DrugRecord>& drugs, const KnowledgeBase& kb,
                                const ScoringWeights& w) {
  RankingMetrics m;
  std::vector<double> pos, neg;
  double rr = 0;
  int r1 = 0, r3 = 0, r5 = 0;
  for (auto& d : drugs) {
    if (d.indications.empty()) continue;
    auto hyps = rank_diseases(d, kb, w);
    std::vector<double> negs;
    for (auto& h : hyps) {
      (h.known ? pos : neg).push_back(h.score);
      if (!h.known) negs.push_back(h.score);
    }
    for (auto& h : hyps) {
      if (!h.known) continue;
      int rank = 1;
      for (double s : negs) rank += s >= h.score;  // pessimistic ties
      m.n_pos++;
      rr += 1.0 / rank;
      r1 += rank <= 1;
      r3 += rank <= 3;
      r5 += rank <= 5;
    }
  }
  double wins = 0;
  for (double p : pos)
    for (double n : neg) wins += p > n ? 1.0 : p == n ? 0.5 : 0.0;
  m.n_pairs = (int)(pos.size() + neg.size());
  if (!pos.empty() && !neg.empty()) m.auroc = wins / ((double)pos.size() * neg.size());
  if (m.n_pos) {
    m.mrr = rr / m.n_pos;
    m.recall1 = (double)r1 / m.n_pos;
    m.recall3 = (double)r3 / m.n_pos;
    m.recall5 = (double)r5 / m.n_pos;
  }
  return m;
}

namespace {

void rule(const char* title) {
  std::printf("\n==================================================================================\n");
  std::printf(" %s\n", title);
  std::printf("==================================================================================\n");
}

double fold(double pred, double obs) { return (pred > 0 && obs > 0) ? pred / obs : nan_v(); }

void print_metrics(const char* label, const RankingMetrics& m) {
  std::printf("  %-34s AUROC %.3f  MRR %.3f  R@1 %.2f  R@3 %.2f  R@5 %.2f\n", label, m.auroc, m.mrr, m.recall1,
              m.recall3, m.recall5);
}

// ---------------------------------------------------------------- PK
void validate_pk(const KnowledgeBase& kb) {
  rule("1. PK model vs published single-dose values (curated drugs, approximate literature data)");
  std::printf("  Acceptance: prediction within 2-fold of observed (common PBPK/PopPK criterion)\n\n");
  std::printf("  %-13s %6s | %9s %9s %6s | %6s %6s | %9s %9s %6s\n", "drug", "dose", "Cmax_pred", "Cmax_obs",
              "fold", "Tmax_p", "Tmax_o", "AUC_pred", "AUC_obs", "fold");
  double sum_log = 0;
  int n = 0, pass = 0;
  for (auto& d : kb.drugs) {
    PKModel m;
    if (!model_from(d, m) || !(d.dose_mg > 0)) continue;
    auto sd = single_dose_metrics(m, d.dose_mg);
    double fc = fold(sd.cmax, d.obs_cmax), fa = fold(sd.auc_inf, d.obs_auc);
    std::printf("  %-13s %4gmg | %9s %9s %6s | %6s %6s | %9s %9s %6s\n", d.name.c_str(), d.dose_mg,
                fmt(sd.cmax).c_str(), fmt(d.obs_cmax).c_str(), fmt(fc, 2).c_str(), fmt(sd.tmax, 2).c_str(),
                fmt(d.obs_tmax, 2).c_str(), fmt(sd.auc_inf).c_str(), fmt(d.obs_auc).c_str(), fmt(fa, 2).c_str());
    if (!std::isnan(fc)) {
      sum_log += std::abs(std::log10(fc));
      n++;
      pass += fc >= 0.5 && fc <= 2.0;
    }
  }
  if (n)
    std::printf("\n  Cmax: %d/%d within 2-fold, geometric mean fold error (GMFE) = %.2f\n", pass, n,
                std::pow(10.0, sum_log / n));
  std::printf("  Note: AUC = F*Dose/CL, and CL is itself usually derived from AUC in the literature, so\n"
              "  AUC agreement mainly checks the integrator; Cmax/Tmax (depend on ka and V) are the real test.\n");
}

// ---------------------------------------------------------------- repurposing
void validate_ranking(const KnowledgeBase& kb, const ValidationOptions& opt) {
  rule("2. Indication recovery on curated drugs (does the score rank known indications first?)");
  std::printf("  %zu drugs x %zu diseases. Positives = labelled/established indications.\n\n", kb.drugs.size(),
              kb.diseases.size());
  std::printf("  %-13s %-58s %s\n", "drug", "known indications (filtered rank)", "top non-indication");
  for (auto& d : kb.drugs) {
    auto hyps = rank_diseases(d, kb);
    std::vector<double> negs;
    for (auto& h : hyps)
      if (!h.known) negs.push_back(h.score);
    std::string ranks;
    for (auto& h : hyps) {
      if (!h.known) continue;
      int rank = 1;
      for (double s : negs) rank += s >= h.score;
      ranks += h.disease + "(" + std::to_string(rank) + ") ";
    }
    std::string top_new = "-";
    for (auto& h : hyps)
      if (!h.known && h.score > 0) {
        top_new = h.disease + " " + fmt(h.score, 2);
        break;
      }
    if (ranks.size() > 58) ranks = ranks.substr(0, 55) + "...";
    std::printf("  %-13s %-58s %s\n", d.name.c_str(), ranks.c_str(), top_new.c_str());
  }

  RankingMetrics full = evaluate_ranking(kb.drugs, kb, {});
  ScoringWeights tonly;
  tonly.pathway = 0;
  tonly.use_exposure = false;
  ScoringWeights noexp;
  noexp.use_exposure = false;
  std::printf("\n  Ablation:\n");
  print_metrics("targets only", evaluate_ranking(kb.drugs, kb, tonly));
  print_metrics("targets + pathways", evaluate_ranking(kb.drugs, kb, noexp));
  print_metrics("targets + pathways + PK exposure", full);

  // permutation test: shuffle which drug owns which target set
  std::mt19937 rng(opt.seed);
  std::vector<std::vector<TargetAction>> sets;
  for (auto& d : kb.drugs) sets.push_back(d.targets);
  int ge = 0;
  double sum = 0;
  for (int p = 0; p < opt.permutations; ++p) {
    std::shuffle(sets.begin(), sets.end(), rng);
    std::vector<DrugRecord> perm = kb.drugs;
    for (size_t i = 0; i < perm.size(); ++i) perm[i].targets = sets[i];
    double a = evaluate_ranking(perm, kb, {}).auroc;
    sum += a;
    ge += a >= full.auroc;
  }
  if (opt.permutations > 0) {
    std::printf("\n  Permutation test (%d shuffles of drug->target assignment, seed %u):\n", opt.permutations,
                opt.seed);
    std::printf("    null mean AUROC %.3f vs observed %.3f, empirical p = %.4f\n", sum / opt.permutations,
                full.auroc, (ge + 1.0) / (opt.permutations + 1.0));
  }
  std::printf("\n  CAVEAT: disease->target annotations were written with knowledge of how these drugs work,\n"
              "  so this is an internal-consistency check (upper bound), not proof of predictive power.\n"
              "  For an independent test, replace data/diseases.csv with associations from an external\n"
              "  source (e.g. Open Targets / DisGeNET) that was not built from drug indications.\n");
}

// ---------------------------------------------------------------- safety
void validate_safety(const KnowledgeBase& kb) {
  rule("3. Known contraindications (should be flagged as direction conflicts)");
  int hit = 0, n = 0;
  for (auto& c : kb.contraindications) {
    const DrugRecord* d = kb.drug(c.drug);
    if (!d) continue;
    n++;
    bool flagged = false;
    double score = 0;
    for (auto& h : rank_diseases(*d, kb))
      if (h.disease == c.disease) {
        flagged = h.safety_flag;
        score = h.score;
      }
    hit += flagged;
    std::printf("  %-13s x %-22s score %+6.2f  %s   (%s)\n", c.drug.c_str(), c.disease.c_str(), score,
                flagged ? "FLAGGED" : "missed ", c.reason.c_str());
  }
  std::printf("  %d/%d known contraindications flagged\n", hit, n);
}

// ---------------------------------------------------------------- PDF rediscovery
void validate_rediscovery(const KnowledgeBase& kb, const ValidationOptions& opt) {
  if (opt.truth_csv.empty()) return;
  CsvTable t;
  try {
    t = read_csv(opt.truth_csv);
  } catch (const std::exception& e) {
    std::printf("\n  (skipping rediscovery: %s)\n", e.what());
    return;
  }
  rule("4. Retrospective rediscovery from monograph PDFs (original label only -> later uses?)");
  std::printf("  Each PDF lists only the drug's ORIGINAL indication(s). A hit means a later-approved or\n"
              "  established use appears in the top-%d NEW hypotheses (non-label, no safety conflict).\n\n",
              opt.top_k);
  int hits = 0, total = 0;
  for (size_t r = 0; r < t.rows.size(); ++r) {
    std::string pdf = opt.pdf_dir + "/" + t.get(r, "pdf");
    std::vector<std::string> warn;
    DrugRecord d;
    std::string method;
    try {
      d = ingest_document(pdf, kb, warn, &method);
    } catch (const std::exception& e) {
      std::printf("  %s: ERROR %s\n", pdf.c_str(), e.what());
      continue;
    }
    auto hyps = rank_diseases(d, kb);
    std::vector<std::string> fresh;
    for (auto& h : hyps)
      if (!h.known && !h.safety_flag && h.score > 0) fresh.push_back(h.disease);

    PKModel m;
    std::string pk = "PK n/a";
    if (model_from(d, m) && d.dose_mg > 0) {
      auto sd = single_dose_metrics(m, d.dose_mg);
      pk = "Cmax pred " + fmt(sd.cmax) + " mg/L";
      if (d.obs_cmax > 0) pk += " vs label " + fmt(d.obs_cmax) + " (fold " + fmt(sd.cmax / d.obs_cmax, 2) + ")";
    }
    std::printf("  %s  [text: %s]  label: %s\n", d.name.c_str(), method.c_str(), join(d.indications, ", ").c_str());
    std::printf("    %s\n", pk.c_str());
    std::string top;
    for (size_t i = 0; i < fresh.size() && (int)i < opt.top_k; ++i) top += (i ? ", " : "") + fresh[i];
    std::printf("    top-%d new hypotheses: %s\n", opt.top_k, top.empty() ? "(none)" : top.c_str());
    for (auto& later : split_nonempty(t.get(r, "later_indications"), ';')) {
      auto it = std::find(fresh.begin(), fresh.end(), later);
      int rank = it == fresh.end() ? -1 : (int)(it - fresh.begin()) + 1;
      bool hit = rank > 0 && rank <= opt.top_k;
      hits += hit;
      total++;
      std::printf("    later use %-32s -> %s\n", later.c_str(),
                  rank < 0 ? "not proposed" : ("rank " + std::to_string(rank) + (hit ? "  HIT" : "")).c_str());
    }
  }
  if (total)
    std::printf("\n  Rediscovery hit rate @%d: %d/%d = %.0f%%\n", opt.top_k, hits, total, 100.0 * hits / total);
}

}  // namespace

int run_validation(const KnowledgeBase& kb, const ValidationOptions& opt) {
  std::printf("PKD validation report  (research use only - not clinical advice)\n");
  validate_pk(kb);
  validate_ranking(kb, opt);
  validate_safety(kb);
  validate_rediscovery(kb, opt);
  std::printf("\n");
  return 0;
}

}  // namespace pkd
