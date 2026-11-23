// pkd - pharmacokinetics (藥物動力學) + drug-repurposing research tool.
// Research / educational use only. Not clinical advice.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "pkd/monograph.hpp"
#include "pkd/pdf_text.hpp"
#include "pkd/pk.hpp"
#include "pkd/repurpose.hpp"
#include "pkd/validate.hpp"

using namespace pkd;

namespace {

struct Args {
  std::vector<std::string> pos;
  std::map<std::string, std::string> opt;
  std::string get(const std::string& k, const std::string& def = "") const {
    auto it = opt.find(k);
    return it == opt.end() ? def : it->second;
  }
  double num(const std::string& k, double def) const { return to_double(get(k), def); }
};

Args parse_args(int argc, char** argv, int from) {
  Args a;
  for (int i = from; i < argc; ++i) {
    std::string s = argv[i];
    if (s.rfind("--", 0) == 0) {
      std::string key = s.substr(2);
      if (i + 1 < argc && std::string(argv[i + 1]).rfind("--", 0) != 0)
        a.opt[key] = argv[++i];
      else
        a.opt[key] = "1";
    } else {
      a.pos.push_back(s);
    }
  }
  return a;
}

void usage() {
  std::printf(
      "pkd - pharmacokinetic modelling + drug repurposing from monograph PDFs\n"
      "      (research/education only, not clinical advice)\n\n"
      "usage:\n"
      "  pkd discover <file.pdf|.txt>... [--top N] [--csv DIR]   ingest monographs, simulate PK, rank new uses\n"
      "  pkd extract  <file.pdf>                                 print extracted text\n"
      "  pkd simulate <drug> [--dose mg] [--tau h] [--doses n] [--hours T] [--csv out.csv]\n"
      "  pkd rank     <drug> [--top N]                           repurposing ranks for a curated drug\n"
      "  pkd validate [--perm N] [--seed S] [--top K] [--pdf-dir DIR] [--truth CSV]\n"
      "  pkd list                                                curated drugs and diseases\n\n"
      "global: --data DIR   knowledge-base folder (default: %s)\n",
      PKD_DATA_DIR);
}

void ascii_plot(const std::vector<Sample>& s, int width = 64, int height = 12) {
  if (s.empty()) return;
  double tmax = s.back().t, cmax = 0;
  for (auto& x : s) cmax = std::max(cmax, x.c);
  if (cmax <= 0) return;
  std::vector<std::string> grid(height, std::string(width, ' '));
  for (auto& x : s) {
    int col = std::min(width - 1, (int)(x.t / tmax * (width - 1)));
    int row = std::min(height - 1, (int)((1 - x.c / cmax) * (height - 1) + 0.5));
    grid[row][col] = '*';
  }
  std::printf("\n  C (mg/L)\n");
  for (int r = 0; r < height; ++r) {
    double level = cmax * (1 - (double)r / (height - 1));
    std::printf("  %9s |%s\n", (r % 3 == 0) ? fmt(level, 3).c_str() : "", grid[r].c_str());
  }
  std::printf("  %9s +%s\n", "", std::string(width, '-').c_str());
  std::printf("  %9s 0%*s%s h\n", "", width - 2, "", fmt(tmax, 3).c_str());
}

void print_pk(const DrugRecord& d, bool show_reported) {
  PKModel m;
  std::string why;
  if (!model_from(d, m, &why) || !(d.dose_mg > 0)) {
    std::printf("  PK: unavailable (%s)\n", why.empty() ? "dose missing" : why.c_str());
    return;
  }
  auto sd = single_dose_metrics(m, d.dose_mg);
  std::printf("  PK single %g mg: Cmax %s mg/L at %s h | AUC0-inf %s mg*h/L | terminal t1/2 %s h\n", d.dose_mg,
              fmt(sd.cmax).c_str(), fmt(sd.tmax, 2).c_str(), fmt(sd.auc_inf).c_str(), fmt(sd.half_life).c_str());
  if (show_reported && (d.obs_cmax > 0 || d.obs_auc > 0)) {
    std::printf("     reported:");
    if (d.obs_cmax > 0) std::printf(" Cmax %s (pred/obs %s)", fmt(d.obs_cmax).c_str(), fmt(sd.cmax / d.obs_cmax, 2).c_str());
    if (d.obs_auc > 0) std::printf(" AUC %s (pred/obs %s)", fmt(d.obs_auc).c_str(), fmt(sd.auc_inf / d.obs_auc, 2).c_str());
    std::printf("\n");
  }
  auto ss = steady_state(m, d.dose_mg, d.tau_h);
  std::printf("  PK steady state q%gh: Cmax,ss %s | Cmin,ss %s | Cavg %s mg/L | accumulation %sx | ~%d doses to 90%%\n",
              d.tau_h, fmt(ss.cmax).c_str(), fmt(ss.cmin).c_str(), fmt(ss.cavg).c_str(),
              fmt(ss.accumulation, 2).c_str(), ss.doses_to_90pct);
  double cov = exposure_coverage(d);
  if (std::isnan(cov))
    std::printf("  Exposure coverage: unknown (no effective concentration) -> neutral weight\n");
  else
    std::printf("  Exposure coverage Css/(Css+Ceff): %.2f  (Ceff %s mg/L)\n", cov, fmt(d.ceff).c_str());
}

void print_hypotheses(const DrugRecord& d, const KnowledgeBase& kb, int top) {
  auto hyps = rank_diseases(d, kb);
  std::printf("\n  %-4s %-32s %6s %6s %6s  %-8s %s\n", "rank", "disease", "score", "target", "path", "status",
              "top evidence");
  int shown = 0;
  for (size_t i = 0; i < hyps.size() && shown < top; ++i) {
    auto& h = hyps[i];
    if (h.score <= 0) break;
    std::string status = h.known ? "label" : h.safety_flag ? "CAUTION" : "NEW";
    std::printf("  %-4zu %-32s %6.3f %6.2f %6.2f  %-8s %s\n", i + 1, h.disease.c_str(), h.score, h.target_score,
                h.pathway_score, status.c_str(), h.evidence.empty() ? "" : h.evidence[0].c_str());
    shown++;
  }
  for (auto& h : hyps) {
    if (!h.safety_flag || h.score > 0) continue;
    std::string why;
    for (auto& e : h.evidence)
      if (e.rfind("CONFLICT", 0) == 0) why = e;
    std::printf("  !    %-32s %6.3f  SAFETY %s\n", h.disease.c_str(), h.score, why.c_str());
  }
}

void describe(const DrugRecord& d) {
  std::printf("  Extracted: MW %s | dose %s mg %s q%gh | F %s | CL %s L/h | V %s L | ka %s /h | fu %s\n",
              fmt(d.mw, 4).c_str(), fmt(d.dose_mg).c_str(), d.route.c_str(), d.tau_h, fmt(d.F, 2).c_str(),
              fmt(d.CL).c_str(), fmt(d.V).c_str(), fmt(d.ka).c_str(), fmt(d.fu, 2).c_str());
  std::string ts;
  for (auto& t : d.targets) ts += t.symbol + "(" + to_string(t.dir) + ") ";
  std::printf("  Targets:   %s\n", ts.empty() ? "(none recognised)" : ts.c_str());
  std::printf("  Label indications: %s\n", d.indications.empty() ? "(none mapped)" : join(d.indications, ", ").c_str());
  for (auto& n : d.notes) std::printf("  note: %s\n", n.c_str());
}

int cmd_discover(const Args& a, const KnowledgeBase& kb) {
  if (a.pos.empty()) {
    usage();
    return 2;
  }
  int top = (int)a.num("top", 8);
  int rc = 0;
  for (auto& path : a.pos) {
    std::vector<std::string> warn;
    std::string method;
    DrugRecord d;
    try {
      d = ingest_document(path, kb, warn, &method);
    } catch (const std::exception& e) {
      std::fprintf(stderr, "error: %s\n", e.what());
      rc = 1;
      continue;
    }
    std::printf("\n=== %s   [%s, text via %s]\n", d.name.c_str(), path.c_str(), method.c_str());
    describe(d);
    for (auto& w : warn) std::printf("  warning: %s\n", w.c_str());
    print_pk(d, true);
    print_hypotheses(d, kb, top);
    if (!a.get("csv").empty()) {
      PKModel m;
      if (model_from(d, m) && d.dose_mg > 0) {
        std::string out = a.get("csv") + "/" + normalize_key(d.name) + "_pk.csv";
        std::ofstream f(out);
        f << "t_h,conc_mg_L\n";
        double T = std::max(5 * d.tau_h, 5 * m.terminal_half_life());
        for (auto& s : simulate(m, d.dose_mg, d.tau_h, (int)std::ceil(T / d.tau_h), T, 0.25))
          f << s.t << "," << s.c << "\n";
        std::printf("  wrote %s\n", out.c_str());
      }
    }
  }
  std::printf("\nHypotheses are computational research leads only - not clinical advice.\n");
  return rc;
}

int cmd_simulate(const Args& a, const KnowledgeBase& kb) {
  if (a.pos.empty()) {
    usage();
    return 2;
  }
  const DrugRecord* src = kb.drug(a.pos[0]);
  if (!src) {
    std::fprintf(stderr, "unknown drug '%s' (see `pkd list`)\n", a.pos[0].c_str());
    return 1;
  }
  DrugRecord d = *src;
  d.dose_mg = a.num("dose", d.dose_mg);
  d.tau_h = a.num("tau", d.tau_h);
  PKModel m;
  model_from(d, m);
  int n = (int)a.num("doses", 1);
  double T = a.num("hours", n > 1 ? n * d.tau_h : std::max(24.0, 5 * m.terminal_half_life()));
  std::printf("%s: %g mg %s every %g h x %d, model F=%s ka=%s/h CL=%s L/h V=%s L\n", d.name.c_str(), d.dose_mg,
              d.route.c_str(), d.tau_h, n, fmt(m.F).c_str(), fmt(m.ka).c_str(), fmt(m.CL).c_str(), fmt(m.V).c_str());
  print_pk(d, true);
  auto s = simulate(m, d.dose_mg, d.tau_h, n, T, T / 400);
  ascii_plot(s);
  if (!a.get("csv").empty()) {
    std::ofstream f(a.get("csv"));
    f << "t_h,conc_mg_L\n";
    for (auto& x : s) f << x.t << "," << x.c << "\n";
    std::printf("wrote %s (%zu rows)\n", a.get("csv").c_str(), s.size());
  }
  return 0;
}

int cmd_rank(const Args& a, const KnowledgeBase& kb) {
  if (a.pos.empty()) {
    usage();
    return 2;
  }
  const DrugRecord* d = kb.drug(a.pos[0]);
  if (!d) {
    std::fprintf(stderr, "unknown drug '%s'\n", a.pos[0].c_str());
    return 1;
  }
  std::printf("=== %s (curated)\n", d->name.c_str());
  describe(*d);
  print_pk(*d, true);
  print_hypotheses(*d, kb, (int)a.num("top", 8));
  return 0;
}

int cmd_list(const KnowledgeBase& kb) {
  std::printf("Curated drugs (%zu):\n", kb.drugs.size());
  for (auto& d : kb.drugs) {
    std::string ts;
    for (auto& t : d.targets) ts += t.symbol + " ";
    std::printf("  %-13s targets: %-28s indications: %s\n", d.name.c_str(), ts.c_str(), join(d.indications, ", ").c_str());
  }
  std::printf("\nDiseases (%zu):\n", kb.diseases.size());
  for (auto& d : kb.diseases) {
    std::string ts;
    for (auto& t : d.targets) ts += t.symbol + " ";
    std::printf("  %-32s %-16s %s\n", d.key.c_str(), d.area.c_str(), ts.c_str());
  }
  std::printf("\nTargets: %zu\n", kb.targets.size());
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    usage();
    return 2;
  }
  std::string cmd = argv[1];
  Args a = parse_args(argc, argv, 2);
  if (cmd == "-h" || cmd == "--help" || cmd == "help") {
    usage();
    return 0;
  }
  if (cmd == "extract") {
    if (a.pos.empty()) {
      usage();
      return 2;
    }
    auto r = extract_pdf_text(a.pos[0]);
    if (!r.ok) {
      std::fprintf(stderr, "error: %s\n", r.error.c_str());
      return 1;
    }
    std::fprintf(stderr, "[method: %s]\n", r.method.c_str());
    std::printf("%s\n", r.text.c_str());
    return 0;
  }

  KnowledgeBase kb;
  try {
    kb = KnowledgeBase::load(a.get("data", PKD_DATA_DIR));
  } catch (const std::exception& e) {
    std::fprintf(stderr, "error loading knowledge base: %s\n", e.what());
    return 1;
  }
  if (cmd == "discover" || cmd == "ingest") return cmd_discover(a, kb);
  if (cmd == "simulate") return cmd_simulate(a, kb);
  if (cmd == "rank") return cmd_rank(a, kb);
  if (cmd == "list") return cmd_list(kb);
  if (cmd == "validate") {
    ValidationOptions o;
    o.permutations = (int)a.num("perm", 1000);
    o.seed = (unsigned)a.num("seed", 42);
    o.top_k = (int)a.num("top", 3);
    o.pdf_dir = a.get("pdf-dir", std::string(PKD_SAMPLES_DIR) + "/pdfs");
    o.truth_csv = a.get("truth", std::string(PKD_SAMPLES_DIR) + "/rediscovery_truth.csv");
    return run_validation(kb, o);
  }
  usage();
  return 2;
}
