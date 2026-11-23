// Minimal self-contained unit tests (no framework dependency).
#include <cmath>
#include <cstdio>
#include <string>

#include "pkd/monograph.hpp"
#include "pkd/pdf_text.hpp"
#include "pkd/pk.hpp"
#include "pkd/repurpose.hpp"

using namespace pkd;

static int g_fail = 0, g_pass = 0;
#define CHECK(cond)                                                     \
  do {                                                                  \
    if (cond) {                                                         \
      ++g_pass;                                                         \
    } else {                                                            \
      ++g_fail;                                                         \
      std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);       \
    }                                                                   \
  } while (0)
#define CHECK_NEAR(a, b, rel) CHECK(std::abs((a) - (b)) <= (rel) * std::max(std::abs((double)(b)), 1e-12))

static void test_pk_vs_analytic() {
  PKModel m;
  m.F = 0.8; m.ka = 1.2; m.CL = 5; m.V = 50;
  for (double t : {0.5, 1.0, 2.0, 5.0, 10.0, 30.0}) {
    double num = 0;
    integrate(m, 100, 1e9, 1, t, [&](double tt, double c) { if (std::abs(tt - t) < 1e-9) num = c; });
    CHECK_NEAR(num, bateman(m, 100, t), 1e-5);
  }
  auto sd = single_dose_metrics(m, 100);
  CHECK_NEAR(sd.auc_inf, m.F * 100 / m.CL, 2e-3);  // AUC = F*D/CL
  double k = m.kel();
  CHECK_NEAR(sd.tmax, std::log(m.ka / k) / (m.ka - k), 0.02);
}

static void test_iv_two_compartment_auc() {
  PKModel m;
  m.iv = true; m.CL = 4; m.V = 20; m.Q = 6; m.V2 = 60;
  auto sd = single_dose_metrics(m, 200);
  CHECK_NEAR(sd.auc_inf, 200 / 4.0, 5e-3);
  CHECK(sd.half_life > kLn2 * 20 / 4);  // terminal phase slower than k10 alone
}

static void test_steady_state() {
  PKModel m;
  m.F = 0.5; m.ka = 1.5; m.CL = 10; m.V = 100;
  auto ss = steady_state(m, 200, 12);
  CHECK_NEAR(ss.cavg, 0.5 * 200 / (10 * 12), 5e-3);  // Css,avg = F*D/(CL*tau)
  double k = m.kel();
  CHECK_NEAR(ss.accumulation, 1.0 / (1 - std::exp(-k * 12)), 0.05);
}

static void test_ka_from_tmax() {
  for (double ka : {1.3, 0.05, 4.0}) {
    double k = 0.1;
    double tmax = std::log(ka / k) / (ka - k);
    CHECK_NEAR(ka_from_tmax(tmax, k), ka, 1e-6);
  }
}

static void test_content_stream() {
  std::string t = text_from_content_stream("BT /F1 10 Tf (Hel) Tj [(lo) -300 (World)] TJ 0 -14 Td (Next \\(x\\)) Tj ET");
  CHECK(t.find("Hello World") != std::string::npos);
  CHECK(t.find("\nNext (x)") != std::string::npos);
}

static void test_pdf_roundtrip() {
  std::vector<std::string> lines = {"Drug Name: Testamol", "Clearance: 600 mL/min (approx)", "Back\\slash ok"};
  for (bool z : {true, false}) {
    std::string pdf = build_text_pdf(lines, z);
    std::string txt = extract_text_builtin(pdf);
    CHECK(txt.find("Drug Name: Testamol") != std::string::npos);
    CHECK(txt.find("600 mL/min (approx)") != std::string::npos);
    CHECK(txt.find("Back\\slash ok") != std::string::npos);
  }
}

static void test_monograph_parse(const KnowledgeBase& kb) {
  std::string txt =
      "Drug Name: Examplolol\n"
      "Indications and Usage:\n"
      "Indicated for pulmonary arterial hypertension.\n\n"
      "Usual dose: 40 mg orally twice daily.\n"
      "Mechanism of Action:\n"
      "Examplolol is a nonselective beta-adrenergic receptor blocking agent.\n\n"
      "Bioavailability: 25%\n"
      "Volume of distribution: approximately 4 L/kg\n"
      "Clearance: 1,000 mL/min\n"
      "Time to peak: 1 to 2 hours\n"
      "Peak plasma concentration (Cmax): 35 ng/mL\n"
      "AUC: 170 ng*h/mL\n"
      "Protein binding: 90%\n";
  std::vector<std::string> warn;
  DrugRecord d = parse_monograph(txt, kb, warn);
  CHECK(d.name == "Examplolol");
  CHECK_NEAR(d.dose_mg, 40, 1e-9);
  CHECK_NEAR(d.tau_h, 12, 1e-9);
  CHECK_NEAR(d.F, 0.25, 1e-9);
  CHECK_NEAR(d.V, 280, 1e-9);
  CHECK_NEAR(d.CL, 60, 1e-9);
  CHECK_NEAR(d.obs_tmax, 1.5, 1e-9);
  CHECK_NEAR(d.obs_cmax, 0.035, 1e-9);
  CHECK_NEAR(d.obs_auc, 0.17, 1e-9);
  CHECK_NEAR(d.fu, 0.1, 1e-9);
  CHECK(!std::isnan(d.ka));
  bool b1 = false, b2 = false;
  for (auto& t : d.targets) {
    b1 |= t.symbol == "ADRB1" && t.dir == Direction::Down;
    b2 |= t.symbol == "ADRB2" && t.dir == Direction::Down;
  }
  CHECK(b1 && b2);
  // longest-match: PAH must not also be read as systemic hypertension
  CHECK(d.indications.size() == 1 && d.indications[0] == "pulmonary_arterial_hypertension");
}

static void test_direction_mining(const KnowledgeBase& kb) {
  auto t = mine_targets("Metformin activates AMP-activated protein kinase (AMPK) through inhibition of "
                        "respiratory chain complex I.", kb);
  bool ampk_up = false, c1_down = false;
  for (auto& x : t) {
    ampk_up |= x.symbol == "PRKAA1" && x.dir == Direction::Up;
    c1_down |= x.symbol == "MT_COMPLEX_I" && x.dir == Direction::Down;
  }
  CHECK(ampk_up);
  CHECK(c1_down);
}

static void test_repurposing(const KnowledgeBase& kb) {
  const DrugRecord* p = kb.drug("Propranolol");
  CHECK(p != nullptr);
  bool asthma_flagged = false;
  for (auto& h : rank_diseases(*p, kb))
    if (h.disease == "asthma") asthma_flagged = h.safety_flag && h.score < 0;
  CHECK(asthma_flagged);

  const DrugRecord* s = kb.drug("Sildenafil");
  auto hs = rank_diseases(*s, kb);
  CHECK(hs[0].disease == "erectile_dysfunction" || hs[0].disease == "pulmonary_arterial_hypertension");
  CHECK(hs[1].disease == "erectile_dysfunction" || hs[1].disease == "pulmonary_arterial_hypertension");
}

int main() {
  KnowledgeBase kb = KnowledgeBase::load(PKD_DATA_DIR);
  test_pk_vs_analytic();
  test_iv_two_compartment_auc();
  test_steady_state();
  test_ka_from_tmax();
  test_content_stream();
  test_pdf_roundtrip();
  test_monograph_parse(kb);
  test_direction_mining(kb);
  test_repurposing(kb);
  std::printf("%d passed, %d failed\n", g_pass, g_fail);
  return g_fail ? 1 : 0;
}
