#include "pkd/monograph.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <optional>
#include <regex>
#include <set>
#include <stdexcept>

#include "pkd/pdf_text.hpp"
#include "pkd/pk.hpp"

namespace pkd {
namespace {

// number, optionally a range "3 to 5" / "3-5" (range -> midpoint unless lower bound requested)
const std::string NUM = R"(([0-9]+(?:\.[0-9]+)?)(?:\s*(?:-|to|–|—)\s*([0-9]+(?:\.[0-9]+)?))?)";

struct Doc {
  std::vector<std::string> raw, low;
};

Doc make_doc(const std::string& text) {
  static const std::regex thousands(R"((\d),(\d{3})(?!\d))");
  std::string t = std::regex_replace(text, thousands, "$1$2");
  Doc d;
  for (auto& l : split(t, '\n')) {
    std::string r = trim(l);
    d.raw.push_back(r);
    d.low.push_back(lower(r));
  }
  return d;
}

bool is_label_line(const std::string& raw) {
  static const std::regex re(R"(^[A-Za-z][A-Za-z0-9 /()\-.,&]{1,48}:)");
  return std::regex_search(raw, re);
}

std::string label_part(const Doc& d, size_t i) {
  if (!is_label_line(d.raw[i])) return "";
  return d.low[i].substr(0, d.low[i].find(':'));
}

// Lines that may hold a field. Pass 1: label appears before the colon of a
// "Label: value" line. Pass 2: label appears anywhere.
std::vector<size_t> candidates(const Doc& d, const std::vector<std::string>& labels,
                               const std::vector<std::string>& exclude = {}, bool label_only = false) {
  std::vector<size_t> out;
  std::set<size_t> seen;
  auto excluded = [&](size_t i) {
    for (auto& e : exclude)
      if (d.low[i].find(e) != std::string::npos) return true;
    return false;
  };
  for (int pass = 0; pass < (label_only ? 1 : 2); ++pass)
    for (auto& lab : labels)
      for (size_t i = 0; i < d.low.size(); ++i) {
        if (seen.count(i) || excluded(i)) continue;
        std::string hay = pass == 0 ? label_part(d, i) : d.low[i];
        if (hay.find(lab) != std::string::npos) {
          out.push_back(i);
          seen.insert(i);
        }
      }
  return out;
}

// Field value: text after the label colon plus continuation lines.
std::string field_text(const Doc& d, size_t i, bool lowercase = true) {
  const auto& L = lowercase ? d.low : d.raw;
  std::string s = L[i];
  if (is_label_line(d.raw[i])) s = s.substr(s.find(':') + 1);
  for (size_t j = i + 1; j < L.size() && j < i + 15; ++j) {
    if (d.raw[j].empty() || is_label_line(d.raw[j])) break;
    s += " " + L[j];
  }
  return trim(s);
}

double range_value(const std::smatch& m, bool lower_bound = false) {
  double a = std::stod(m[1].str());
  if (m[2].matched && !lower_bound) return (a + std::stod(m[2].str())) / 2;
  return a;
}

std::optional<double> parse_time_h(const std::string& s) {
  static const std::regex re(NUM + R"(\s*(hours?|hrs?|h|minutes?|mins?|days?)\b)");
  std::smatch m;
  if (!std::regex_search(s, m, re)) return {};
  double v = range_value(m);
  std::string u = m[3].str();
  if (u.rfind("min", 0) == 0) v /= 60;
  if (u.rfind("day", 0) == 0) v *= 24;
  return v;
}

std::optional<double> parse_volume_L(const std::string& s) {
  static const std::regex re(NUM + R"(\s*(l/kg|ml/kg|liters?|litres?|ml|l)\b)");
  std::smatch m;
  if (!std::regex_search(s, m, re)) return {};
  double v = range_value(m);
  std::string u = m[3].str();
  if (u == "l/kg") v *= kBodyWeightKg;
  else if (u == "ml/kg") v *= kBodyWeightKg / 1000;
  else if (u == "ml") v /= 1000;
  return v;
}

std::optional<double> parse_clearance_Lh(const std::string& s) {
  static const std::regex re(NUM + R"(\s*(ml/min/kg|l/h/kg|l/hr/kg|ml/min|l/min|l/hours?|l/hr|l/h)\b)");
  std::smatch m;
  if (!std::regex_search(s, m, re)) return {};
  double v = range_value(m);
  std::string u = m[3].str();
  if (u == "ml/min/kg") v *= 0.06 * kBodyWeightKg;
  else if (u == "l/h/kg" || u == "l/hr/kg") v *= kBodyWeightKg;
  else if (u == "ml/min") v *= 0.06;
  else if (u == "l/min") v *= 60;
  return v;
}

std::optional<double> parse_conc_mgL(const std::string& s, bool lower_bound = false) {
  static const std::regex re(NUM + R"(\s*(ng/ml|mcg/ml|µg/ml|μg/ml|ug/ml|mg/ml|mg/l|mcg/l|ug/l|ng/l)\b)");
  std::smatch m;
  if (!std::regex_search(s, m, re)) return {};
  double v = range_value(m, lower_bound);
  std::string u = m[3].str();
  if (u == "ng/ml" || u == "mcg/l" || u == "ug/l") v *= 1e-3;
  else if (u == "mg/ml") v *= 1e3;
  else if (u == "ng/l") v *= 1e-6;
  return v;
}

std::optional<double> parse_auc_mghL(const std::string& s) {
  static const std::regex re(NUM +
                             R"(\s*(ng|mcg|µg|μg|ug|mg)\s*(?:\*|·|x|\.|-|•)?\s*(?:hours?|hrs?|h)\s*/\s*(ml|l)\b)");
  std::smatch m;
  if (!std::regex_search(s, m, re)) return {};
  double v = range_value(m);
  std::string mass = m[3].str(), vol = m[4].str();
  double mg = mass == "ng" ? 1e-6 : mass == "mg" ? 1.0 : 1e-3;
  double L = vol == "ml" ? 1e-3 : 1.0;
  return v * mg / L;
}

std::optional<double> parse_percent(const std::string& s) {
  static const std::regex re(NUM + R"(\s*%)");
  std::smatch m;
  if (!std::regex_search(s, m, re)) return {};
  return range_value(m) / 100.0;
}

std::optional<double> parse_fraction_or_percent(const std::string& s) {
  if (auto p = parse_percent(s)) return p;
  static const std::regex re(R"((?:^|[^0-9.])(0?\.[0-9]+|1\.0)(?![0-9]))");
  std::smatch m;
  if (std::regex_search(s, m, re)) return std::stod(m[1].str());
  return {};
}

std::optional<double> parse_dose_mg(const std::string& s) {
  static const std::regex re(R"(([0-9]+(?:\.[0-9]+)?)\s*(mg|mcg|g)\b)");
  std::smatch m;
  if (!std::regex_search(s, m, re)) return {};
  double v = std::stod(m[1].str());
  if (m[2].str() == "g") v *= 1000;
  if (m[2].str() == "mcg") v /= 1000;
  return v;
}

std::optional<double> parse_number(const std::string& s) {
  static const std::regex re(R"(([0-9]+(?:\.[0-9]+)?))");
  std::smatch m;
  if (!std::regex_search(s, m, re)) return {};
  return std::stod(m[1].str());
}

std::optional<double> parse_interval_h(const std::string& s) {
  static const std::regex every(R"(every\s+([0-9]+(?:\.[0-9]+)?)\s*(hours?|hrs?|h|days?)\b)");
  std::smatch m;
  if (std::regex_search(s, m, every)) {
    double v = std::stod(m[1].str());
    return m[2].str().rfind("day", 0) == 0 ? v * 24 : v;
  }
  auto has = [&](std::initializer_list<const char*> keys) {
    for (auto k : keys)
      if (s.find(k) != std::string::npos) return true;
    return false;
  };
  if (has({"twice daily", "twice a day", "b.i.d", " bid"})) return 12.0;
  if (has({"three times daily", "three times a day", "t.i.d", " tid"})) return 8.0;
  if (has({"four times daily", "four times a day", "q.i.d", " qid"})) return 6.0;
  if (has({"once weekly", "once a week", "weekly"})) return 168.0;
  if (has({"once daily", "once a day", "daily", "q.d"})) return 24.0;
  return {};
}

template <class Parse>
std::optional<double> find_field(const Doc& d, const std::vector<std::string>& labels, Parse parse,
                                 size_t* line = nullptr, const std::vector<std::string>& exclude = {}) {
  for (size_t i : candidates(d, labels, exclude)) {
    if (auto v = parse(field_text(d, i))) {
      if (line) *line = i;
      return v;
    }
  }
  return {};
}

bool word_boundary(const std::string& s, size_t pos, size_t len) {
  bool before = pos == 0 || !std::isalnum((unsigned char)s[pos - 1]);
  bool after = pos + len >= s.size() || !std::isalnum((unsigned char)s[pos + len]);
  return before && after;
}

struct Span {
  size_t pos, len;
  std::string key;
};

// Keep longest matches; drop a match if a strictly longer accepted match overlaps it.
std::vector<Span> resolve_overlaps(std::vector<Span> hits) {
  std::stable_sort(hits.begin(), hits.end(), [](const Span& a, const Span& b) { return a.len > b.len; });
  std::vector<Span> acc;
  for (auto& h : hits) {
    bool drop = false;
    for (auto& a : acc)
      if (h.pos < a.pos + a.len && a.pos < h.pos + h.len && a.len > h.len) drop = true;
    if (!drop) acc.push_back(h);
  }
  std::sort(acc.begin(), acc.end(), [](const Span& a, const Span& b) { return a.pos < b.pos; });
  return acc;
}

void find_all(const std::string& hay, const std::string& needle, const std::string& key, std::vector<Span>& out) {
  if (needle.size() < 3) return;
  size_t p = 0;
  while ((p = hay.find(needle, p)) != std::string::npos) {
    if (word_boundary(hay, p, needle.size())) out.push_back({p, needle.size(), key});
    p += needle.size();
  }
}

Direction direction_near(const std::string& low, const Span& hit, const std::vector<Span>& all_hits) {
  // sentence containing the hit
  size_t b = 0;
  for (size_t i = hit.pos; i > 0; --i)
    if ((low[i - 1] == '.' || low[i - 1] == ';') && (i >= low.size() || low[i] == ' ')) {
      b = i;
      break;
    }
  size_t e = low.size();
  for (size_t i = hit.pos + hit.len; i < low.size(); ++i)
    if ((low[i] == '.' || low[i] == ';') && (i + 1 >= low.size() || low[i + 1] == ' ')) {
      e = i;
      break;
    }
  struct Kw {
    const char* stem;
    Direction dir;
  };
  static const Kw kws[] = {{"inhibit", Direction::Down},  {"antagon", Direction::Down},
                           {"block", Direction::Down},    {"suppress", Direction::Down},
                           {"agonist", Direction::Up},    {"activat", Direction::Up},
                           {"stimulat", Direction::Up},   {"enhanc", Direction::Up},
                           {"potentiat", Direction::Up},  {"opener", Direction::Up}};
  double best = 1e18;
  Direction dir = Direction::Any;
  for (auto& kw : kws) {
    std::string stem = kw.stem;
    size_t p = b;
    while ((p = low.find(stem, p)) != std::string::npos && p < e) {
      bool inside_target_name = false;
      for (auto& h : all_hits)
        if (p >= h.pos && p < h.pos + h.len) inside_target_name = true;
      bool is_antagonist = stem == "agonist" && p >= 3 && low.compare(p - 3, 3, "ant") == 0;
      if (!inside_target_name && !is_antagonist) {
        double dist = std::abs((double)p - (double)hit.pos);
        if (dist < best) {
          best = dist;
          dir = kw.dir;
        }
      }
      p += stem.size();
    }
  }
  return dir;
}

}  // namespace

std::vector<TargetAction> mine_targets(const std::string& text, const KnowledgeBase& kb) {
  std::string low = lower(text);
  std::vector<Span> hits;
  for (auto& t : kb.targets) {
    for (auto& syn : t.synonyms) find_all(low, lower(syn), t.symbol, hits);
    // gene symbol written verbatim (case-sensitive); short symbols are too ambiguous
    if (t.symbol.size() >= 4 && t.symbol.find('_') == std::string::npos) {
      size_t p = 0;
      while ((p = text.find(t.symbol, p)) != std::string::npos) {
        if (word_boundary(text, p, t.symbol.size())) hits.push_back({p, t.symbol.size(), t.symbol});
        p += t.symbol.size();
      }
    }
  }
  auto acc = resolve_overlaps(hits);
  std::vector<TargetAction> out;
  std::set<std::string> seen;
  for (auto& h : acc) {
    if (seen.count(h.key)) continue;
    seen.insert(h.key);
    TargetAction ta;
    ta.symbol = h.key;
    ta.dir = direction_near(low, h, acc);
    ta.action = "mined: \"" + text.substr(h.pos, h.len) + "\"";
    out.push_back(ta);
  }
  return out;
}

std::vector<std::string> mine_indications(const std::string& text, const KnowledgeBase& kb) {
  std::string low = lower(text);
  std::vector<Span> hits;
  for (auto& d : kb.diseases) {
    std::vector<std::string> names = d.synonyms;
    std::string k = d.key;
    std::replace(k.begin(), k.end(), '_', ' ');
    names.push_back(k);
    for (auto& n : names) find_all(low, lower(n), d.key, hits);
  }
  std::vector<std::string> out;
  for (auto& h : resolve_overlaps(hits))
    if (std::find(out.begin(), out.end(), h.key) == out.end()) out.push_back(h.key);
  return out;
}

DrugRecord parse_monograph(const std::string& text, const KnowledgeBase& kb, std::vector<std::string>& warnings) {
  Doc d = make_doc(text);
  DrugRecord r;

  // ---- name
  for (size_t i : candidates(d, {"drug name", "generic name", "proprietary name", "name"}, {}, true)) {
    std::string v = field_text(d, i, false);
    if (!v.empty()) {
      r.name = v.substr(0, 60);
      break;
    }
  }
  if (r.name.empty()) {
    for (auto& l : d.raw)
      if (!l.empty()) {
        r.name = l.substr(0, 40);
        break;
      }
    warnings.push_back("no 'Drug Name:' field; using first line as name");
  }

  // ---- numeric fields
  if (auto v = find_field(d, {"molecular weight", "mol. wt", "molecular mass"}, parse_number)) r.mw = *v;

  size_t dose_line = 0;
  if (auto v = find_field(d, {"recommended dose", "usual dose", "starting dose", "dose", "dosage"}, parse_dose_mg,
                          &dose_line)) {
    r.dose_mg = *v;
    std::string ft = field_text(d, dose_line);
    if (auto tau = parse_interval_h(ft)) r.tau_h = *tau;
    if (ft.find("intravenous") != std::string::npos || ft.find(" iv ") != std::string::npos) r.route = "iv";
  }
  if (std::isnan(r.tau_h))
    if (auto v = find_field(d, {"dosing interval", "frequency", "dosage and administration"}, parse_interval_h))
      r.tau_h = *v;

  if (auto v = find_field(d, {"absolute bioavailability", "oral bioavailability", "bioavailability"},
                          parse_fraction_or_percent))
    r.F = std::clamp(*v, 0.0, 1.0);

  if (auto v = find_field(d, {"half-life", "half life", "t1/2", "t½"}, parse_time_h)) r.half_life_h = *v;

  size_t v_line = 0, cl_line = 0;
  if (auto v = find_field(d, {"volume of distribution", "vss", "vd/f", "vd", "apparent volume"}, parse_volume_L,
                          &v_line))
    r.V = *v;
  if (auto v = find_field(d,
                          {"total body clearance", "systemic clearance", "total clearance", "plasma clearance",
                           "apparent clearance", "apparent oral clearance", "cl/f", "clearance"},
                          parse_clearance_Lh, &cl_line, {"renal clearance", "creatinine"}))
    r.CL = *v;
  if (auto v = find_field(d, {"absorption rate constant", "absorption rate"}, parse_number)) r.ka = *v;
  if (auto v = find_field(d, {"tmax", "t max", "time to peak", "time to maximum"}, parse_time_h)) r.obs_tmax = *v;
  if (auto v = find_field(d,
                          {"cmax", "c max", "peak plasma concentration", "maximum plasma concentration",
                           "peak serum concentration", "peak concentration"},
                          [](const std::string& s) { return parse_conc_mgL(s); }))
    r.obs_cmax = *v;
  if (auto v = find_field(d, {"auc", "area under"}, parse_auc_mghL)) r.obs_auc = *v;
  if (auto v = find_field(d, {"protein binding", "protein bound", "bound to plasma proteins"}, parse_percent))
    r.fu = std::clamp(1.0 - *v, 0.0, 1.0);
  if (auto v = find_field(d,
                          {"therapeutic range", "therapeutic concentration", "effective concentration",
                           "target concentration"},
                          [](const std::string& s) { return parse_conc_mgL(s, true); }))
    r.ceff = *v;

  // ---- mechanism -> targets
  auto mech = candidates(d, {"mechanism of action", "mechanism", "pharmacodynamics"}, {}, true);
  std::string mech_text;
  if (!mech.empty()) {
    mech_text = field_text(d, mech[0], false);
  } else {
    mech_text = text;
    warnings.push_back("no 'Mechanism of Action' section; mining targets from the whole document (less reliable)");
  }
  auto explicit_t = candidates(d, {"targets", "molecular target"}, {}, true);
  if (!explicit_t.empty()) {
    static const std::regex re(R"(([A-Za-z0-9_]{2,})\s*\(([^)]+)\))");
    std::string ft = field_text(d, explicit_t[0], false);
    for (std::sregex_iterator it(ft.begin(), ft.end(), re), end; it != end; ++it) {
      TargetAction ta;
      ta.symbol = (*it)[1].str();
      for (auto& c : ta.symbol) c = (char)std::toupper((unsigned char)c);
      ta.action = (*it)[2].str();
      ta.dir = parse_direction(ta.action);
      r.targets.push_back(ta);
    }
  }
  for (auto& t : mine_targets(mech_text, kb)) {
    bool dup = false;
    for (auto& e : r.targets) dup |= e.symbol == t.symbol;
    if (!dup) r.targets.push_back(t);
  }
  if (r.targets.empty()) warnings.push_back("no known molecular target recognised; repurposing will rely on nothing");

  // ---- indications
  auto ind = candidates(d, {"indications and usage", "indications", "indicated for", "approved uses"}, {}, true);
  std::string ind_text;
  if (!ind.empty()) {
    ind_text = field_text(d, ind[0], false);
  } else {
    for (size_t i = 0; i < d.low.size(); ++i)
      if (d.low[i].find("indicated") != std::string::npos) ind_text += " " + d.raw[i];
    if (ind_text.empty()) warnings.push_back("no indications section found");
  }
  r.indications = mine_indications(ind_text, kb);

  // ---- derivations
  bool cl_app = r.CL > 0 && (d.low[cl_line].find("cl/f") != std::string::npos ||
                             d.low[cl_line].find("apparent") != std::string::npos);
  bool v_app = r.V > 0 && (d.low[v_line].find("/f") != std::string::npos ||
                           d.low[v_line].find("apparent") != std::string::npos);
  if ((cl_app || v_app) && std::isnan(r.F)) {
    r.F = 1.0;
    r.notes.push_back("apparent (oral) CL/F or V/F reported and F unknown -> F set to 1 (exposure-equivalent)");
    if (cl_app != v_app) warnings.push_back("only one of CL and V is apparent (/F); exposure may be biased");
  } else if (cl_app || v_app) {
    if (cl_app) r.CL *= r.F;
    if (v_app) r.V *= r.F;
    r.notes.push_back("apparent CL/F or V/F converted to systemic values using F");
  }
  if (std::isnan(r.V) && r.CL > 0 && r.half_life_h > 0) {
    r.V = r.CL * r.half_life_h / kLn2;
    r.notes.push_back("V derived from CL x t1/2 / ln2 = " + fmt(r.V) + " L");
  }
  if (std::isnan(r.CL) && r.V > 0 && r.half_life_h > 0) {
    r.CL = r.V * kLn2 / r.half_life_h;
    r.notes.push_back("CL derived from V x ln2 / t1/2 = " + fmt(r.CL) + " L/h");
  }
  if (r.CL > 0 && r.V > 0 && r.half_life_h > 0) {
    double implied = kLn2 * r.V / r.CL;
    double ratio = implied / r.half_life_h;
    if (ratio > 2 || ratio < 0.5)
      warnings.push_back("stated t1/2 " + fmt(r.half_life_h) + " h vs CL/V-implied " + fmt(implied) +
                         " h (multi-compartment or terminal-phase t1/2?)");
  }
  if (std::isnan(r.F)) {
    r.F = lower(r.route) == "iv" ? 1.0 : 1.0;
    if (lower(r.route) != "iv") warnings.push_back("bioavailability not found; assumed F = 1 (exposure may be overestimated)");
  }
  if (std::isnan(r.ka) && lower(r.route) != "iv") {
    if (r.obs_tmax > 0 && r.CL > 0 && r.V > 0) {
      r.ka = ka_from_tmax(r.obs_tmax, r.CL / r.V);
      if (!std::isnan(r.ka)) r.notes.push_back("ka solved from Tmax " + fmt(r.obs_tmax) + " h -> " + fmt(r.ka) + " /h");
    }
    if (std::isnan(r.ka)) {
      r.ka = 1.0;
      warnings.push_back("absorption rate unknown; assumed ka = 1 /h");
    }
  }
  if (std::isnan(r.tau_h)) {
    r.tau_h = 24;
    warnings.push_back("dosing interval not found; assumed every 24 h");
  }
  if (std::isnan(r.dose_mg)) warnings.push_back("dose not found; PK simulation unavailable");
  if (!(r.CL > 0) || !(r.V > 0)) warnings.push_back("CL and V not both available; PK simulation unavailable");
  return r;
}

DrugRecord ingest_document(const std::string& path, const KnowledgeBase& kb, std::vector<std::string>& warnings,
                           std::string* method) {
  PdfExtractResult ex = extract_pdf_text(path);
  if (!ex.ok) throw std::runtime_error(path + ": " + ex.error);
  if (method) *method = ex.method;
  DrugRecord r = parse_monograph(ex.text, kb, warnings);
  r.source = path;
  return r;
}

}  // namespace pkd
