#include "pkd/pk.hpp"

#include <algorithm>
#include <cmath>

namespace pkd {

double PKModel::lambda_z() const {
  double k10 = kel();
  double beta = k10;
  if (Q > 0 && V2 > 0) {
    double k12 = Q / V, k21 = Q / V2;
    double s = k10 + k12 + k21;
    beta = (s - std::sqrt(s * s - 4 * k10 * k21)) / 2;
  }
  if (!iv && ka > 0) return std::min(beta, ka);  // flip-flop kinetics if ka < beta
  return beta;
}

double PKModel::terminal_half_life() const { return kLn2 / lambda_z(); }

void integrate(const PKModel& m, double dose, double tau, int n_doses, double t_end, const Observer& obs) {
  const bool two = m.Q > 0 && m.V2 > 0;
  const double kel = m.kel();
  const double k12 = two ? m.Q / m.V : 0, k21 = two ? m.Q / m.V2 : 0;
  const double ka = m.iv ? 0 : m.ka;
  const double fastest = std::max({kel + k12, k21, ka, 1e-6});
  const double hmax = std::min(0.1, 0.1 / fastest);

  double g = 0, c = 0, p = 0;  // amounts (mg) in gut, central, peripheral
  auto deriv = [&](double G, double C, double P, double& dG, double& dC, double& dP) {
    dG = -ka * G;
    dC = ka * G - (kel + k12) * C + k21 * P;
    dP = k12 * C - k21 * P;
  };
  int given = 0;
  auto give = [&] {
    if (m.iv)
      c += dose;
    else
      g += m.F * dose;
    ++given;
  };

  double t = 0;
  give();
  obs(0, c / m.V);
  while (t < t_end - 1e-12) {
    double next_dose = given < n_doses ? given * tau : INFINITY;
    double seg_end = std::min(t_end, next_dose);
    double span = seg_end - t;
    int n = std::max(1, (int)std::ceil(span / hmax));
    double h = span / n;
    for (int i = 0; i < n; ++i) {
      double g1, c1, p1, g2, c2, p2, g3, c3, p3, g4, c4, p4;
      deriv(g, c, p, g1, c1, p1);
      deriv(g + h / 2 * g1, c + h / 2 * c1, p + h / 2 * p1, g2, c2, p2);
      deriv(g + h / 2 * g2, c + h / 2 * c2, p + h / 2 * p2, g3, c3, p3);
      deriv(g + h * g3, c + h * c3, p + h * p3, g4, c4, p4);
      g += h / 6 * (g1 + 2 * g2 + 2 * g3 + g4);
      c += h / 6 * (c1 + 2 * c2 + 2 * c3 + c4);
      p += h / 6 * (p1 + 2 * p2 + 2 * p3 + p4);
      t = (i == n - 1) ? seg_end : t + h;
      obs(t, c / m.V);
    }
    if (given < n_doses && std::abs(t - next_dose) < 1e-9) {
      give();
      if (m.iv) obs(t, c / m.V);
    }
  }
}

std::vector<Sample> simulate(const PKModel& m, double dose, double tau, int n_doses, double t_end, double dt_out) {
  std::vector<Sample> out;
  double next = 0;
  integrate(m, dose, tau, n_doses, t_end, [&](double t, double c) {
    if (t + 1e-9 >= next) {
      out.push_back({t, c});
      next = t + dt_out;
    }
  });
  return out;
}

SingleDoseMetrics single_dose_metrics(const PKModel& m, double dose) {
  double lz = m.lambda_z();
  double t_end = std::min(std::max(24.0, 8 * kLn2 / lz), 24.0 * 90);
  double cmax = 0, tmax = 0, auc = 0, pt = 0, pc = 0, clast = 0;
  bool first = true;
  integrate(m, dose, 1e9, 1, t_end, [&](double t, double c) {
    if (c > cmax) {
      cmax = c;
      tmax = t;
    }
    if (!first && t > pt) auc += (t - pt) * (c + pc) / 2;
    first = false;
    pt = t;
    pc = c;
    clast = c;
  });
  auc += clast / lz;
  return {cmax, tmax, auc, kLn2 / lz};
}

SteadyStateMetrics steady_state(const PKModel& m, double dose, double tau) {
  double lz = m.lambda_z();
  double th = kLn2 / lz;
  int n = (int)std::ceil(7 * th / tau) + 1;
  n = std::clamp(n, 2, 500);
  double t0 = (n - 1) * tau, t_end = n * tau;
  double cmax = 0, cmin = INFINITY, auc = 0, pt = -1, pc = 0;
  integrate(m, dose, tau, n, t_end, [&](double t, double c) {
    if (t + 1e-9 < t0) return;
    cmax = std::max(cmax, c);
    cmin = std::min(cmin, c);
    if (pt >= 0 && t > pt) auc += (t - pt) * (c + pc) / 2;
    pt = t;
    pc = c;
  });
  SingleDoseMetrics sd = single_dose_metrics(m, dose);
  SteadyStateMetrics ss;
  ss.cmax = cmax;
  ss.cmin = cmin;
  ss.cavg = auc / tau;
  ss.accumulation = sd.cmax > 0 ? cmax / sd.cmax : nan_v();
  ss.doses_to_90pct = std::max(1, (int)std::ceil(std::log(10.0) / (lz * tau)));
  return ss;
}

double bateman(const PKModel& m, double dose, double t) {
  double k = m.kel();
  if (m.iv) return dose / m.V * std::exp(-k * t);
  if (std::abs(m.ka - k) < 1e-12) return m.F * dose * k * t * std::exp(-k * t) / m.V;
  return m.F * dose * m.ka / (m.V * (m.ka - k)) * (std::exp(-k * t) - std::exp(-m.ka * t));
}

double ka_from_tmax(double tmax, double k) {
  if (!(tmax > 0) || !(k > 0)) return nan_v();
  auto tm = [&](double ka) {
    if (std::abs(ka - k) < 1e-9 * k) return 1.0 / k;
    return std::log(ka / k) / (ka - k);
  };
  double lo = k * 1e-4, hi = k * 1e4;  // tm() is decreasing in ka
  if (tmax > tm(lo) || tmax < tm(hi)) return nan_v();
  for (int i = 0; i < 200; ++i) {
    double mid = std::sqrt(lo * hi);
    if (tm(mid) > tmax)
      lo = mid;
    else
      hi = mid;
  }
  return std::sqrt(lo * hi);
}

bool model_from(const DrugRecord& d, PKModel& m, std::string* why) {
  if (!(d.CL > 0) || !(d.V > 0)) {
    if (why) *why = "clearance and/or volume of distribution unavailable";
    return false;
  }
  m = PKModel{};
  m.CL = d.CL;
  m.V = d.V;
  m.Q = d.Q;
  m.V2 = d.V2;
  m.iv = lower(d.route) == "iv";
  m.F = std::isnan(d.F) ? 1.0 : d.F;
  m.ka = std::isnan(d.ka) ? 1.0 : d.ka;
  return true;
}

}  // namespace pkd
