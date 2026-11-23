#pragma once
#include <functional>
#include <vector>

#include "pkd/knowledge_base.hpp"

namespace pkd {

// Linear compartmental model:
//   gut --ka--> central (V) <--Q--> peripheral (V2)
//                  |CL
// Oral doses enter the gut scaled by F; IV doses enter the central compartment.
struct PKModel {
  double F = 1, ka = 1, CL = 1, V = 10, Q = 0, V2 = 0;
  bool iv = false;
  double kel() const { return CL / V; }
  double lambda_z() const;  // terminal rate constant, 1/h
  double terminal_half_life() const;
};

struct Sample {
  double t, c;  // h, mg/L
};

struct SingleDoseMetrics {
  double cmax, tmax, auc_inf, half_life;
};

struct SteadyStateMetrics {
  double cmax, cmin, cavg, accumulation;
  int doses_to_90pct;
};

using Observer = std::function<void(double t, double c)>;

// RK4 integration of the ODE system with dose events every tau hours.
void integrate(const PKModel& m, double dose_mg, double tau_h, int n_doses, double t_end, const Observer& obs);
std::vector<Sample> simulate(const PKModel& m, double dose_mg, double tau_h, int n_doses, double t_end,
                             double dt_out);

SingleDoseMetrics single_dose_metrics(const PKModel& m, double dose_mg);
SteadyStateMetrics steady_state(const PKModel& m, double dose_mg, double tau_h);

// Closed-form one-compartment solution (Bateman function), used for verification.
double bateman(const PKModel& m, double dose_mg, double t);

// Solve tmax = ln(ka/k)/(ka-k) for ka. Returns NaN if no solution.
double ka_from_tmax(double tmax_h, double k);

// Build a model from a drug record; returns false (with reason) if CL or V is missing.
bool model_from(const DrugRecord& d, PKModel& m, std::string* why = nullptr);

}  // namespace pkd
