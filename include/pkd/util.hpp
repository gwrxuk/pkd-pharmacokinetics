#pragma once
#include <limits>
#include <string>
#include <vector>

namespace pkd {

inline double nan_v() { return std::numeric_limits<double>::quiet_NaN(); }
constexpr double kLn2 = 0.69314718055994530942;
constexpr double kBodyWeightKg = 70.0;  // reference adult for per-kg values

std::string trim(const std::string& s);
std::string lower(std::string s);
std::vector<std::string> split(const std::string& s, char delim);
std::vector<std::string> split_nonempty(const std::string& s, char delim);  // trimmed
std::string join(const std::vector<std::string>& v, const std::string& sep);
std::string normalize_key(const std::string& s);  // "Type 2 Diabetes" -> "type_2_diabetes"
bool read_file(const std::string& path, std::string& out);
bool ends_with(const std::string& s, const std::string& suffix);
double to_double(const std::string& s, double def);
std::string fmt(double v, int sig = 3);  // "n/a" for NaN

struct CsvTable {
  std::vector<std::string> header;
  std::vector<std::vector<std::string>> rows;
  int col(const std::string& name) const;  // throws if missing
  std::string get(size_t row, const std::string& name) const;
};
// Comma-separated, no quoting; '#' lines and blank lines are skipped.
CsvTable read_csv(const std::string& path);

}  // namespace pkd
