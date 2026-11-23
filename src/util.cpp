#include "pkd/util.hpp"

#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace pkd {

std::string trim(const std::string& s) {
  size_t b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  size_t e = s.find_last_not_of(" \t\r\n");
  return s.substr(b, e - b + 1);
}

std::string lower(std::string s) {
  for (auto& c : s) c = (char)std::tolower((unsigned char)c);
  return s;
}

std::vector<std::string> split(const std::string& s, char delim) {
  std::vector<std::string> out;
  std::string cur;
  for (char c : s) {
    if (c == delim) {
      out.push_back(cur);
      cur.clear();
    } else {
      cur += c;
    }
  }
  out.push_back(cur);
  return out;
}

std::vector<std::string> split_nonempty(const std::string& s, char delim) {
  std::vector<std::string> out;
  for (auto& p : split(s, delim)) {
    std::string t = trim(p);
    if (!t.empty()) out.push_back(t);
  }
  return out;
}

std::string join(const std::vector<std::string>& v, const std::string& sep) {
  std::string out;
  for (size_t i = 0; i < v.size(); ++i) {
    if (i) out += sep;
    out += v[i];
  }
  return out;
}

std::string normalize_key(const std::string& s) {
  std::string out;
  bool last_us = false;
  for (unsigned char c : s) {
    if (std::isalnum(c)) {
      out += (char)std::tolower(c);
      last_us = false;
    } else if (!out.empty() && !last_us) {
      out += '_';
      last_us = true;
    }
  }
  while (!out.empty() && out.back() == '_') out.pop_back();
  return out;
}

bool read_file(const std::string& path, std::string& out) {
  std::ifstream f(path, std::ios::binary);
  if (!f) return false;
  std::ostringstream ss;
  ss << f.rdbuf();
  out = ss.str();
  return true;
}

bool ends_with(const std::string& s, const std::string& suffix) {
  return s.size() >= suffix.size() && s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

double to_double(const std::string& s, double def) {
  std::string t = trim(s);
  if (t.empty()) return def;
  try {
    return std::stod(t);
  } catch (...) {
    return def;
  }
}

std::string fmt(double v, int sig) {
  if (std::isnan(v)) return "n/a";
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.*g", sig, v);
  return buf;
}

int CsvTable::col(const std::string& name) const {
  for (size_t i = 0; i < header.size(); ++i)
    if (header[i] == name) return (int)i;
  throw std::runtime_error("CSV column not found: " + name);
}

std::string CsvTable::get(size_t row, const std::string& name) const {
  int c = col(name);
  if (row >= rows.size() || c >= (int)rows[row].size()) return "";
  return rows[row][c];
}

CsvTable read_csv(const std::string& path) {
  std::string content;
  if (!read_file(path, content)) throw std::runtime_error("cannot read " + path);
  CsvTable t;
  std::istringstream in(content);
  std::string line;
  while (std::getline(in, line)) {
    std::string tl = trim(line);
    if (tl.empty() || tl[0] == '#') continue;
    std::vector<std::string> cells;
    for (auto& c : split(tl, ',')) cells.push_back(trim(c));
    if (t.header.empty())
      t.header = cells;
    else
      t.rows.push_back(cells);
  }
  return t;
}

}  // namespace pkd
