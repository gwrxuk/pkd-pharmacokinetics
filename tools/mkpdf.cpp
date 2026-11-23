// mkpdf: convert plain-text monographs into simple text PDFs.
//   mkpdf in.txt out.pdf
//   mkpdf <dir-with-txt> <out-dir>      (converts every *.txt)
#include <cstdio>
#include <filesystem>
#include <string>

#include "pkd/pdf_text.hpp"
#include "pkd/util.hpp"

namespace fs = std::filesystem;

static bool convert(const fs::path& in, const fs::path& out) {
  std::string text;
  if (!pkd::read_file(in.string(), text)) {
    std::fprintf(stderr, "cannot read %s\n", in.string().c_str());
    return false;
  }
  std::vector<std::string> lines;
  for (auto& l : pkd::split(text, '\n')) {
    if (!l.empty() && l.back() == '\r') l.pop_back();
    lines.push_back(l);
  }
  if (!pkd::write_text_pdf(out.string(), lines, true)) {
    std::fprintf(stderr, "cannot write %s\n", out.string().c_str());
    return false;
  }
  std::printf("%s -> %s\n", in.string().c_str(), out.string().c_str());
  return true;
}

int main(int argc, char** argv) {
  if (argc != 3) {
    std::fprintf(stderr, "usage: mkpdf in.txt out.pdf | mkpdf in_dir out_dir\n");
    return 2;
  }
  fs::path in = argv[1], out = argv[2];
  if (fs::is_directory(in)) {
    fs::create_directories(out);
    bool ok = true;
    for (auto& e : fs::directory_iterator(in))
      if (e.path().extension() == ".txt") ok &= convert(e.path(), out / e.path().stem().concat(".pdf"));
    return ok ? 0 : 1;
  }
  return convert(in, out) ? 0 : 1;
}
