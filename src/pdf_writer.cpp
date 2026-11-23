// Tiny PDF 1.4 writer: Helvetica text pages with FlateDecode content streams.
#include <zlib.h>

#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "pkd/pdf_text.hpp"

namespace pkd {
namespace {

// UTF-8 -> WinAnsi bytes (subset)
std::string to_winansi(const std::string& s) {
  std::string out;
  for (size_t i = 0; i < s.size();) {
    unsigned char c = s[i];
    unsigned cp = 0;
    int len = 1;
    if (c < 0x80) {
      cp = c;
    } else if ((c >> 5) == 6 && i + 1 < s.size()) {
      cp = ((c & 0x1F) << 6) | (s[i + 1] & 0x3F);
      len = 2;
    } else if ((c >> 4) == 14 && i + 2 < s.size()) {
      cp = ((c & 0x0F) << 12) | ((s[i + 1] & 0x3F) << 6) | (s[i + 2] & 0x3F);
      len = 3;
    } else {
      cp = '?';
      len = (c >> 3) == 30 ? 4 : 1;
    }
    i += len;
    if (cp < 0x80)
      out += (char)cp;
    else if (cp >= 0xA0 && cp <= 0xFF)
      out += (char)cp;
    else if (cp == 0x2013) out += (char)0x96;
    else if (cp == 0x2014) out += (char)0x97;
    else if (cp == 0x2018 || cp == 0x2019) out += (char)0x92;
    else if (cp == 0x201C || cp == 0x201D) out += (char)0x94;
    else if (cp == 0x2022) out += (char)0x95;
    else if (cp == 0x03BC) out += (char)0xB5;  // Greek mu -> micro sign
    else out += '?';
  }
  return out;
}

std::string escape_pdf(const std::string& s) {
  std::string out;
  for (char c : s) {
    if (c == '(' || c == ')' || c == '\\') out += '\\';
    out += c;
  }
  return out;
}

std::vector<std::string> wrap(const std::vector<std::string>& lines, size_t width) {
  std::vector<std::string> out;
  for (const auto& l : lines) {
    if (l.size() <= width) {
      out.push_back(l);
      continue;
    }
    std::string rest = l;
    while (rest.size() > width) {
      size_t cut = rest.rfind(' ', width);
      if (cut == std::string::npos || cut == 0) cut = width;
      out.push_back(rest.substr(0, cut));
      rest = rest.substr(cut + (rest[cut] == ' ' ? 1 : 0));
    }
    out.push_back(rest);
  }
  return out;
}

}  // namespace

std::string build_text_pdf(const std::vector<std::string>& in_lines, bool compress) {
  const size_t kLinesPerPage = 54;
  std::vector<std::string> lines = wrap(in_lines, 92);
  if (lines.empty()) lines.push_back("");
  size_t npages = (lines.size() + kLinesPerPage - 1) / kLinesPerPage;

  std::vector<std::string> objs;  // objs[k] is object k+1 body
  objs.push_back("<< /Type /Catalog /Pages 2 0 R >>");
  std::string kids;
  for (size_t p = 0; p < npages; ++p) kids += std::to_string(4 + 2 * p) + " 0 R ";
  objs.push_back("<< /Type /Pages /Kids [" + kids + "] /Count " + std::to_string(npages) + " >>");
  objs.push_back("<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>");

  for (size_t p = 0; p < npages; ++p) {
    std::string content = "BT\n/F1 10 Tf\n13 TL\n54 760 Td\n";
    for (size_t i = p * kLinesPerPage; i < std::min(lines.size(), (p + 1) * kLinesPerPage); ++i)
      content += "(" + escape_pdf(to_winansi(lines[i])) + ") Tj T*\n";
    content += "ET\n";
    std::string body, dict;
    if (compress) {
      uLongf cap = compressBound(content.size());
      std::string z(cap, '\0');
      compress2((Bytef*)&z[0], &cap, (const Bytef*)content.data(), content.size(), 9);
      z.resize(cap);
      body = z;
      dict = "<< /Length " + std::to_string(body.size()) + " /Filter /FlateDecode >>";
    } else {
      body = content;
      dict = "<< /Length " + std::to_string(body.size()) + " >>";
    }
    objs.push_back("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Resources << /Font << /F1 3 0 R >> >> "
                   "/Contents " + std::to_string(5 + 2 * p) + " 0 R >>");
    objs.push_back(dict + "\nstream\n" + body + "\nendstream");
  }

  std::string pdf = "%PDF-1.4\n%\xE2\xE3\xCF\xD3\n";
  std::vector<size_t> offsets;
  for (size_t k = 0; k < objs.size(); ++k) {
    offsets.push_back(pdf.size());
    pdf += std::to_string(k + 1) + " 0 obj\n" + objs[k] + "\nendobj\n";
  }
  size_t xref = pdf.size();
  pdf += "xref\n0 " + std::to_string(objs.size() + 1) + "\n0000000000 65535 f \n";
  char buf[32];
  for (size_t off : offsets) {
    std::snprintf(buf, sizeof(buf), "%010zu 00000 n \n", off);
    pdf += buf;
  }
  pdf += "trailer\n<< /Size " + std::to_string(objs.size() + 1) + " /Root 1 0 R >>\nstartxref\n" +
         std::to_string(xref) + "\n%%EOF\n";
  return pdf;
}

bool write_text_pdf(const std::string& path, const std::vector<std::string>& lines, bool compress) {
  std::ofstream f(path, std::ios::binary);
  if (!f) return false;
  std::string pdf = build_text_pdf(lines, compress);
  f.write(pdf.data(), (std::streamsize)pdf.size());
  return (bool)f;
}

}  // namespace pkd
