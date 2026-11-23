// Minimal, dependency-light PDF text extraction.
// Handles the common case of text-based PDFs with FlateDecode content streams
// and simple (single-byte) fonts. Complex PDFs (CID/Type0 fonts, encrypted)
// are delegated to `pdftotext` when it is installed.
#include "pkd/pdf_text.hpp"

#include <zlib.h>

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "pkd/util.hpp"

namespace pkd {

bool inflate_zlib(const std::string& in, std::string& out) {
  out.clear();
  z_stream zs{};
  if (inflateInit(&zs) != Z_OK) return false;
  zs.next_in = (Bytef*)in.data();
  zs.avail_in = (uInt)in.size();
  std::vector<char> buf(1 << 16);
  for (;;) {
    zs.next_out = (Bytef*)buf.data();
    zs.avail_out = (uInt)buf.size();
    int ret = inflate(&zs, Z_NO_FLUSH);
    size_t have = buf.size() - zs.avail_out;
    out.append(buf.data(), have);
    if (ret == Z_STREAM_END) break;
    if (ret != Z_OK) {
      inflateEnd(&zs);
      return !out.empty();  // keep partial output from slightly corrupt streams
    }
    if (zs.avail_in == 0 && have == 0) break;
  }
  inflateEnd(&zs);
  return true;
}

namespace {

// WinAnsi / Latin-1 byte -> UTF-8
void append_pdf_byte(std::string& out, unsigned char c) {
  if (c == '\n' || c == '\t') {
    out += (char)c;
    return;
  }
  if (c < 0x20) return;
  if (c < 0x80) {
    out += (char)c;
    return;
  }
  switch (c) {
    case 0x91: case 0x92: out += '\''; return;
    case 0x93: case 0x94: out += '"'; return;
    case 0x95: out += '*'; return;
    case 0x96: case 0x97: out += '-'; return;
    case 0xA0: out += ' '; return;
    default: break;
  }
  if (c >= 0xA0) {  // Latin-1 code point -> 2-byte UTF-8
    out += (char)(0xC0 | (c >> 6));
    out += (char)(0x80 | (c & 0x3F));
  }
}

enum class Kind { Num, Str, Arr, Name, Other };
struct Operand {
  Kind kind;
  double num = 0;
  std::string str;  // decoded bytes (Str), concatenated text (Arr), name
};

bool is_delim(char c) {
  return std::isspace((unsigned char)c) || c == '(' || c == ')' || c == '<' || c == '>' || c == '[' ||
         c == ']' || c == '{' || c == '}' || c == '/' || c == '%';
}

std::string read_literal(const std::string& s, size_t& i) {  // s[i] == '('
  std::string out;
  int depth = 1;
  ++i;
  while (i < s.size() && depth > 0) {
    char c = s[i];
    if (c == '\\' && i + 1 < s.size()) {
      char n = s[++i];
      switch (n) {
        case 'n': out += '\n'; break;
        case 'r': out += '\r'; break;
        case 't': out += '\t'; break;
        case 'b': out += '\b'; break;
        case 'f': out += '\f'; break;
        case '(': out += '('; break;
        case ')': out += ')'; break;
        case '\\': out += '\\'; break;
        case '\r':
          if (i + 1 < s.size() && s[i + 1] == '\n') ++i;
          break;
        case '\n': break;
        default:
          if (n >= '0' && n <= '7') {
            int v = n - '0';
            for (int k = 0; k < 2 && i + 1 < s.size() && s[i + 1] >= '0' && s[i + 1] <= '7'; ++k)
              v = v * 8 + (s[++i] - '0');
            out += (char)(v & 0xFF);
          } else {
            out += n;
          }
      }
      ++i;
      continue;
    }
    if (c == '(') depth++;
    if (c == ')') {
      depth--;
      if (depth == 0) {
        ++i;
        break;
      }
    }
    out += c;
    ++i;
  }
  return out;
}

std::string read_hex(const std::string& s, size_t& i) {  // s[i] == '<'
  ++i;
  std::string digits;
  while (i < s.size() && s[i] != '>') {
    if (std::isxdigit((unsigned char)s[i])) digits += s[i];
    ++i;
  }
  if (i < s.size()) ++i;
  if (digits.size() % 2) digits += '0';
  std::string out;
  for (size_t k = 0; k + 1 < digits.size(); k += 2) out += (char)std::strtol(digits.substr(k, 2).c_str(), nullptr, 16);
  return out;
}

void ensure_newline(std::string& out) {
  if (!out.empty() && out.back() != '\n') out += '\n';
}
void append_text(std::string& out, const std::string& bytes) {
  for (unsigned char c : bytes) append_pdf_byte(out, c);
}

}  // namespace

std::string text_from_content_stream(const std::string& s) {
  std::string out;
  std::vector<Operand> st;
  size_t i = 0;
  auto last_num = [&](size_t back) -> double {  // back=0 -> last operand
    if (st.size() > back && st[st.size() - 1 - back].kind == Kind::Num) return st[st.size() - 1 - back].num;
    return 0;
  };
  while (i < s.size()) {
    char c = s[i];
    if (std::isspace((unsigned char)c)) {
      ++i;
      continue;
    }
    if (c == '%') {
      while (i < s.size() && s[i] != '\n' && s[i] != '\r') ++i;
      continue;
    }
    if (c == '(') {
      st.push_back({Kind::Str, 0, read_literal(s, i)});
      continue;
    }
    if (c == '<') {
      if (i + 1 < s.size() && s[i + 1] == '<') {
        i += 2;
        continue;
      }
      st.push_back({Kind::Str, 0, read_hex(s, i)});
      continue;
    }
    if (c == '>') {
      ++i;
      continue;
    }
    if (c == '[') {
      ++i;
      Operand arr{Kind::Arr, 0, ""};
      while (i < s.size() && s[i] != ']') {
        char a = s[i];
        if (std::isspace((unsigned char)a)) {
          ++i;
        } else if (a == '(') {
          arr.str += read_literal(s, i);
        } else if (a == '<') {
          arr.str += read_hex(s, i);
        } else if (std::isdigit((unsigned char)a) || a == '-' || a == '+' || a == '.') {
          size_t b = i;
          while (i < s.size() && !is_delim(s[i])) ++i;
          double kern = std::atof(s.substr(b, i - b).c_str());
          if (kern < -200) arr.str += ' ';  // large negative kern = word gap
        } else {
          ++i;
        }
      }
      if (i < s.size()) ++i;
      st.push_back(arr);
      continue;
    }
    if (c == ']' || c == '{' || c == '}') {
      ++i;
      continue;
    }
    if (c == '/') {
      size_t b = ++i;
      while (i < s.size() && !is_delim(s[i])) ++i;
      st.push_back({Kind::Name, 0, s.substr(b, i - b)});
      continue;
    }
    if (std::isdigit((unsigned char)c) || c == '-' || c == '+' || c == '.') {
      size_t b = i;
      while (i < s.size() && !is_delim(s[i])) ++i;
      st.push_back({Kind::Num, std::atof(s.substr(b, i - b).c_str()), ""});
      continue;
    }
    // operator
    size_t b = i;
    while (i < s.size() && !is_delim(s[i])) ++i;
    if (i == b) {
      ++i;
      continue;
    }
    std::string op = s.substr(b, i - b);
    if (op == "Tj") {
      if (!st.empty() && st.back().kind == Kind::Str) append_text(out, st.back().str);
    } else if (op == "TJ") {
      if (!st.empty() && st.back().kind == Kind::Arr) append_text(out, st.back().str);
    } else if (op == "'" || op == "\"") {
      ensure_newline(out);
      if (!st.empty() && st.back().kind == Kind::Str) append_text(out, st.back().str);
    } else if (op == "T*" || op == "Tm" || op == "ET") {
      ensure_newline(out);
    } else if (op == "Td" || op == "TD") {
      double ty = last_num(0), tx = last_num(1);
      if (ty != 0)
        ensure_newline(out);
      else if (tx > 0 && !out.empty() && out.back() != ' ' && out.back() != '\n')
        out += ' ';
    } else if (op == "BI") {  // inline image: skip to EI
      size_t e = s.find("EI", i);
      i = (e == std::string::npos) ? s.size() : e + 2;
    }
    st.clear();
  }
  return out;
}

std::string extract_text_builtin(const std::string& data) {
  std::string out;
  size_t pos = 0;
  while ((pos = data.find("stream", pos)) != std::string::npos) {
    if (pos >= 3 && data.compare(pos - 3, 3, "end") == 0) {
      pos += 6;
      continue;
    }
    size_t start = pos + 6;
    if (start >= data.size() || (data[start] != '\r' && data[start] != '\n')) {
      pos = start;
      continue;
    }
    if (data[start] == '\r') ++start;
    if (start < data.size() && data[start] == '\n') ++start;
    size_t end = data.find("endstream", start);
    if (end == std::string::npos) break;
    size_t obj = data.rfind("obj", pos);
    std::string dict = (obj == std::string::npos || pos - obj > 4096) ? "" : data.substr(obj, pos - obj);
    std::string raw = data.substr(start, end - start);
    while (!raw.empty() && (raw.back() == '\n' || raw.back() == '\r')) raw.pop_back();
    pos = end + 9;

    static const char* skip_markers[] = {"/Image", "/FontFile", "/Length1", "/XRef", "/ObjStm", "/Metadata",
                                         "/DCTDecode", "/JPXDecode", "/CCITTFaxDecode", "/JBIG2Decode",
                                         "/ASCII85Decode", "/LZWDecode", "/RunLengthDecode"};
    bool skip = false;
    for (auto m : skip_markers)
      if (dict.find(m) != std::string::npos) skip = true;
    if (skip) continue;

    std::string content;
    if (dict.find("/FlateDecode") != std::string::npos) {
      if (!inflate_zlib(raw, content)) continue;
    } else {
      content = raw;
    }
    if (content.find("BT") == std::string::npos) continue;
    std::string t = text_from_content_stream(content);
    if (!t.empty()) {
      ensure_newline(out);
      out += t;
    }
  }
  return out;
}

double text_quality(const std::string& text) {
  if (text.empty()) return 0;
  size_t good = 0, letters = 0;
  for (unsigned char c : text) {
    if (std::isalnum(c) || std::isspace(c) || std::ispunct(c) || c >= 0x80) good++;
    if (std::isalpha(c)) letters++;
  }
  double q = (double)good / text.size();
  if (letters < 50) q *= (double)letters / 50.0;
  return q;
}

namespace {
std::string shell_quote(const std::string& s) {
  std::string out = "'";
  for (char c : s) {
    if (c == '\'')
      out += "'\\''";
    else
      out += c;
  }
  return out + "'";
}

bool run_pdftotext(const std::string& path, std::string& out) {
  std::string cmd = "pdftotext -enc UTF-8 " + shell_quote(path) + " - 2>/dev/null";
  FILE* p = popen(cmd.c_str(), "r");
  if (!p) return false;
  char buf[8192];
  size_t n;
  out.clear();
  while ((n = fread(buf, 1, sizeof(buf), p)) > 0) out.append(buf, n);
  int rc = pclose(p);
  return rc == 0 && !trim(out).empty();
}
}  // namespace

PdfExtractResult extract_pdf_text(const std::string& path, bool allow_external) {
  PdfExtractResult r;
  std::string data;
  if (!read_file(path, data)) {
    r.error = "cannot read file: " + path;
    return r;
  }
  if (data.compare(0, 5, "%PDF-") != 0) {  // treat as plain text monograph
    r.ok = true;
    r.text = data;
    r.method = "plain-text";
    return r;
  }
  if (data.find("/Encrypt") != std::string::npos && !allow_external) {
    r.error = "encrypted PDF (needs pdftotext)";
    return r;
  }
  std::string builtin = extract_text_builtin(data);
  double qb = text_quality(builtin);
  if (qb >= 0.9 && data.find("/Encrypt") == std::string::npos) {
    r.ok = true;
    r.text = builtin;
    r.method = "builtin";
    return r;
  }
  if (allow_external) {
    std::string ext;
    if (run_pdftotext(path, ext) && text_quality(ext) >= qb) {
      r.ok = true;
      r.text = ext;
      r.method = "pdftotext";
      return r;
    }
  }
  if (!trim(builtin).empty()) {
    r.ok = true;
    r.text = builtin;
    r.method = "builtin(low-quality)";
    return r;
  }
  r.error = "no extractable text (scanned image PDF? try OCR first, or install poppler's pdftotext)";
  return r;
}

}  // namespace pkd
