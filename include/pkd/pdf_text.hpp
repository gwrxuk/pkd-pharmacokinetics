#pragma once
#include <string>
#include <vector>

namespace pkd {

struct PdfExtractResult {
  bool ok = false;
  std::string text;
  std::string method;  // "builtin", "pdftotext", "plain-text"
  std::string error;
};

// Extract text from a PDF. Tries the built-in parser (FlateDecode streams,
// Tj/TJ/'/" operators); if that yields little readable text (e.g. CID fonts),
// falls back to the external `pdftotext` (poppler) when available.
PdfExtractResult extract_pdf_text(const std::string& path, bool allow_external = true);

// Built-in pieces, exposed for testing.
std::string extract_text_builtin(const std::string& pdf_bytes);
std::string text_from_content_stream(const std::string& content);
bool inflate_zlib(const std::string& in, std::string& out);
double text_quality(const std::string& text);  // 0..1, fraction of "normal" chars

// Minimal PDF writer (Helvetica, WinAnsi) used to generate sample monographs.
std::string build_text_pdf(const std::vector<std::string>& lines, bool compress = true);
bool write_text_pdf(const std::string& path, const std::vector<std::string>& lines, bool compress = true);

}  // namespace pkd
