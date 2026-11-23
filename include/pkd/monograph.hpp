#pragma once
#include <string>
#include <vector>

#include "pkd/knowledge_base.hpp"

namespace pkd {

// Parse a drug monograph / label (plain text) into a DrugRecord:
//  * PK parameters with unit conversion (mL/min, L/kg, ng/mL, minutes, ranges ...)
//  * targets + direction mined from the "Mechanism of Action" section using the
//    target synonym dictionary in data/targets.csv
//  * label indications mapped to disease keys in data/diseases.csv
// Missing parameters are derived where pharmacokinetically valid
// (V = CL*t1/2/ln2, ka from Tmax, ...) and every derivation is recorded in notes.
DrugRecord parse_monograph(const std::string& text, const KnowledgeBase& kb, std::vector<std::string>& warnings);

// Extract text from a PDF (or read a .txt) and parse it.
// Throws std::runtime_error if no text could be extracted.
DrugRecord ingest_document(const std::string& path, const KnowledgeBase& kb, std::vector<std::string>& warnings,
                           std::string* method = nullptr);

// Exposed for testing.
std::vector<TargetAction> mine_targets(const std::string& mechanism_text, const KnowledgeBase& kb);
std::vector<std::string> mine_indications(const std::string& text, const KnowledgeBase& kb);

}  // namespace pkd
