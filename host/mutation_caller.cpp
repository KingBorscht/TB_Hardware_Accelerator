// CPU-side mutation caller for the TB Hardware Accelerator (SCOPE v2).
//
// This program is intentionally independent of PYNQ/DMA libraries. A board
// driver converts FPGA traceback output into alignments.tsv, allowing the same
// caller to be tested on a workstation with recorded FPGA fixtures.
//
// Build: c++ -std=c++17 -O2 -Wall -Wextra -pedantic mutation_caller.cpp -o tb_mutation_caller
// Run:   ./tb_mutation_caller references.tsv catalogue.tsv alignments.tsv [calls.tsv]

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace tb {

enum class Allele { Ref, Alt, Other, Ambiguous, NotCovered };

struct ReferenceWindow {
    int amplicon_id;
    std::string name;
    std::string sequence;       // gene-oriented H37Rv sequence
    int genome_start_1based;
    int genome_step;            // +1 or -1 in gene orientation
    std::string reference_version;
};

struct CatalogueVariant {
    std::string id;
    int amplicon_id;
    int reference_offset;       // zero-based; left anchor for insertions
    std::string reference;      // '-' represents an insertion
    std::string alternate;      // '-' represents a deletion
    std::string gene;
    std::string annotation;
    std::string catalogue_version;
};

struct Alignment {
    std::string read_id;
    int amplicon_id;
    bool aligned;
    int score;
    int query_begin;
    int query_end;
    int reference_begin;
    int reference_end;
    std::string cigar;          // complete M/I/D traceback path
    std::string read_sequence;  // normalized, gene-oriented read
};

struct CigarOp { char code; int length; };

struct Event {
    int offset;                 // zero-based reference start or insertion anchor
    std::string reference;
    std::string alternate;
};

struct ReadObservations {
    std::vector<Event> events;
    std::set<int> covered_offsets;
    std::set<int> ambiguous_offsets;
};

struct Counts { int ref = 0, alt = 0, other = 0, ambiguous = 0, not_covered = 0; };

[[noreturn]] void fail(const std::string& message) { throw std::runtime_error(message); }

std::vector<std::string> split_tab(const std::string& line) {
    std::vector<std::string> fields;
    std::stringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) fields.push_back(field);
    return fields;
}

int as_int(const std::string& text, const std::string& context) {
    try {
        size_t consumed = 0;
        const int value = std::stoi(text, &consumed);
        if (consumed != text.size()) fail("invalid integer in " + context + ": " + text);
        return value;
    } catch (const std::exception&) {
        fail("invalid integer in " + context + ": " + text);
    }
}

std::string uppercase_dna(std::string value, const std::string& context) {
    for (char& base : value) {
        base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
        if (base != 'A' && base != 'C' && base != 'G' && base != 'T' && base != 'N' && base != '-') {
            fail("invalid DNA base in " + context);
        }
    }
    return value;
}

std::vector<CigarOp> parse_cigar(const std::string& cigar) {
    if (cigar.empty()) fail("empty CIGAR");
    std::vector<CigarOp> operations;
    int length = 0;
    bool saw_digit = false;
    for (char character : cigar) {
        if (std::isdigit(static_cast<unsigned char>(character))) {
            saw_digit = true;
            length = length * 10 + (character - '0');
            continue;
        }
        if (!saw_digit || length == 0 || (character != 'M' && character != 'I' && character != 'D')) {
            fail("invalid CIGAR: " + cigar);
        }
        operations.push_back({character, length});
        length = 0;
        saw_digit = false;
    }
    if (saw_digit) fail("CIGAR ends with a length: " + cigar);
    return operations;
}

std::vector<ReferenceWindow> read_references(const std::string& path) {
    std::ifstream input(path);
    if (!input) fail("cannot open reference file: " + path);
    std::vector<ReferenceWindow> references;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tab(line);
        if (fields.size() != 6) fail("references.tsv line " + std::to_string(line_number) + " needs 6 fields");
        ReferenceWindow item {as_int(fields[0], "references.tsv"), fields[1], uppercase_dna(fields[2], "reference"),
                              as_int(fields[3], "references.tsv"), as_int(fields[4], "references.tsv"), fields[5]};
        if (item.genome_step != 1 && item.genome_step != -1) fail("reference genome_step must be +1 or -1");
        references.push_back(std::move(item));
    }
    return references;
}

std::vector<CatalogueVariant> read_catalogue(const std::string& path) {
    std::ifstream input(path);
    if (!input) fail("cannot open catalogue file: " + path);
    std::vector<CatalogueVariant> variants;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tab(line);
        if (fields.size() != 8) fail("catalogue.tsv line " + std::to_string(line_number) + " needs 8 fields");
        CatalogueVariant item {fields[0], as_int(fields[1], "catalogue.tsv"), as_int(fields[2], "catalogue.tsv"),
                               uppercase_dna(fields[3], "catalogue reference"), uppercase_dna(fields[4], "catalogue alternate"),
                               fields[5], fields[6], fields[7]};
        if (item.reference == "-" && item.alternate == "-") fail("catalogue variant cannot be '-' to '-'");
        variants.push_back(std::move(item));
    }
    return variants;
}

std::vector<Alignment> read_alignments(const std::string& path) {
    std::ifstream input(path);
    if (!input) fail("cannot open alignment file: " + path);
    std::vector<Alignment> alignments;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') continue;
        const auto fields = split_tab(line);
        if (fields.size() != 10) fail("alignments.tsv line " + std::to_string(line_number) + " needs 10 fields");
        const std::string aligned_text = fields[2];
        if (aligned_text != "0" && aligned_text != "1") fail("alignment aligned field must be 0 or 1");
        Alignment item {fields[0], as_int(fields[1], "alignments.tsv"), aligned_text == "1", as_int(fields[3], "alignments.tsv"),
                        as_int(fields[4], "alignments.tsv"), as_int(fields[5], "alignments.tsv"),
                        as_int(fields[6], "alignments.tsv"), as_int(fields[7], "alignments.tsv"),
                        fields[8], uppercase_dna(fields[9], "read sequence")};
        alignments.push_back(std::move(item));
    }
    return alignments;
}

ReadObservations walk_traceback(const Alignment& alignment, const ReferenceWindow& reference) {
    if (!alignment.aligned) return {};
    int query = alignment.query_begin;
    int ref = alignment.reference_begin;
    if (query < 0 || ref < 0) fail("negative traceback endpoint for read " + alignment.read_id);
    ReadObservations result;
    for (const CigarOp& op : parse_cigar(alignment.cigar)) {
        if (op.code == 'M') {
            for (int i = 0; i < op.length; ++i, ++query, ++ref) {
                if (query >= static_cast<int>(alignment.read_sequence.size()) || ref >= static_cast<int>(reference.sequence.size())) {
                    fail("M operation exceeds sequence bounds for read " + alignment.read_id);
                }
                result.covered_offsets.insert(ref);
                const char read_base = alignment.read_sequence[query];
                const char reference_base = reference.sequence[ref];
                if (read_base == 'N') {
                    result.ambiguous_offsets.insert(ref);
                } else if (read_base != reference_base) {
                    result.events.push_back({ref, std::string(1, reference_base), std::string(1, read_base)});
                }
            }
        } else if (op.code == 'D') {
            if (ref + op.length > static_cast<int>(reference.sequence.size())) fail("D operation exceeds reference bounds");
            result.events.push_back({ref, reference.sequence.substr(ref, op.length), "-"});
            for (int i = 0; i < op.length; ++i) result.covered_offsets.insert(ref + i);
            ref += op.length;
        } else { // I; the event is anchored to the preceding reference base, or 0 at the window start.
            if (query + op.length > static_cast<int>(alignment.read_sequence.size())) fail("I operation exceeds read bounds");
            result.events.push_back({std::max(0, ref - 1), "-", alignment.read_sequence.substr(query, op.length)});
            query += op.length;
        }
    }
    if (query != alignment.query_end || ref != alignment.reference_end) {
        fail("CIGAR does not consume stated endpoints for read " + alignment.read_id);
    }
    return result;
}

bool overlaps(const Event& event, const CatalogueVariant& variant) {
    const int event_last = event.offset + std::max(1, static_cast<int>(event.reference == "-" ? 0 : event.reference.size())) - 1;
    const int variant_last = variant.reference_offset + std::max(1, static_cast<int>(variant.reference == "-" ? 0 : variant.reference.size())) - 1;
    return event.offset <= variant_last && variant.reference_offset <= event_last;
}

Allele classify(const CatalogueVariant& variant, const ReadObservations& observed) {
    for (const Event& event : observed.events) {
        if (event.offset == variant.reference_offset && event.reference == variant.reference && event.alternate == variant.alternate) {
            return Allele::Alt;
        }
        if (overlaps(event, variant)) return Allele::Other;
    }
    if (observed.ambiguous_offsets.count(variant.reference_offset)) return Allele::Ambiguous;
    if (observed.covered_offsets.count(variant.reference_offset)) return Allele::Ref;
    return Allele::NotCovered;
}

const char* text(Allele allele) {
    switch (allele) {
        case Allele::Ref: return "REF";
        case Allele::Alt: return "ALT";
        case Allele::Other: return "OTHER_ALLELE";
        case Allele::Ambiguous: return "AMBIGUOUS";
        case Allele::NotCovered: return "NOT_COVERED";
    }
    return "UNKNOWN";
}

void add_count(Counts& counts, Allele allele) {
    switch (allele) {
        case Allele::Ref: ++counts.ref; break;
        case Allele::Alt: ++counts.alt; break;
        case Allele::Other: ++counts.other; break;
        case Allele::Ambiguous: ++counts.ambiguous; break;
        case Allele::NotCovered: ++counts.not_covered; break;
    }
}

std::string status(const Counts& counts) {
    const int eligible = counts.ref + counts.alt;
    if (eligible < 50) return "INSUFFICIENT_COVERAGE";
    const double fraction = static_cast<double>(counts.alt) / eligible;
    if (fraction >= 0.70) return "RESISTANT_SUPPORTED";
    if (fraction >= 0.10) return "MIXED";
    return "NOT_SUPPORTED";
}

void run(const std::string& references_path, const std::string& catalogue_path, const std::string& alignments_path, std::ostream& output) {
    std::map<int, ReferenceWindow> references;
    for (auto& item : read_references(references_path)) {
        if (!references.emplace(item.amplicon_id, std::move(item)).second) fail("duplicate amplicon in references.tsv");
    }
    const auto catalogue = read_catalogue(catalogue_path);
    std::map<std::pair<std::string, std::string>, Allele> evidence; // (catalogue ID, read ID)

    for (const Alignment& alignment : read_alignments(alignments_path)) {
        const auto reference_it = references.find(alignment.amplicon_id);
        if (reference_it == references.end()) fail("alignment references unknown amplicon");
        const ReadObservations observed = walk_traceback(alignment, reference_it->second);
        for (const CatalogueVariant& variant : catalogue) {
            if (variant.amplicon_id != alignment.amplicon_id) continue;
            const auto key = std::make_pair(variant.id, alignment.read_id);
            if (!evidence.emplace(key, classify(variant, observed)).second) fail("duplicate read/catalogue evidence: " + variant.id + "/" + alignment.read_id);
        }
    }

    output << "catalogue_id\tamplicon_id\tgene\tannotation\th37rv_position\treference\talternate\teligible_coverage\talt_count\tallele_fraction\tstatus\tref_count\tother_count\tambiguous_count\tnot_covered_count\treference_version\tcatalogue_version\n";
    for (const CatalogueVariant& variant : catalogue) {
        Counts counts;
        for (const auto& entry : evidence) if (entry.first.first == variant.id) add_count(counts, entry.second);
        const int eligible = counts.ref + counts.alt;
        const double fraction = eligible == 0 ? 0.0 : static_cast<double>(counts.alt) / eligible;
        const ReferenceWindow& reference = references.at(variant.amplicon_id);
        const int genomic_position = reference.genome_start_1based + reference.genome_step * variant.reference_offset;
        output << variant.id << '\t' << variant.amplicon_id << '\t' << variant.gene << '\t' << variant.annotation << '\t'
               << genomic_position << '\t' << variant.reference << '\t' << variant.alternate << '\t' << eligible << '\t'
               << counts.alt << '\t' << fraction << '\t' << status(counts) << '\t' << counts.ref << '\t'
               << counts.other << '\t' << counts.ambiguous << '\t' << counts.not_covered << '\t'
               << reference.reference_version << '\t' << variant.catalogue_version << '\n';
    }
}

} // namespace tb

int main(int argc, char** argv) {
    if (argc != 4 && argc != 5) {
        std::cerr << "Usage: " << argv[0] << " references.tsv catalogue.tsv alignments.tsv [calls.tsv]\n";
        return 2;
    }
    try {
        if (argc == 5) {
            std::ofstream output(argv[4]);
            if (!output) tb::fail("cannot open output file");
            tb::run(argv[1], argv[2], argv[3], output);
        } else {
            tb::run(argv[1], argv[2], argv[3], std::cout);
        }
    } catch (const std::exception& error) {
        std::cerr << "mutation_caller: " << error.what() << '\n';
        return 1;
    }
}
