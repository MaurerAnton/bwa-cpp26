#include <bwa/io/seq_io.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/vector.hpp>
#include <cstring>

namespace bwa::io {

// Fast quality score conversion
uint8_t phred33_to_score(char c) noexcept {
    return static_cast<uint8_t>(std::max(0, static_cast<int>(c) - 33));
}

char score_to_phred33(uint8_t score) noexcept {
    return static_cast<char>(std::min(126, static_cast<int>(score) + 33));
}

// Parse FASTQ quality line to scores
void parse_quality(std::string_view qual, core::Vector<uint8_t>& scores) {
    scores.resize(qual.size());
    for (size_t i = 0; i < qual.size(); ++i) {
        scores[i] = phred33_to_score(qual[i]);
    }
}

// Encode scores to FASTQ quality
void encode_quality(std::span<const uint8_t> scores, core::PmrString& qual) {
    qual.clear();
    qual.reserve(scores.size());
    for (uint8_t s : scores) {
        qual.kputc(score_to_phred33(s));
    }
}

} // namespace bwa::io