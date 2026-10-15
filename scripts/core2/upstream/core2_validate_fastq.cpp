// Exhaustive technical validation of four full SRA-exported FASTQ streams.
// Checks every spot ID, sequence/quality length, and FASTQ record framing.
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
int main(int argc, char **argv) {
  try {
    if (argc != 4) throw std::runtime_error("usage: validate_fastq directory run expected_records");
    const std::string dir = argv[1], accession = argv[2];
    const unsigned long long expected = std::stoull(argv[3]);
    const std::vector<size_t> lengths{8, 24, 50, 50};
    std::cout << "{\"streams\":[";
    for (size_t i = 0; i < 4; ++i) {
      const std::string name = accession + "_" + std::to_string(i+1) + ".fastq";
      std::ifstream f(dir + "/" + name);
      if (!f) throw std::runtime_error("cannot open " + name);
      std::string header, sequence, plus, quality;
      unsigned long long count = 0;
      while (std::getline(f, header)) {
        ++count;
        if (!std::getline(f, sequence) || !std::getline(f, plus) || !std::getline(f, quality))
          throw std::runtime_error("truncated record in " + name);
        if (header != "@" + accession + "." + std::to_string(count))
          throw std::runtime_error("unexpected spot ID in " + name + " record " + std::to_string(count));
        if (plus != "+" || sequence.size() != lengths[i] || quality.size() != lengths[i])
          throw std::runtime_error("invalid record structure/length in " + name);
        if (sequence.find_first_not_of("ACGTN") != std::string::npos)
          throw std::runtime_error("invalid base in " + name);
        for (unsigned char q : quality) if (q < 33 || q > 126)
          throw std::runtime_error("invalid FASTQ quality in " + name);
      }
      if (f.bad() || count != expected) throw std::runtime_error("read error or unexpected record count in " + name);
      if (i) std::cout << ",";
      std::cout << "{\"file\":\"" << name << "\",\"records\":" << count
                << ",\"sequence_length\":" << lengths[i] << "}";
    }
    std::cout << "],\"every_spot_id_checked\":true,\"synchronized\":true,\"sequence_quality_lengths_checked\":true}\n";
  } catch (const std::exception &e) { std::cerr << e.what() << "\n"; return 1; }
}
