#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
using namespace std;

class Plink2Reader {
private:
    std::ifstream pgen_file;
    std::ifstream pvar_file;
    std::ifstream psam_file;

    uint32_t variant_count;
    uint32_t sample_count;
    uint64_t file_size;

public:
    Plink2Reader(const std::string& pgen_path,
        const std::string& pvar_path,
        const std::string& psam_path) {

        // Open files
        pgen_file.open(pgen_path, std::ios::binary);
        pvar_file.open(pvar_path);
        psam_file.open(psam_path);

        if (!pgen_file.is_open() || !pvar_file.is_open() || !psam_file.is_open()) {
            throw std::runtime_error("Failed to open one or more PLINK2 files");
        }

        // Read header from pgen file
        readHeader();
    }

    ~Plink2Reader() {
        if (pgen_file.is_open()) pgen_file.close();
        if (pvar_file.is_open()) pvar_file.close();
        if (psam_file.is_open()) psam_file.close();
    }

private:
    void readHeader() {
        // Read magic numbers (first 2 bytes should be 0x6c, 0x1b)
        char magic[2];
        pgen_file.read(magic, 2);
        if (magic[0] != 0x6c || magic[1] != 0x1b) {
            throw std::runtime_error("Invalid PGEN file format");
        }

        // Read mode byte
        char mode;
        pgen_file.read(&mode, 1);

        // Read variant and sample counts
        pgen_file.read(reinterpret_cast<char*>(&variant_count), 4);
        pgen_file.read(reinterpret_cast<char*>(&sample_count), 4);

        // Get file size
        pgen_file.seekg(0, std::ios::end);
        file_size = pgen_file.tellg();
        pgen_file.seekg(11); // Return to start of genotype data
    }

public:
    void readGenotypesChunk(std::vector<std::vector<int>>& genotypes, uint32_t start_variant, uint32_t end_variant, uint32_t start_sample, uint32_t end_sample) {
        if (end_variant > variant_count || end_sample > sample_count) {
            throw std::out_of_range("Requested chunk is out of range");
        }

        uint32_t num_variants = end_variant - start_variant;
        uint32_t num_samples = end_sample - start_sample;

        genotypes.resize(num_samples, std::vector<int>(num_variants));

        // Calculate the starting position in the file
        uint64_t start_pos = 11 + (start_variant * sample_count + start_sample) / 4;
        pgen_file.seekg(start_pos);

        for (uint32_t variant = start_variant; variant < end_variant; ++variant) {
            for (uint32_t sample = start_sample; sample < end_sample; ++sample) {
                // Read 2-bit genotype (0,1,2, or 3 for missing)
                uint8_t byte;
                pgen_file.read(reinterpret_cast<char*>(&byte), 1);
                int genotype = byte & 0x03; // Extract first genotype from byte
                genotypes[sample - start_sample][variant - start_variant] = (genotype == 3) ? -1 : genotype; // -1 for missing
            }
        }
    }

    void readVariantInfoChunk(std::vector<std::string>& variant_ids, uint32_t start_variant, uint32_t end_variant) {
        if (end_variant > variant_count) {
            throw std::out_of_range("Requested chunk is out of range");
        }

        std::string line;
        // Skip header line in .pvar
        std::getline(pvar_file, line);

        // Skip to the start variant
        for (uint32_t i = 0; i < start_variant; ++i) {
            std::getline(pvar_file, line);
        }

        // Read the chunk of variants
        for (uint32_t i = start_variant; i < end_variant; ++i) {
            std::getline(pvar_file, line);
            std::string id = line.substr(line.find('\t') + 1);
            id = id.substr(0, id.find('\t'));
            variant_ids.push_back(id);
        }
    }

    void readSampleInfoChunk(std::vector<std::string>& sample_ids, uint32_t start_sample, uint32_t end_sample) {
        if (end_sample > sample_count) {
            throw std::out_of_range("Requested chunk is out of range");
        }

        std::string line;
        // Skip header line in .psam
        std::getline(psam_file, line);

        // Skip to the start sample
        for (uint32_t i = 0; i < start_sample; ++i) {
            std::getline(psam_file, line);
        }

        // Read the chunk of samples
        for (uint32_t i = start_sample; i < end_sample; ++i) {
            std::getline(psam_file, line);
            std::string id = line.substr(0, line.find('\t'));
            sample_ids.push_back(id);
        }
    }
};

// Example usage
int main() {
    try {
        Plink2Reader reader("plink2.pgen", "plink2.pvar", "plink2.psam");

        // Read a chunk of sample IDs
        std::vector<std::string> sample_ids;
        reader.readSampleInfoChunk(sample_ids, 0, 10); // Read first 10 samples

        // Read a chunk of variant IDs
        std::vector<std::string> variant_ids;
        reader.readVariantInfoChunk(variant_ids, 0, 100); // Read first 100 variants

        // Read a chunk of genotype data
        std::vector<std::vector<int>> genotypes;
        reader.readGenotypesChunk(genotypes, 0, 100, 0, 10); // Read first 100 variants for first 10 samples

        // Print some example data
        std::cout << "Samples: " << sample_ids.size() << "\n";
        std::cout << "Variants: " << variant_ids.size() << "\n";
        std::cout << "Genotypes for sample 0:\n";
        for (int i = 0; i < (int)genotypes[0].size(); ++i) {
            std::cout << genotypes[0][i] << " ";
        }
        std::cout << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}