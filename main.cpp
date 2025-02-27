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
    void readGenotypes(std::vector<std::vector<int>>& genotypes) {
        genotypes.resize(sample_count, std::vector<int>(variant_count));

        // This is a simplified version - actual implementation depends on storage mode
        for (uint32_t variant = 0; variant < variant_count; ++variant) {
            for (uint32_t sample = 0; sample < sample_count; ++sample) {
                // Read 2-bit genotype (0,1,2, or 3 for missing)
                uint8_t byte;
                pgen_file.read(reinterpret_cast<char*>(&byte), 1);
                int genotype = byte & 0x03; // Extract first genotype from byte
                genotypes[sample][variant] = (genotype == 3) ? -1 : genotype; // -1 for missing

                // Note: This assumes simplest format (mode 1).
                // Real implementation would need to handle different compression modes
            }
        }
    }

    void readVariantInfo(std::vector<std::string>& variant_ids) {
        std::string line;
        // Skip header line in .pvar
        std::getline(pvar_file, line);

        while (std::getline(pvar_file, line)) {
            std::string id = line.substr(line.find('\t') + 1);
            id = id.substr(0, id.find('\t'));
            variant_ids.push_back(id);
        }
    }

    void readSampleInfo(std::vector<std::string>& sample_ids) {
        std::string line;
        // Skip header line in .psam
        std::getline(psam_file, line);

        while (std::getline(psam_file, line)) {
            std::string id = line.substr(0, line.find('\t'));
            sample_ids.push_back(id);
        }
    }
};

// Example usage
int main() {
    try {
        Plink2Reader reader("data2.pgen", "data2.pvar", "data2.psam");

        // Read sample IDs
        std::vector<std::string> sample_ids;
        reader.readSampleInfo(sample_ids);

        // Read variant IDs
        std::vector<std::string> variant_ids;
        reader.readVariantInfo(variant_ids);

        // Read genotype data
        std::vector<std::vector<int>> genotypes;
        reader.readGenotypes(genotypes);

        // Print some example data
        std::cout << "Samples: " << sample_ids.size() << "\n";
        std::cout << "Variants: " << variant_ids.size() << "\n";
        std::cout << "Genotypes for sample 0:\n";
        std::cout << "Genotypes size: " << genotypes[0].size() << std::endl;

        cout << genotypes.size() << endl;


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