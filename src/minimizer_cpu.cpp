#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cstdint>
#include <cstring>

#define HASH_WIDTH 64
#define DATA_SIZE_BYTES (256 * 1024 * 1024)  
#define NUM_WORDS (DATA_SIZE_BYTES / sizeof(uint64_t))

inline uint64_t bfc_hash_64(uint64_t key, uint64_t mask) {
    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
}


inline void parallel_hash_calc(uint64_t input_data, uint64_t &output_data) {
    uint64_t mask = ~0ULL;
    output_data = bfc_hash_64(input_data, mask);
}

void minimizer_cpu(const uint64_t *input_minimizers, uint64_t *output_hashes, unsigned int data_size_words) {
    for (unsigned int idx = 0; idx < data_size_words; idx++) {
        uint64_t input_data = input_minimizers[idx];
        uint64_t output_data;
        parallel_hash_calc(input_data, output_data);
        output_hashes[idx] = output_data;
    }
}

int main() {
    std::cout << "=====================================\n";
    std::cout << "  Test CPU : Minimizer Hash \n";
    std::cout << "=====================================\n";

    std::vector<uint64_t> minimizers(NUM_WORDS);
    std::vector<uint64_t> hashes(NUM_WORDS);

    std::cout << " Génération de données aléatoires...\n";
    std::mt19937_64 gen(42);
    std::uniform_int_distribution<uint64_t> dis(0, UINT64_MAX);
    for (auto &v : minimizers) {
        v = dis(gen);
    }

    auto start = std::chrono::high_resolution_clock::now();
    minimizer_cpu(minimizers.data(), hashes.data(), NUM_WORDS);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;

    double time_s = elapsed.count();
    double total_hashes = static_cast<double>(NUM_WORDS);
    double hashes_per_s = total_hashes / time_s;
    double mhash_per_s = hashes_per_s / 1e6;
    double mb_per_s = (DATA_SIZE_BYTES / (1024.0 * 1024.0)) / time_s;

    std::cout << "[CPU] Terminé ✅\n";
    std::cout << "Temps total : " << time_s << " s\n";
    std::cout << "Débit : " << mhash_per_s << " M hash/s\n";
    std::cout << "Vitesse : " << mb_per_s << " MB/s\n";
    return 0;
}

