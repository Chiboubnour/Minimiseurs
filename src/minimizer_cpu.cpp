#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
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
    std::cout << "  Test CPU : Minimizer Hash\n";
    std::cout << "=====================================\n";

    int nTests = 64;

    uint64_t nElements = (32ull * 1024 * 1024) / sizeof(uint64_t);

    std::vector<uint64_t> minimizers(nElements);
    std::vector<uint64_t> hashes(nElements);

    std::mt19937_64 gen(42);
    std::uniform_int_distribution<uint64_t> dis(0, UINT64_MAX);
    for (auto &v : minimizers) {
        v = dis(gen);
    }

    std::vector<double> durations_us;
    durations_us.reserve(nTests);

    for (int tt = 0; tt < nTests; tt++) {
        auto start = std::chrono::high_resolution_clock::now();
        
        minimizer_cpu(minimizers.data(), hashes.data(), nElements);
        minimizers[rand() % nElements] = hashes[rand() % nElements];
        
        auto end = std::chrono::high_resolution_clock::now();
        double time_us = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
        durations_us.push_back(time_us);
    }

    std::sort(durations_us.begin(), durations_us.end());
    double time_min = durations_us.front() / 1e6;       
    double time_max = durations_us.back() / 1e6;
    double time_median = durations_us[durations_us.size() / 2] / 1e6;
    double time_mean = (std::accumulate(durations_us.begin(), durations_us.end(), 0.0) / durations_us.size()) / 1e6;

    double total_hashes = nElements;
    double hashes_per_s = total_hashes / time_median;
    double mhash_per_s  = hashes_per_s / 1e6;
    double mb_per_s     = (nElements * sizeof(uint64_t)) / (1024.0 * 1024.0) / time_median;

    std::cout << "Terminé ✅  "
              << " | t_min : " << time_min << " s"
              << " | t_med : " << time_median << " s"
              << " | t_max : " << time_max << " s"
              << " | moy : " << time_mean << " s"
              << " | Débit : " << mhash_per_s << " M hash/s"
              << " | Vitesse : " << mb_per_s << " MB/s"
              << std::endl;

    return 0;
}



/*
data (3:0) <= read(3:0)
write(3:0) <= data(3:0)

data (3:0) <=      read(3:0)
write(3:0) <= ? ? ? data(0)
data (2:0) <=      data(3:1)

data (6:3) <= read(3:0)
write(3:0) <= data(3:0)
data (2:0) <= data(6:4)

data (6:3) <= read(3:0)
write(3:0) <= data(3:0)
data (2:0) <= data(6:4)

data (6:3) <= read(3:0)
write(3:0) <= data(3:0)
data (2:0) <= data(6:4)

data (6:3) <= read(3:0)
write(3:0) <= data(3:0)
data (2:0) <= data(6:4)

data (6:3) <= read(3:0)
write(3:0) <= data(3:0)
data (2:0) <= data(6:4)
*/