#include <iostream>
#include <vector>
#include <cstdint>
#include <random>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <omp.h>

#define SMER_SIZE 56
#define WINDOW_SIZE 8

static inline uint64_t mask_right(int numbits) {
    return (numbits >= 64) ? ~0ULL : ((1ULL << numbits) - 1ULL);
}

static inline uint64_t bfc_hash_64(uint64_t key, uint64_t mask) {
    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
}

static inline uint8_t decode_enc_from_packed(const uint64_t* packed, uint64_t idx_base) {
    uint64_t word = packed[idx_base / 8];
    int j = (int)(idx_base % 8);
    uint8_t byte = (uint8_t)((word >> (8 * j)) & 0xFFu);
    return (byte >> 1) & 0x3u;
}

void minimizer_sequence_cpu_omp(const uint64_t* packed_sequence,
                                uint64_t n_bases,
                                std::vector<uint64_t>& out_minimizers,
                                double& t_smer_s,
                                double& t_hash_s,
                                double& t_dedup_s,
                                size_t& n_smers_out)
{
    constexpr int smer_bases = SMER_SIZE / 2;
    const uint64_t HASH_MASK = mask_right(SMER_SIZE);

    out_minimizers.clear();
    if (n_bases == 0 || n_bases < (uint64_t)smer_bases) {
        t_smer_s = t_hash_s = t_dedup_s = 0.0; n_smers_out = 0; return;
    }

    std::vector<uint64_t> vmins;
    vmins.reserve(n_bases);

    uint64_t current_smer = 0;
    uint64_t cur_inv_smer = 0;

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < smer_bases - 1; i++) {
        uint8_t enc = decode_enc_from_packed(packed_sequence, (uint64_t)i);
        current_smer <<= 2;
        current_smer |= (uint64_t)enc;
        cur_inv_smer >>= 2;
        uint64_t comp = (uint64_t)((0x2u ^ enc) & 0x3u);
        cur_inv_smer &= mask_right(SMER_SIZE - 2);
        cur_inv_smer |= (comp << (SMER_SIZE - 2));
    }

    const uint64_t last_pos = n_bases - 1;
    for (uint64_t i = (uint64_t)(smer_bases - 1); i <= last_pos; i++) {
        uint8_t enc = decode_enc_from_packed(packed_sequence, i);
        current_smer <<= 2;
        current_smer |= (uint64_t)enc;
        cur_inv_smer = (cur_inv_smer >> 2) |
                       (((uint64_t)((0x2u ^ enc) & 0x3u)) << (SMER_SIZE - 2));
        uint64_t vmin = (current_smer < cur_inv_smer) ? current_smer : cur_inv_smer;
        vmins.push_back(vmin & HASH_MASK);
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    t_smer_s = std::chrono::duration<double>(t1 - t0).count();
    n_smers_out = vmins.size();

    std::vector<uint64_t> vhashes(vmins.size());
    auto t2 = std::chrono::high_resolution_clock::now();
    #pragma omp parallel for schedule(static)
    for (ptrdiff_t i = 0; i < (ptrdiff_t)vmins.size(); i++) {
        vhashes[(size_t)i] = bfc_hash_64(vmins[(size_t)i], HASH_MASK);
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    t_hash_s = std::chrono::duration<double>(t3 - t2).count();

    if (vhashes.empty()) { t_dedup_s = 0.0; return; }

    uint64_t buffer[WINDOW_SIZE] = {0};
    size_t idx = 0;
    size_t preload = std::min((size_t)WINDOW_SIZE, vhashes.size());
    for (size_t p = 0; p < preload; p++) {
        buffer[p] = vhashes[idx++];
    }

    uint64_t lastElement = (uint64_t)(-1);
    auto t4 = std::chrono::high_resolution_clock::now();
    while (idx < vhashes.size()) {
        uint64_t vhash = vhashes[idx++];
        uint64_t minz = vhash;
        for (size_t p = 0; p < preload; p++) {
            minz = (buffer[p] < minz) ? buffer[p] : minz;
        }
        for (size_t p = 0; p + 1 < preload; p++) buffer[p] = buffer[p + 1];
        if (preload < (size_t)WINDOW_SIZE) buffer[preload++] = vhash; else buffer[WINDOW_SIZE - 1] = vhash;
        if (lastElement != minz) { out_minimizers.push_back(minz); lastElement = minz; }
    }
    auto t5 = std::chrono::high_resolution_clock::now();
    t_dedup_s = std::chrono::duration<double>(t5 - t4).count();
}

static std::vector<uint64_t> make_random_packed_sequence(uint64_t n_bases, uint32_t seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> enc_dis(0, 3);
    const uint64_t n_words = (n_bases + 7) / 8;
    std::vector<uint64_t> data(n_words, 0);
    for (uint64_t i = 0; i < n_bases; i++) {
        uint8_t enc = (uint8_t)enc_dis(rng);
        uint8_t byte = (uint8_t)(enc << 1); 
        uint64_t &w = data[i / 8];
        int j = (int)(i % 8);
        w &= ~((uint64_t)0xFF << (8 * j));
        w |= ((uint64_t)byte) << (8 * j);
    }
    return data;
}

int main(int argc, char** argv) {
    uint64_t n_bases = 32ull * 1024 * 1024; // 32M bases
    int threads = 0; 
    if (argc >= 2) {
        n_bases = std::strtoull(argv[1], nullptr, 10);
        if (n_bases < (SMER_SIZE/2)) n_bases = (SMER_SIZE/2);
    }
    if (argc >= 3) {
        threads = std::atoi(argv[2]);
    }

    if (threads > 0) omp_set_num_threads(threads);
    int max_threads = omp_get_max_threads();

    std::cout << "=====================================\n";
    std::cout << "  CPU Minimizer OMP (sequence -> hash)\n";
    std::cout << "  Bases   : " << n_bases << "\n";
    std::cout << "  Threads : " << max_threads << "\n";
    std::cout << "=====================================\n";

    auto packed = make_random_packed_sequence(n_bases, 42);

    std::vector<uint64_t> out;
    double t_smer=0.0, t_hash=0.0, t_dedup=0.0; size_t n_smers=0;
    auto t_all0 = std::chrono::high_resolution_clock::now();
    minimizer_sequence_cpu_omp(packed.data(), n_bases, out, t_smer, t_hash, t_dedup, n_smers);
    auto t_all1 = std::chrono::high_resolution_clock::now();

    double dt_s = std::chrono::duration<double>(t_all1 - t_all0).count();
    double hashes = (double)out.size();
    double mhash_s = hashes / 1e6 / dt_s;
    double mbin_s = (n_bases /* ~ bytes */) / (1024.0 * 1024.0) / dt_s;

 

    uint64_t smer_bases = SMER_SIZE/2;
    uint64_t n_smer_expected = (n_bases >= smer_bases) ? (n_bases - smer_bases + 1) : 0;
    std::cout << "S-mers       : " << n_smers << "\n";
    std::cout << "Minimizers   : " << out.size() << "\n";
    std::cout << "t_smer       : " << t_smer << " s\n";
    std::cout << "t_hash       : " << t_hash << " s\n";
    std::cout << "t_dedup      : " << t_dedup << " s\n";

    std::cout << "Temps        : " << dt_s << " s\n";
    std::cout << "Debit        : " << mhash_s << " Mhash/s\n";
    std::cout << "Vitesse      : " << mbin_s << " MB/s \n";

    if (!out.empty()) {
        std::cout << "Exemples: ";
        for (size_t i = 0; i < std::min<size_t>(5, out.size()); i++) {
            std::cout << std::hex << out[i] << (i+1<5?" ":"\n") << std::dec;
        }
    }

    int nTests = 64;
    std::vector<double> t_total_v; t_total_v.reserve(nTests);
    std::vector<double> t_smer_v;  t_smer_v.reserve(nTests);
    std::vector<double> t_hash_v;  t_hash_v.reserve(nTests);
    std::vector<double> t_dedup_v; t_dedup_v.reserve(nTests);

    for (int tt = 0; tt < nTests; ++tt) {
        std::vector<uint64_t> out_i;
        double ts=0, th=0, td=0; size_t nsm=0;
        auto tr0 = std::chrono::high_resolution_clock::now();
        minimizer_sequence_cpu_omp(packed.data(), n_bases, out_i, ts, th, td, nsm);
        auto tr1 = std::chrono::high_resolution_clock::now();
        double ttot = std::chrono::duration<double>(tr1 - tr0).count();
        t_total_v.push_back(ttot);
        t_smer_v.push_back(ts);
        t_hash_v.push_back(th);
        t_dedup_v.push_back(td);
        if (!out_i.empty()) packed[0] ^= out_i[out_i.size()/2];
    }

    auto stats = [](std::vector<double>& v){
        std::sort(v.begin(), v.end());
        double vmin=v.front(), vmax=v.back();
        double vmed=v[v.size()/2];
        double vmean = std::accumulate(v.begin(), v.end(), 0.0)/v.size();
        return std::tuple<double,double,double,double>(vmin,vmed,vmax,vmean);
    };

    auto [tot_min, tot_med, tot_max, tot_mean]   = stats(t_total_v);
    auto [smer_min, smer_med, smer_max, smer_mean] = stats(t_smer_v);
    auto [hash_min, hash_med, hash_max, hash_mean] = stats(t_hash_v);
    auto [ded_min,  ded_med,  ded_max,  ded_mean ] = stats(t_dedup_v);

    std::cout << "\n===== Bench multi-tests (" << nTests << ") =====\n";
    std::cout << "Total : min=" << tot_min <<  " s  med=" << tot_med <<  " s  max=" << tot_max << " s   mean=" << tot_mean << " s\n";
    std::cout << "s-mer : min=" << smer_min << " s  med=" << smer_med << " s  max=" << smer_max << " s  mean=" << smer_mean << " s\n";
    std::cout << "hash  : min=" << hash_min << " s  med=" << hash_med << " s  max=" << hash_max << " s  mean=" << hash_mean << " s\n";
    std::cout << "dedup : min=" << ded_min  << " s  med=" << ded_med  << " s  max=" << ded_max  << " s  mean=" << ded_mean  << " s\n";

    return 0;
}
