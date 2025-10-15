#include "xparameters.h"
#include "xil_printf.h"
#include "xil_cache.h"
#include "xminimizer.h"
#include "xtime_l.h"
#include <random>
#include <cstdint>
#include <stdio.h>

#define DATA_SIZE_BYTES (256 * 1024 * 1024)
#define NUM_WORDS      (DATA_SIZE_BYTES / sizeof(uint64_t))  

static uint64_t minimizers[NUM_WORDS] __attribute__((aligned(64)));
static uint64_t hashes_hw[NUM_WORDS]  __attribute__((aligned(64)));
static uint64_t hashes_sw[NUM_WORDS]  __attribute__((aligned(64)));

uint64_t bfc_hash_64(uint64_t key, uint64_t mask) {
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

void minimizer_arm(const uint64_t *input_minimizers, uint64_t *output_hashes, unsigned int data_size_words) {
    uint64_t output_data;
    for (unsigned int idx = 0; idx < data_size_words; idx++) {
        uint64_t input_data = input_minimizers[idx];
        parallel_hash_calc(input_data, output_data);
        output_hashes[idx] = output_data;
    }
}

void generate_random_minimizers(uint64_t *data, size_t num_words) {
    std::mt19937_64 gen(42);
    std::uniform_int_distribution<uint64_t> dis(0, UINT64_MAX);
    for (size_t i = 0; i < num_words; i++) {
        data[i] = dis(gen);
    }
}

int main() {
    xil_printf("\n=====================================\n\r");
    xil_printf("   🚀 Test Comparatif FPGA vs ARM\n\r");
    xil_printf("   Build : %s - %s\n\r", __DATE__, __TIME__);
    xil_printf("=====================================\n\r");

    // === Génération des minimizers ===
    xil_printf("[SW] Génération de %d MB de minimizers aléatoires...\n\r", DATA_SIZE_BYTES / (1024 * 1024));
    generate_random_minimizers(minimizers, NUM_WORDS);

    // === Test logiciel (ARM) ===
    xil_printf("\n[SW] Lancement version ARM...\n\r");
    XTime t1, t2;
    XTime_GetTime(&t1);
    minimizer_arm(minimizers, hashes_sw, NUM_WORDS);
    XTime_GetTime(&t2);

    uint64_t cycles_sw = (uint64_t)(t2 - t1);
    double freq_hz = 100e6; 
    double time_sw_s = (double)cycles_sw / freq_hz;
    double sw_hashes_per_s = NUM_WORDS / time_sw_s;
    double sw_hashes_per_ms = sw_hashes_per_s / 1e6;
    double sw_mb_per_s = (DATA_SIZE_BYTES / (1024.0 * 1024.0)) / time_sw_s;

    xil_printf("[SW] Terminé ✅\n\r");
    xil_printf("⏱ cycles : %llu\n\r", (unsigned long long)cycles_sw);
    printf("⏱ temps  : %.3f s\n\r", time_sw_s);
    printf("⚡ Débit  : %.2f Mhash/s\n\r", sw_hashes_per_ms);
    printf("🚀 Vitesse : %.2f MB/s\n\r", sw_mb_per_s);

    // === Test matériel (FPGA) ===
    xil_printf("\n[HW] Initialisation IP Minimizer...\n\r");
    XMinimizer accel;
    XMinimizer_Config *cfg = XMinimizer_LookupConfig(XPAR_MINIMIZER_0_DEVICE_ID);
    if (!cfg || XMinimizer_CfgInitialize(&accel, cfg) != XST_SUCCESS) {
        xil_printf("[ERREUR] Init IP\n\r");
        return XST_FAILURE;
    }
    xil_printf("[IP] OK ✅\n\r");

    // Flush cache
    Xil_DCacheFlushRange((UINTPTR)minimizers, DATA_SIZE_BYTES);
    Xil_DCacheFlushRange((UINTPTR)hashes_hw, DATA_SIZE_BYTES);

    XMinimizer_Set_input_minimizers(&accel, (u64)minimizers);
    XMinimizer_Set_output_hashes(&accel, (u64)hashes_hw);
    XMinimizer_Set_data_size_words(&accel, NUM_WORDS);

    xil_printf("[HW] Lancement accélérateur FPGA...\n\r");
    XTime t3, t4;
    XTime_GetTime(&t3);
    XMinimizer_Start(&accel);
    while (!XMinimizer_IsDone(&accel)) { }
    XTime_GetTime(&t4);

    uint64_t cycles_hw = (uint64_t)(t4 - t3);
    double time_hw_s = (double)cycles_hw / freq_hz;
    double hw_hashes_per_s = NUM_WORDS / time_hw_s;
    double hw_hashes_per_ms = hw_hashes_per_s / 1e6;
    double hw_mb_per_s = (DATA_SIZE_BYTES / (1024.0 * 1024.0)) / time_hw_s;

    xil_printf("[HW] Terminé ✅\n\r");

    // === Résumé comparatif ===
    xil_printf("\n=========== 🔍 COMPARAISON ===========\n\r");
    printf("ARM  : %.3f s | %.2f Mhash/s | %.2f MB/s\n", time_sw_s, sw_hashes_per_ms, sw_mb_per_s);
    printf("FPGA : %.3f s | %.2f Mhash/s | %.2f MB/s\n", time_hw_s, hw_hashes_per_ms, hw_mb_per_s);
    xil_printf("=====================================\n\r");

    xil_printf("\nTest terminé avec succès ✅\n");
    return 0;
}
