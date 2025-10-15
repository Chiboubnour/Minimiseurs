#include <ap_int.h>
#include <hls_stream.h>
#include <cstdint>

#define HASH_WIDTH 64

inline ap_uint<64> bfc_hash_64(ap_uint<64> key, ap_uint<64> mask) {
    #pragma HLS INLINE
    key = (~key + (key << 21)) & mask;
    key = key ^ (key >> 24);
    key = ((key + (key << 3)) + (key << 8)) & mask;
    key = key ^ (key >> 14);
    key = ((key + (key << 2)) + (key << 4)) & mask;
    key = key ^ (key >> 28);
    key = (key + (key << 31)) & mask;
    return key;
}

void parallel_hash_calc(ap_uint<64> input_data, ap_uint<64>& output_data) {
    #pragma HLS INLINE
    ap_uint<64> mask = ~0ULL;
    output_data = bfc_hash_64(input_data, mask);
}

void minimizer(
    const ap_uint<64>* input_minimizers,
    ap_uint<64>* output_hashes,
    unsigned int data_size_words
) {
    #pragma HLS INTERFACE m_axi port=input_minimizers offset=slave bundle=gmem_in max_read_burst_length=16
    #pragma HLS INTERFACE m_axi port=output_hashes offset=slave bundle=gmem_out max_write_burst_length=16
    #pragma HLS INTERFACE s_axilite port=data_size_words
    #pragma HLS INTERFACE s_axilite port=return

    const unsigned int BURST_SIZE = 16;

    for (unsigned int burst_idx = 0; burst_idx < data_size_words; burst_idx += BURST_SIZE) {
        unsigned int remaining = data_size_words - burst_idx;
        unsigned int current_burst = (remaining < BURST_SIZE) ? remaining : BURST_SIZE;

        for (unsigned int word_idx = 0; word_idx < current_burst; word_idx++) {
            #pragma HLS PIPELINE II=1

            unsigned int idx = burst_idx + word_idx;
            ap_uint<64> input_data = input_minimizers[idx];
            ap_uint<64> output_data;

            parallel_hash_calc(input_data, output_data);
            output_hashes[idx] = output_data;
        }
    }
}

