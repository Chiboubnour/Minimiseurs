/*
    Entrée : packed_sequence[] (64 bits / mot)
	Chaque mot contient 8 nucléotides encodés sur 1 octet (8 bits/base)
	m_axi → FPGA

        ┌─────────────────────────────┐
        │ thread_reader_v2             │
        │ - Décode chaque base         │
        │ - Encodage 2 bits + validité │
        │ - 8 bases → 24 bits          │
        └───────────┬─────────────────┘
                    │ Flux : ap_uint<24>
                    │ (8 bases encodées en 3 bits chacune)
                    ▼
        ┌─────────────────────────────┐
        │ thread_smer_v2               │
        │ - Construit s-mer (forward)  │
        │ - Construit reverse complement│
        │ - Prend la plus petite (canonique)
        │ - Calcule hash  → smer_size bits
        └───────────┬─────────────────┘
                    │ Flux : ap_uint<smer_size>
                    │ (smer_size = 2 × taille_smer_en_bases)
                    ▼
        ┌─────────────────────────────┐
        │ thread_dedup_v2              │
        │ - Fenêtre glissante (kmer-s) │
        │ - Trouve minimizer           │
        │ - Évite doublons             │
        └───────────┬─────────────────┘
                    │ Flux : ap_uint<smer_size>
                    ▼
        ┌─────────────────────────────┐
        │ thread_store_v2              │
        │ - Stocke hash en mémoire     │
        │ - Compte nMinimizers         │
        └───────────┬─────────────────┘
                    │
                    ▼
       tab_hash[] (64 bits / hash)
       nMinizrs (64 bits)
       m_axi → DDR


 #include <ap_int.h>
 #include <hls_stream.h>
 #include <cstdint>

 #define SMER_SIZE 56
 #define WINDOW_SIZE 8
 #define DATA_DEPTH 1024
 #define MEM_UNIT 64

 inline ap_uint<2> nucl_encode(char nucl) {
     #pragma HLS INLINE
     switch (nucl) {
         case 'A': return 0;
         case 'C': return 1;
         case 'G': return 2;
         case 'T': return 3;
         default : return 0;
     }
 }

 inline ap_uint<64> min_v1(const ap_uint<64> a, const ap_uint<64> b) {
     #pragma HLS INLINE
     return (a < b) ? a : b;
 }

 inline ap_uint<64> mask_right(int numbits) {
     #pragma HLS INLINE
     return (numbits >= MEM_UNIT) ? ~0ULL : ((1ULL << numbits) - 1ULL);
 }

 inline ap_uint<64> hash_u64(ap_uint<64> key, ap_uint<64> mask) {
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

 void thread_reader_v2(
     const ap_uint<64>* packed_sequence,
     ap_uint<64> n_bases,
     hls::stream< ap_uint<24> >& stream_o
 ) {
     const int n_words = (int)((n_bases + 7) / 8);

     for (int i = 0; i < n_words; ++i) {
     #pragma HLS PIPELINE II=1
         const ap_uint<64> word_8b = packed_sequence[i];

         ap_uint<24> word_3b = 0;
         bool all_valid = true;

         for (int j = 0; j < 8; ++j) {
         #pragma HLS UNROLL
             const int idx = i * 8 + j;

             ap_uint<1> valid = (idx < (int)n_bases) ? ap_uint<1>(1) : ap_uint<1>(0);
             ap_uint<8> c     = valid ? word_8b.range(8*(j+1)-1, 8*j) : ap_uint<8>(0);

             const ap_uint<2> enc = (c >> 1) & 0x3;
             //const ap_uint<2> enc = nucl_encode(c);
             word_3b.range(3*j+1, 3*j) = enc;  // 2 bits
             word_3b[3*j+2]            = valid;

             all_valid &= (bool)valid;
         }

         stream_o.write(word_3b);

         // Si ce mot contient des invalids, on s'arrête
         if (!all_valid) break;
     }

     if ((n_bases & 7) == 0) {
     #pragma HLS PIPELINE II=1
         stream_o.write( (ap_uint<24>)0 );
     }
 }

 void thread_smer_v2(
     hls::stream< ap_uint<24> >& stream_i,
     ap_uint<64> n_bases,
     hls::stream< ap_uint<SMER_SIZE> >& stream_o
 ) {
     constexpr int smer = SMER_SIZE / 2;
     const ap_uint<64> HASH_MASK = mask_right(SMER_SIZE);

     ap_uint<SMER_SIZE> current_smer = 0;
     ap_uint<SMER_SIZE> cur_inv_smer = 0;
     ap_uint<24> word_24b = 0;

     for (int i = 0; i < smer - 1; i++) {
     #pragma HLS PIPELINE II=1
         if ((i & 0x07) == 0){
             word_24b = stream_i.read();
         }
         const ap_uint<2> c_nucl = word_24b.range(1,0);

         current_smer <<= 2;
         current_smer(1,0) = c_nucl;

         cur_inv_smer >>= 2;
         cur_inv_smer(SMER_SIZE-1, SMER_SIZE-2) = (0x2 ^ c_nucl);

         word_24b >>= 3;
     }

     const int last_pos = (int)n_bases - 1;
     for (int i = smer - 1; i <= last_pos; i++) {
     #pragma HLS PIPELINE II=1
         if ((i & 0x07) == 0) word_24b = stream_i.read();

         const ap_uint<2> c_nucl = word_24b.range(1,0);
         const bool valid = word_24b[2];

         current_smer <<= 2;
         current_smer(1,0) = c_nucl;

         cur_inv_smer = (cur_inv_smer >> 2) | ((ap_uint<SMER_SIZE>)((0x2 ^ c_nucl)) << (SMER_SIZE-2));

         const ap_uint<SMER_SIZE> vmin  = min_v1(current_smer, cur_inv_smer);
         const ap_uint<SMER_SIZE> vhash = hash_u64(vmin, HASH_MASK);

         word_24b >>= 3;

         if (valid) {
             stream_o.write(vhash);
         } else {
             stream_o.write(0x00);
             return;
         }
     }

     stream_o.write(0x00);
 }

 void thread_dedup_v2(
     hls::stream< ap_uint<SMER_SIZE> >& stream_i,
     hls::stream< ap_uint<SMER_SIZE> >& stream_o
 ) {
     ap_uint<SMER_SIZE> buffer[WINDOW_SIZE];
     for (int p = 0; p < WINDOW_SIZE; p++) {
     #pragma HLS PIPELINE II=1
         buffer[p] = stream_i.read();
     }

     ap_uint<SMER_SIZE> lastElement = (ap_uint<SMER_SIZE>)(-1);

     while (true) {
     #pragma HLS PIPELINE II=1
         const ap_uint<SMER_SIZE> vhash = stream_i.read();
         if (vhash == 0) { stream_o.write(0x00); break; }

         ap_uint<SMER_SIZE> minz = vhash;
         for (int p = 0; p < WINDOW_SIZE; p++) {   //WINDOW_SIZE comparaisons à chaque cycle !!!
         #pragma HLS UNROLL
             minz = (buffer[p] < minz) ? buffer[p] : minz;
         }
         for (int p = 0; p < WINDOW_SIZE - 1; p++) {
         #pragma HLS UNROLL
             buffer[p] = buffer[p+1];
         }
         buffer[WINDOW_SIZE-1] = vhash;

         if (lastElement != minz) {
             stream_o.write(minz);
             lastElement = minz;
         }
     }
 }

 void thread_store_v2(
     hls::stream< ap_uint<SMER_SIZE> >& stream_i,
     ap_uint<64>* tab_hash,
     ap_uint<64>* nElements
 ) {
     int cnt = 0;
     while (true) {
     #pragma HLS PIPELINE II=1
         const ap_uint<SMER_SIZE> vHash = stream_i.read();
         if (vHash == 0) break;
         tab_hash[cnt++] = vHash;
     }
     *nElements = cnt;
 }

 #define KMER 32
 #define SMER 28

 void minimizer(
     const ap_uint<64>* packed_sequence,
     ap_uint<64> n,
     ap_uint<64>* tab_hash,
     ap_uint<64>* nMinizrs
 ) {
     #pragma HLS INTERFACE mode=m_axi port=packed_sequence
     #pragma HLS INTERFACE mode=m_axi port=tab_hash
     #pragma HLS INTERFACE mode=s_axilite port=nMinizrs
     #pragma HLS INTERFACE  mode=s_axilite port=n
     #pragma HLS INTERFACE mode=s_axilite port=return bundle=control
     #pragma HLS DATAFLOW

     constexpr int z = KMER - SMER;

     hls::stream< ap_uint<24>, DATA_DEPTH > fifo_1;
     hls::stream< ap_uint<SMER_SIZE>, DATA_DEPTH > fifo_2;
     hls::stream< ap_uint<SMER_SIZE>, DATA_DEPTH > fifo_3;

     thread_reader_v2(packed_sequence, n, fifo_1);
     thread_smer_v2(fifo_1, n, fifo_2);
     thread_dedup_v2(fifo_2, fifo_3);
     thread_store_v2(fifo_3, tab_hash, nMinizrs);
 }
*/
