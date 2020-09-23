
#include "sparsepp/sparsepp/spp.h"
#include <unordered_map>
#include <zlib.h>
#include <fstream>
#include "bam_genes_aux.hpp"
#include "bam_genes.hpp"
#include "tokenizer.hpp"
#include "sbam_merge.hpp"
#include "gzstream.hpp"
#include "htslib/htslib/hts_endian.h"
using namespace gwsc;

int main(int argc, char * argv[]){
    if(argc != 3) return 1;
    std::string bamin = argv[1];
    std::string bamout = argv[2];
    BamReader bin;
    BamDetail read;
    bin.set_bam(bamin);

    bin.set_threads(4);

    samFile * bam_out = sam_open(bamout.c_str(), "wb");
    hts_set_threads(bam_out, 4);
    if(sam_hdr_write(bam_out, bin.header()) < 0) {
        std::cerr << "Error writing header\n";
        exit(1);
    }

    size_t kept = 0, removed = 0;
    while(bin.next(read.b) != nullptr){
        uint32_t * cig = bam_get_cigar(read.b);
        auto bam = read.b;
        bool bad = false;
        for(size_t j = 0; j < bam->core.n_cigar; j+=3){
            CigarElement c1(cig[j]);
            CigarElement c2(cig[j + 1]);
            CigarElement c3(cig[j + 2]);
            if(c1.op == Cigar::REF_SKIP && c2.op == Cigar::DEL && c3.op == Cigar::REF_SKIP){
                bad = true;
            }
        }
        if(bad) removed++;
        else{
            kept++;
            if(sam_write1(bam_out, bin.header(), bam) < 0){
                std::cout << "error writing\n";
                exit(0);
            }
        }
    }
    sam_close(bam_out);
    std::cout << "Done kept = " << kept << " bad = " << removed << "\n";
    return 1;
}
