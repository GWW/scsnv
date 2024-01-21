#pragma once
/*
Copyright (c) 2018-2020 Gavin W. Wilson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "pbase.hpp"
#include <exception>
#include "bam_genes.hpp"
#include "parallel_hashmap/parallel_hashmap/phmap.h"
#include "fasta.hpp"
#include "pileup_worker.hpp"

namespace gwsc{

class ProgPileup : public ProgBase {
    public:
        using BarcodeHash = phmap::flat_hash_map<std::string, unsigned int>;
        argagg::parser parser() const {
            argagg::parser argparser {{
                { "help", {"-h", "--help"},
                  "shows this help message", 0},
                { "cellranger", {"-c", "--cellranger"},
                  "Indicates the merged bam file is from cell ranger", 0},
                { "txidx", {"-i", "--index"},
                  "Transcript Index", 1},
                { "reference", {"-r", "--ref"},
                  "Reference Genome File", 1},
                { "out", {"-o", "--output"},
                  "Output file prefix, ie. mutations.txt", 1},
                { "threads", {"-t", "--threads"},
                  "Number of processor threads (Default 1)", 1},
                { "rthreads", {"-x", "--read-threads"},
                  "Number of reader threads (Default 1)", 1},
                { "snvlist", {"-s", "--snvs"},
                "Only quantify this list of tab separated SNVs (including reference only bases) with a header chorm, pos (zero based), REF, ALT", 1},
                { "passed", {"-p", "--passed"},
                 "Only process reads from the barcodes in this file", 1},
                { "dups", {"-d", "--dups"},
                  "Count PCR duplicates", 0},
                { "minedge", {"--min-edge"},
                  "Minimum distance at least 1 read is from the end of a read (Default 5)", 1},
                { "splicewin", {"--splice-win"},
                  "Window around a potential SNV to look for the nearest splice donor / acceptor (Default 10)", 1},
                { "minbarcodes", {"--min-barcodes"},
                  "Minimum number of unique barcodes to examine position and to use for base rate calculations (Default 5)", 1},
                { "minalt", {"--min-alt"},
                  "Minimum alternative unique barcodes to count as a SNP (Default 2)", 1},
                { "minqual", {"--min-qual"},
                  "Minimum base quality to count position (Default 20)", 1},
                { "minaf", {"--min-af"},
                  "Minimum allele fraction to count as SNP (Default 0.01)", 1},
                { "genes", {"--read-genes"},
                  "Number of genes to read for paralell processing. Larger values use more memory (Default 500)", 1},
                { "library", {"-l", "--library"},
                  "libary type (V2)", 1},
              }};
            return argparser;
        }

        std::string usage() const {
            return "scsnv pileup -i <transcript index prefix> -r <genome.fa> -b <barcode_counts.txt.gz> -o <output> in.bam";
        }

        void load();

        int run();

    private:

        std::vector<uint32_t>    tids_;
        std::vector<uint32_t>    pos_;
        std::vector<uint8_t>     refs_;
        //std::vector<BarcodeRate> barcode_rates_;
        struct BBout{
            std::vector<uint32_t> barcodes;
            std::vector<uint32_t> snps;
            std::vector<uint16_t> plus;
            std::vector<uint16_t> minus;
        };


        void write_h5f_(const std::string & out, std::vector<BBout> & bouts);

        void read_passed_();
        void parse_targets_();

        template <typename T, typename P>
        int run_wrap_();

        unsigned int filter_func(const std::string & barcode, unsigned int fno){
            (void)fno;
            auto it = bchash_.find(barcode);
            if(it != bchash_.end()){
                return it->second;
            }
            return std::numeric_limits<unsigned int>::max();
        }

        Fastas                        genome_;
        std::string                   bam_file_;
        std::vector<PositionCoverage> coverage_;
        std::vector<std::string>      barcodes_;
        BarcodeHash                   bchash_;
        std::vector<std::string>      bmap_;
        TargetFinder::trefs           targets_;

        std::string      index_;
        std::string      lib_type_;
        std::string      passed_;
        std::string      out_;
        std::string      ref_;
        std::string      snvlist_;


        //unsigned int min_coverage_ = 20;
        unsigned int min_alternative_ = 2;
        unsigned int min_qual_ = 20;
        unsigned int min_barcodes_ = 5;
        unsigned int min_edge_ = 5;
        unsigned int splice_win_ = 10;
        unsigned int read_genes_ = 500;
        double       min_af_ = 0.01;

        unsigned int     threads_ = 1;
        unsigned int     rthreads_ = 1;
        bool             dups_ = false;
        bool             cellranger_ = false;

};
}

