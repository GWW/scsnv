#pragma once
/*
Copyright (c) 2018-2019 Gavin W. Wilson

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

#include "base.hpp"
#include <exception>
#include "../muts/bam_genes.hpp"
#include "../util/fasta.hpp"
#include "../muts/pileup_worker.hpp"

namespace gwsc{

class ProgPileup : public ProgBase {
    public:
        argagg::parser parser() const {
            argagg::parser argparser {{
                { "help", {"-h", "--help"},
                  "shows this help message", 0},
                { "cellranger", {"-c", "--cellranger"},
                  "Indicates the merged or collapsed bam file is from cell ranger", 0},
                { "txidx", {"-i", "--index"},
                  "Transcript Index", 1},
                { "reference", {"-r", "--ref"},
                  "Reference Genome File", 1},
                { "out", {"-o", "--output"},
                  "Output file prefix, ie. mutations.txt", 1},
                { "threads", {"-t", "--threads"},
                  "Number of processor threads (Default 1)", 1},
                { "barcodes", {"-b", "--barcodes"},
                  "Barcode counts file", 1},
                //{ "passed", {"-p", "--passed"},
                // "Only process reads from the barcodes in this file", 1},
                { "dups", {"-d", "--dups"},
                  "Count PCR duplicates", 0},
                //{ "mincov", {"--min-coverage"},
                //  "Minimum read coverage to count as a SNP and for overall base rate calculations (Default 20)", 1},
                { "minbarcodes", {"--min-barcodes"},
                  "Minimum number of unique barcodes to examine position and to use for base rate calculations (Default 15)", 1},
                { "minalt", {"--min-alt"},
                  "Minimum alternative unique barcodes to count as a SNP (Default 10)", 1},
                { "minqual", {"--min-qual"},
                  "Minimum base quality to count position (Default 20)", 1},
                { "minaf", {"--min-af"},
                  "Minimum allele fraction to count as SNP (Default 0.05)", 1},
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
        std::vector<BarcodeRate> barcode_rates_;
        struct BBout{
            std::vector<uint32_t> barcodes;
            std::vector<uint32_t> snps;
            std::vector<uint16_t> plus;
            std::vector<uint16_t> minus;
        };


        template <typename T>
        int run_wrap_();

        template <typename TL>
        void write_h5f_(const std::string & out, TL & bc, std::vector<BBout> & bouts);

        Fastas                        genome_;
        std::vector<std::string>      bam_files_;
        std::vector<PositionCoverage> coverage_;

        std::string      index_;
        std::string      lib_type_;
        std::string      bc_counts_;
        std::string      out_;
        std::string      ref_;


        //unsigned int min_coverage_ = 20;
        unsigned int min_alternative_ = 10;
        unsigned int min_qual_ = 20;
        unsigned int min_barcodes_ = 15;
        double min_af_ = 0.05;

        unsigned int     threads_ = 1;
        bool             dups_ = false;
        bool             cellranger_ = false;

};

}
