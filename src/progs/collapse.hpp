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

namespace gwsc{

class ProgCollapse : public ProgBase {
    public:
        argagg::parser parser() const {
            argagg::parser argparser {{
                { "help", {"-h", "--help"},
                  "shows this help message", 0},
                { "cellranger", {"-c", "--cellranger"},
                  "Bam file is a Cell Ranger pos sorted BAM file only a single positional argument is accepted and it must point to the cell ranger bam file", 0},
                { "txidx", {"-i", "--index"},
                  "Transcript Index", 1},
                { "reference", {"-r", "--ref"},
                  "Reference Genome File", 1},
                { "umi_map", {"-u", "--umi-map"},
                  "Summary UMI map ie. sample/summary_umi_map.gz (not needed for Cell Ranger)", 1},
                { "out", {"-o", "--output"},
                  "Output file prefix, ie. sample/merged", 1},
                { "barcodes", {"-b", "--barcodes"},
                  "Barcode counts file", 1},
                { "threads", {"-t", "--threads"},
                  "Number of processor threads (Default 1)", 1},
                { "bam_write", {"-w", "--bam-write"},
                  "Number of writer threads to use when emitting sorted bam files (Default 1)", 1},
                { "library", {"-l", "--library"},
                  "libary type (V2)", 1},
              }};
            return argparser;
        }

        std::string usage() const {
            return "scmap collapse -u <umi_map> -r <genome.fa> -i <transcript index prefix> -o <out prefix> -b <barcodes> <tmp bam prefix> <tmp bam prefix 2> ... <tmp bam prefix N>";
        }

        void load();

        int run();

    private:
        template <typename T>
        int run_wrap_();

        std::vector<std::string> bam_files_;
        Fastas                   genome_;

        std::string      index_;
        std::string      umi_map_;
        std::string      lib_type_;
        std::string      bc_counts_;
        std::string      out_;
        std::string      ref_;
        unsigned int     bam_write_threads_ = 1;
        unsigned int     threads_ = 1;
        bool             cellranger_ = false;

};

}
