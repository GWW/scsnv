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
#include "fasta.hpp"

namespace gwsc{

class ProgCollapse : public ProgBase {
    public:
        argagg::parser parser() const {
            argagg::parser argparser {{
                { "help", {"-h", "--help"},
                  "shows this help message", 0},
                { "txidx", {"-i", "--index"},
                  "Transcript Index", 1},
                { "reference", {"-r", "--ref"},
                  "Reference Genome File", 1},
                { "out", {"-o", "--output"},
                  "Output file prefix, ie. sample/merged", 1},
                { "barcodes", {"-b", "--barcodes"},
                  "Raw Barcode counts file", 1},
                { "threads", {"-t", "--threads"},
                  "Number of processor threads (Default 1)", 1},
                { "reads", {"-m", "--reads"},
                  "Maximum reads to read from a gene (in millions) (Default 10)", 1},
                { "bam_write", {"-w", "--bam-write"},
                  "Number of writer threads to use when emitting sorted bam files (Default 1)", 1},
                { "library", {"-l", "--library"},
                  "libary type (V2)", 1},
              }};
            return argparser;
        }

        std::string usage() const {
            return "scsnv collapse -r <genome.fa> -i <transcript index prefix> -o <out prefix> -b <barcodes> <tmp bam prefix> <tmp bam prefix 2> ... <tmp bam prefix N>";
        }

        void load();

        int run();

    private:
        template <typename T>
        int run_wrap_();

        std::string      bam_file_;
        Fastas           genome_;

        std::string      index_;
        std::string      umi_map_;
        std::string      lib_type_;
        std::string      bc_counts_;
        std::string      out_;
        std::string      ref_;
        uint64_t         max_reads_ = 10000000;
        unsigned int     bam_write_threads_ = 1;
        unsigned int     threads_ = 1;

};

}
