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
#include "../scsnv/build.hpp"
#include "../scsnv/index.hpp"
#include <exception>

namespace gwsc{

class ProgQuant : public ProgBase {
    public:
        argagg::parser parser() const {
            argagg::parser argparser {{
                { "help", {"-h", "--help"},
                  "shows this help message", 0},
                { "txidx", {"-i", "--index"},
                  "Transcript Index", 1},
                { "threads", {"-t", "--threads"},
                  "Number of threads", 1},
                { "no_bam", {"--no-bam"},
                  "Disable Writing data required to correct the UMI of bam files", 0},
                { "cgroups", {"-c", "--count-groups"},
                  "Gene Groups for cell quantification", 1},
                { "molecules", {"-m", "--min-molecules"},
                  "Minimum number of barcode cDNA molecules", 1},
                { "library", {"-l", "--library"},
                  "libary type (V2)", 1},
                { "output", {"-o", "--output"},
                  "Output file prefix, ie. sample/summary", 1},
              }};
            return argparser;
        }

        std::string usage()  const {
            return "scsnv quant -i <transcript index prefix> -o <out prefix> <map prefix 1> ... <map prefix N>";
        }

        void load() {
            tx_idx_ = args_["txidx"].as<std::string>();
            out_prefix_ = args_["output"].as<std::string>();
            threads_ = args_["threads"].as<unsigned int>(1);
            min_molecules_ = args_["molecules"].as<unsigned int>(1);
            lib_type_ = args_["library"].as<std::string>("V2");
            gene_groups_ = args_["cgroups"].as<std::string>("");
            bam_ = !args_["no_bam"];

            if(args_.pos.size() == 0){
                throw std::runtime_error("Missing mapped file prefix(es)");
            }
            for(auto & f : args_.pos){
                std::string d = f;
                while(d.rbegin() != d.rend() && *d.rbegin() == '\\') d.pop_back();
                prefixes_.push_back(d);
            }
            std::sort(prefixes_.begin(), prefixes_.end(), ReadStringCmp());
        }

        int run();

    private:
        std::string              tx_idx_;
        std::string              gene_groups_;
        std::string              out_prefix_;
        std::string              lib_type_;
        std::vector<std::string> prefixes_;
        FastqPairs               fastqs_;
        TXIndex                  txi_;
        unsigned int             threads_;
        unsigned int             min_molecules_;
        bool                     bam_ = false;
};

}
