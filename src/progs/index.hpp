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
#include <exception>

namespace gwsc{

class ProgIndex : public ProgBase {
    public:
        argagg::parser parser() const {
            argagg::parser argparser {{
                { "gtf", {"-g", "--gtf"},
                  "GTF file (required)", 1},
                { "ref", {"-r", "--ref"},
                  "Genome fasta file (required)", 1},
                { "introns", {"--retained-introns"},
                  "Keep retained introns, note alignments to these will be considered exonic", 0},
                { "bwabuild", {"--skip-build"},
                  "Skip automatically building the BWA index (not recommended)", 0},
                { "tlen", {"-l", "--min-length"},
                  "Minimum Transcript Length [100]", 1},
                { "help", {"-h", "--help"},
                  "shows this help message", 0},
              }};
            return argparser;
        }

        std::string usage() const {
            return "scsnv index -g genes.gtf -r genome.fa -l 100 out_prefix";
        }

        void load() {
            gtf_ = args_["gtf"].as<std::string>();
            ref_ = args_["ref"].as<std::string>();
            min_length_ = args_["tlen"].as<unsigned int>(min_length_);
            retained_introns_ = args_["introns"];
            bwa_build_ = !args_["bwabuild"];
            if(args_.pos.size() != 1){
                throw std::runtime_error("Missing the output option");
            }
            out_ = args_.as<std::string>(0);
        }

        int run() {
            TXIndexBuild idx(out_, bwa_build_);
            idx.build(ref_, gtf_, min_length_, retained_introns_);
            return EXIT_SUCCESS;
        }

    private:
        std::string gtf_;
        std::string ref_;
        std::string out_;
        unsigned int min_length_ = 100;
        bool         retained_introns_ = false;
        bool         bwa_build_ = true;
};

}
