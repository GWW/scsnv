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
#include "index.hpp"
#include "tokenizer.hpp"
#include "fasta.hpp"
#include "pileup_aux.hpp"
#include "parallel_hashmap/parallel_hashmap/phmap.h"
#include <exception>

namespace gwsc{

class ProgAccuracy : public ProgBase {
    public:
        argagg::parser parser() const;
        std::string usage() const {
            return "scsnv accuracy -g genome_bam.bam -b passed_barcodes.txt -o out_file scsnv_bam.bam cellranger_bam.bam";
        }

        int run();
        void load();

    private:
        template <typename B>
        int run_();
        void read_passed_(unsigned int blength);

        phmap::flat_hash_map<std::string, unsigned int>    bchash_;
        TXIndex                  idx_;
        std::string              lib_;
        std::string              bcin_;
        std::string              outp_;
        std::string              isnvs_;
        std::string              iref_;
        std::string              iprefix_;
        std::string              iname_;
        std::string              bam_;
        std::string              aligners_;
        std::vector<std::string> barcodes_;
        Fastas                   genome_;
        unsigned int             threads_;
        unsigned int             window_;
};


}
