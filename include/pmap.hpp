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
#include "reader.hpp"
#include <exception>

namespace gwsc{

class ProgMap : public ProgBase {
    public:
        argagg::parser parser() const;
        std::string usage() const;
        void load();
        int run();

    private:
        template <typename T>
        int run_wrap_();

        void write_progress_();

        std::string              tx_idx_;
        std::string              genome_idx_;
        std::string              out_prefix_;
        std::string              bc_counts_;
        std::string              lib_type_;
        std::string              tmp_bam_;
        std::string              gene_groups_;
        std::vector<std::string> dirs_;
        FastqPairs               fastqs_;
        TXIndex                  txi_;
        double                   dust_;
        size_t                   written_ = 0;
        size_t                   lwritten_ = 0;
        size_t                   total_ = 0;
        size_t                   start_ = 0;
        size_t                   downsample_ = 0;
        size_t                   seed_ = 0;
        unsigned int             threads_;
        unsigned int             qthreads_;
        unsigned int             min_overhang_;
        unsigned int             bam_per_thread_;
        unsigned int             bam_per_file_;
        unsigned int             bam_write_threads_;
        bool                     bam_;
        //bool                     write_tags_;
        bool                     internal_ = false;
};

}
