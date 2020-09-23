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
#include "align_aux.hpp"
#include "index.hpp"
#include "genome_align.hpp"
#include "bwa/bwamem.h"

namespace gwsc {

class TranscriptAlign{
    TranscriptAlign( const TranscriptAlign& ) = delete;
    TranscriptAlign& operator=(const TranscriptAlign&) = delete;
    public:
        TranscriptAlign(const TXIndex & index, StrandMode smode) : idx_(index), smode_(smode) {

        }

        ~TranscriptAlign(){
            if(bidx_ != nullptr)
                bwa_idx_destroy(bidx_);
            if(args_ != nullptr)
                free(args_);
        }

        void unload(){
            if(bidx_ != nullptr)
                bwa_idx_destroy(bidx_);
            if(args_ != nullptr)
                free(args_);

            bidx_ = nullptr;
            args_ = nullptr;
        }

        void load(const std::string & prefix, unsigned int min_overhang);
        void set_gidx(const GenomeAlign * gidx) {
            gidx_ = gidx;
        }
        void align(AlignGroup & ad, const std::string & seq, unsigned int len) const;
        void get(AlignGroup & ad, size_t i, const std::string & seq, unsigned int len) const;
        void project(AlignData & a, CigarString & tc) const;

        AlignScore ascore;
        


    private:
        void make_align_(AlignGroup & ad, mem_aln_t & a) const;
        void trim_splices_(AlignData & a) const;

        const TXIndex              & idx_;
        bwaidx_t                   * bidx_  = nullptr;
        const GenomeAlign         * gidx_  = nullptr;
        mem_opt_t                  * args_ = nullptr;
        unsigned int                 min_overhang_;
        StrandMode                   smode_;
};

}
