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
#include "bwa/bwamem.h"

namespace gwsc {

class GenomeAlign{
    GenomeAlign( const GenomeAlign& ) = delete;
    GenomeAlign& operator=(const GenomeAlign&) = delete;
    public:
        GenomeAlign(const TXIndex & index, StrandMode smode) : idx_(index), smode_(smode) {

        }

        ~GenomeAlign(){
            if(bidx_ != nullptr)
                bwa_idx_destroy(bidx_);
            if(args_ != nullptr)
                free(args_);
        }
        void load(const std::string & prefix, unsigned int min_overhang);
        void align(AlignGroup & ad, const std::string & seq, unsigned int len) const;
        void get(AlignGroup & ad, size_t i, const std::string & seq, unsigned int len) const;
        void rescore(AlignData & a, const std::string & seq) const;
        void verify(AlignData & a, const std::string & seq, const std::string & msg = "") const;
        const bwaidx_t * gidx() const {
            return bidx_;
        }

        void unload(){
            if(bidx_ != nullptr)
                bwa_idx_destroy(bidx_);
            if(args_ != nullptr)
                free(args_);

            bidx_ = nullptr;
            args_ = nullptr;
        }

        AlignScore ascore;

    private:
        void annotate_(AlignData & a, TXIndex::Ref::itree::intervalVector & overlaps, std::vector<AlignBases> & bases) const;
        void make_align_(AlignGroup & ad, mem_aln_t & a) const;
        void calculate_bases_(AlignData & a, const GeneEntry & g, AlignBases & b) const;
        void trim_splices_(AlignData & a, const TXIndex::Ref & ref,
                TXIndex::Ref::itree::intervalVector & overlaps) const;

        const TXIndex              & idx_;
        bwaidx_t                   * bidx_  = nullptr;
        mem_opt_t                  * args_ = nullptr;
        unsigned int                 min_overhang_;
        StrandMode                   smode_;
};


}
