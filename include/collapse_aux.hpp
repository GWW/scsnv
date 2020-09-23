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

#include <vector>
#include <utility>
#include <numeric>
#include <algorithm>
#include <iostream>
#include "htslib/htslib/hts.h"
#include "htslib/htslib/sam.h"
#include "sequence.hpp"
#include "fasta.hpp"
#include "bam_genes.hpp"
#include "aux.hpp"
#include "align_aux.hpp"
#include <mutex>

namespace gwsc {

struct CollapseInsertion {
    CollapseInsertion(unsigned int lft = 0) : lft(lft)
    {

    }

    bool operator<(const CollapseInsertion & i) const {
        return lft < i.lft || (lft == i.lft && len() < i.len()) || (lft == i.lft && len() == i.len() && bases < i.bases);
    }

    bool operator==(const CollapseInsertion & i) const {
        return lft == i.lft && len() == i.len() && bases == i.bases;
    }

    unsigned int len() const {
        return bases.size();
    }

    unsigned int          lft;
    std::string           bases;
    std::string           quals;
};

using CollapseSplice = std::pair<unsigned int, unsigned int>;

struct ReadContig{
    std::vector<CollapseInsertion>   ins;

    size_t       index;
    size_t       island;
    unsigned int lft;
    unsigned int rgt;
    unsigned int qlft;
    unsigned int qrgt;

    unsigned int rlen() const {
        return rgt - lft + 1;
    }

    unsigned int qlen() const {
        return qrgt - qlft + 1;
    }

    CigarString   cig;
    ReadContig  * next = nullptr;
    ReadContig  * prev = nullptr;

    friend std::ostream & operator<<(std::ostream & o, const ReadContig & c){
        o   << "ref = " << c.lft << " - " << c.rgt << " query = " << c.qlft << " " << c.qrgt 
            << " cig = " << c.cig << " next = " << (c.next != nullptr) 
            << " prev = " << (c.prev != nullptr) << " idx = " << c.index;
        if(!c.ins.empty()){
            o << " Insertions =  ";
            for(size_t i = 0; i < c.ins.size(); i++){
                o << " " << c.ins[i].lft << " L = " << c.ins[i].len();
            }
        }
        return o; 
    }
};

struct SortBamDetailTidPos{
    uint32_t max_tid = 0;
    bool operator()(const BamDetail * lhs, const BamDetail * rhs) const {

        int32_t tid = lhs->b->core.tid == -1 ? max_tid : lhs->b->core.tid;
        uint64_t h1 = (static_cast<uint64_t>(tid) << 32) | (lhs->b->core.pos+1)<<1 | bam_is_rev(lhs->b);
        tid = rhs->b->core.tid == -1 ? max_tid : rhs->b->core.tid;
        uint64_t h2 = (static_cast<uint64_t>(tid) << 32) | (rhs->b->core.pos+1)<<1 | bam_is_rev(rhs->b);
        return h1 < h2;
    }
};

class ReadIsland{
    public:
        void merge(std::vector<BamDetail*> & umis, std::string & fbases, std::string & fquals, 
                   CigarString & fcigar, std::vector<CollapseSplice> & splices); 
                   //std::vector<unsigned int> & coverage);
        void debug();


        std::vector<ReadContig*>        contigs;
        //std::vector<CollapseSplice>     splices;
        std::vector<unsigned int>       starts;
        std::vector<unsigned int>       ends;
        std::vector<CollapseInsertion>  ins;
        std::vector<size_t>             conns;
        std::string                     bases;
        std::string                     quals;
        std::string                     cbases;
        std::string                     cquals;
        std::vector<unsigned int>       bcounts;
        std::string                     qmax;
        std::vector<uint8_t>            iquals;
        std::vector<unsigned int>       icounts;
        std::vector<unsigned int>       icoverage;
        size_t                          N;
        size_t                          L;
        unsigned int                    lft;
        unsigned int                    rgt;
        unsigned int                    ibases;
    private:
        void merge_bases_(std::vector<BamDetail *> & umis);
        void merge_ins_();
};



class CollapsedBamWriter{
    public:
        CollapsedBamWriter(const std::string & out, unsigned int bam_write_threads, const bam_hdr_t * bh);
        ~CollapsedBamWriter();

        void start();
        void join();
        void operator()();

        std::vector<BamDetail*>  collapsed;
        size_t                   wrote = 0;
        size_t                   cwrote = 0;

    private:
        std::thread              thread_;
        samFile                * bam_out_;
        bam_hdr_t              * bh_;
};

class BamOutputBuffer {
    public:
        ~BamOutputBuffer(){
            for(auto b : buffer_){
                delete b;
                //bam_destroy1(b);
            }
            buffer_.clear();
        }

        void get(std::vector<BamDetail*> & buffer, size_t N){
            std::lock_guard<std::mutex> lock(mutex_);
            while(N-- > 0){
                if(buffer_.empty()){
                    for(size_t i = 0; i < 100; i++){
                        allocated_++;
                        buffer.push_back(new BamDetail());
                    }
                }else{
                    buffer.push_back(buffer_.back());
                    buffer_.pop_back();
                }
            }
        }

        void ret(std::vector<BamDetail*> & buffer){
            std::lock_guard<std::mutex> lock(mutex_);
            while(!buffer.empty()){
                buffer_.push_back(buffer.back());
                buffer.pop_back();
            }
        }

    private:
        std::mutex                mutex_;
        std::vector<BamDetail*>   buffer_;
        size_t                    allocated_;
};

// Sort barcodes by their index
//From https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
inline void sort_reads_index(T start, T end, std::vector<size_t> & idx) {
    // initialize original index locations
    idx.clear();
    idx.resize(std::distance(start, end));
    std::iota(idx.begin(), idx.end(), 0);

    sort(idx.begin(), idx.end(),
        [&](size_t i1, size_t i2) {
            auto & d1 = *(*(start + i1));
            auto & d2 = *(*(start + i2));
            return std::tie(d1.gid, d1.hash, d1.pos) < std::tie(d2.gid, d2.hash, d2.pos);
        });

    //return idx;
}


}
