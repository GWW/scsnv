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

#include "bam_genes.hpp"
#include "../scsnv/aux.hpp"
#include "../scsnv/align_aux.hpp"

#include <thread>
#include <list>
#include <array>

namespace gwsc {

struct PileupRead {
    void init(BamDetail & d);

    bool next();

    bool between() const {
        return it == cigar.end();
    }

    uint8_t base() const {
        //     "=ACMGRSVTWYHKDBN";
        return "4014244434444444"[bam_seqi(bam_get_seq(d->b), qpos)] - '0';
    }

    uint8_t qual() const {
        return *(bam_get_qual(d->b) + qpos);
    }

    int tid() const {
        return d->b->core.tid;
    }

    int qlen() const {
        return d->b->core.l_qseq;
    }

    BamDetail                   * d;
    CigarString                   cigar;
    CigarString::const_iterator   it;
    unsigned int                  cpos;
    unsigned int                  qpos;
    unsigned int                  rpos;
    unsigned int                  rend;
    unsigned int                  ibases;
    unsigned int                  dbases;
    unsigned int                  abases;
    unsigned int                  NM;
    char                          strand;
};

struct PileupOut{
    PileupOut(){

    }

    PileupOut(PileupRead & r) 
        : 
            d(r.d), rpos(r.rpos), qpos(r.qpos), 
            base(r.base()), qual(r.qual()), ibases(r.ibases), 
            dbases(r.dbases), abases(r.abases), NM(r.NM), rev(r.strand == '-')
    {

    }

    BamDetail *       d;
    unsigned int      rpos;
    unsigned int      qpos;
    uint8_t           base;
    uint8_t           qual;
    unsigned int      ibases;
    unsigned int      dbases;
    unsigned int      abases;
    unsigned int      NM;
    bool              rev;
};

//template <typename T>
class Pileup {
    public:
        using const_iterator = std::vector<PileupOut>::const_iterator;

        Pileup(bool use_dups = false) 
            : dups_(use_dups) {

        }

        ~Pileup(){
            for(auto & r : reads_){
                delete r;
            }
            reads_.clear();
        }

        void build_reads(BamBuffer::rpair range);
        bool next(std::vector<PileupOut> & out);

        int                tid = -1;
        unsigned int       pos = 0;

    private:
        void get_reads_();
        void pileup_reads_(std::vector<PileupOut> & out);

        std::list<PileupRead*>       buffer_;
        std::vector<PileupRead*>     reads_;
        size_t                       idx_ = 0;
        size_t                       count_ = 0;
        int                          ctid_ = -1;
        unsigned int                 cpos_ = 0;
        bool                         dups_;
};

}
