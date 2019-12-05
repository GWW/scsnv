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

#include "../muts/bam_genes.hpp"
#include "pileup.hpp"
#include "../util/fasta.hpp"
#include "../sparsepp/sparsepp/spp.h"
#include <thread>

namespace gwsc {

struct PositionCoverage {
    PositionCoverage(int tid = -1, unsigned int pos = 0, unsigned int pcoverage = 0, unsigned int mcoverage = 0) 
        :tid(tid), pos(pos), pcoverage(pcoverage), mcoverage(mcoverage)
    {
    }

    bool operator<(const PositionCoverage & p) const {
        return tid < p.tid || (tid == p.tid && pos < p.pos);
    }
    int32_t tid = -1;
    uint32_t pos = 0;
    uint32_t pcoverage = 0;
    uint32_t mcoverage = 0;
};

struct BarcodeCount{
    BarcodeCount(uint32_t barcode)
        : barcode(barcode)
    {
        for(size_t i = 0; i < 4; i++) pbases[i] = mbases[i] = 0;
    }

    void inc_base(uint8_t b, uint32_t count, bool rev) {
        if(rev){
            mbases[b] = std::min((uint32_t)std::numeric_limits<uint16_t>::max(), (uint32_t)mbases[b] + count);
        }else{
            pbases[b] = std::min((uint32_t)std::numeric_limits<uint16_t>::max(), (uint32_t)pbases[b] + count);
        }
    }

    uint32_t barcode;
    uint16_t pbases[4];
    uint16_t mbases[4];
};

struct BarcodeRate {
    static const size_t BTOTAL = 24;
    BarcodeRate()
    {
        for(size_t i = 0; i < BTOTAL; i++) prates[i] = mrates[i] = 0;
    }

    BarcodeRate & operator+=(const BarcodeRate & rhs){
        for(size_t i = 0; i < BTOTAL; i++) {
            prates[i] += rhs.prates[i];
            mrates[i] += rhs.mrates[i];
        }
        return *this;
    }
    uint32_t prates[BTOTAL];
    uint32_t mrates[BTOTAL];
};

struct PositionCount{
    void reset() {
        bases = {};
        tid = -1;
        pos = 0;
        coverage = 0;
        barcodes = 0;
        ambig = 0;
        ref = 'N';
        refi = 4;
        pbase = '-';
        mbase = '-';
        bcounts.clear();
    }

    bool operator<(const PositionCount & p) const {
        return tid < p.tid || (tid == p.tid && pos < p.pos);
    }

    struct BaseCount{
        unsigned int         p_count = 0;
        //unsigned int         p_NM = 0;
        //unsigned int         p_gaps = 0;
        //unsigned int         p_lens = 0;
        unsigned int         m_count = 0;
        //unsigned int         m_NM = 0;
        //unsigned int         m_gaps = 0;
        //unsigned int         m_lens = 0;
        unsigned int         t_barcodes = 0;
    };

    std::array<BaseCount, 4>  bases;
    std::vector<BarcodeCount> bcounts;
    int                       tid = -1;
    unsigned int              pos = 0;
    unsigned int              coverage = 0;
    unsigned int              barcodes = 0;
    unsigned int              ambig = 0;
    char                      ref = 'N';
    char                      pbase = '-';
    char                      mbase = '-';
    uint8_t                   refi = 4;
};

class PileupWorker {
    PileupWorker( const PileupWorker& ) = delete;
    PileupWorker& operator=(const PileupWorker&) = delete;

    enum SBaseType {
        PINTRON=1,
        PEXON=2,
        PNC=4,
        P5UTR=8,
        P3UTR=16,
        MINTRON=32,
        MEXON=64,
        MNC=128,
        M5UTR=256,
        M3UTR=512
    };

    public:
        using cb_bhash = std::function<AlignSummary::bint(std::string &)>;
        PileupWorker(const TXIndex & txidx, const Fastas & genome, cb_bhash cb, bool use_dups, bool cellranger) 
            : 
            pup_(use_dups), bhash_(cb), txidx_(txidx), genome_(genome),
            cellranger_(cellranger) 
        {
        }

        ~PileupWorker(){

        }

        void start(){
            thread_ = std::thread(std::ref(*this));
        }

        void set_buffer(BamBuffer * buffer){
            buffer_ = buffer;
        }

        void join() {
            thread_.join();
        }

        void operator()();

        void set_params(unsigned int min_alternative,
                unsigned int min_qual, unsigned int min_barcodes, double min_af) {

            //min_coverage_ = min_coverage;
            min_alternative_ = min_alternative;
            min_qual_ = min_qual;
            min_barcodes_ = min_barcodes;
            min_af_ = min_af;
        }


        std::vector<PositionCount>    positions;
        std::vector<PositionCoverage> coverage;
        std::vector<BarcodeRate>      barcode_rates;
        BarcodeRate                   total_rates;


        size_t pcount = 0;
        size_t plus_bases = 0;
        size_t minus_bases = 0;
        size_t bases = 0;
        size_t pbases = 0;
        size_t reads = 0;

    private:

        void count_barcodes_(PositionCount & p);

        unsigned int build_genes_();

        std::thread                                                   thread_;
        std::vector<size_t>                                           sidx_;
        Pileup                                                        pup_;
        cb_bhash                                                      bhash_;
        std::vector<uint16_t>                                         btypes_;

        const TXIndex                                               & txidx_;
        const Fastas                                                & genome_;
        BamBuffer                                                   * buffer_;
        std::vector<PileupOut>                                        out_;
        std::vector<uint32_t>                                         gids_;
        std::string                                                   tmp_;
        spp::sparse_hash_map<uint32_t, std::array<uint16_t, 8>>       bcounter_;

        size_t                                                        min_alternative_ = 10;
        size_t                                                        min_qual_ = 20;
        size_t                                                        min_barcodes_ = 15;
        double                                                        min_af_ = 0.01;
        bool                                                          cellranger_;
};

}
