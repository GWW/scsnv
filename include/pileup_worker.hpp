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

#include "bam_genes.hpp"
#include "pileup.hpp"
#include "fasta.hpp"
#include "parallel_hashmap/parallel_hashmap/phmap.h"
#include "pileup_aux.hpp"
#include <thread>

namespace gwsc {

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
        PileupWorker(size_t NB, const TXIndex & txidx, const Fastas & genome, bool use_dups) 
            : 
            bcoverage(NB), bbases(NB), pup_(use_dups), slocs_(txidx), txidx_(txidx), 
            genome_(genome) 
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

        void reset(int32_t tid, uint32_t pos){
            pcount = 0;
            coverage.clear();
            targets_.rewind(tid, pos);
        }

        void set_debug(bool debug){
            debug_ = debug;
        }

        void process_range(BamBuffer::rpair range);

        void set_params(unsigned int min_alternative,
                unsigned int min_qual, unsigned int min_barcodes, unsigned int min_edge, unsigned int min_splice,
                double min_af, TargetFinder::trefs * refs = nullptr) {

            //min_coverage_ = min_coverage;
            min_alternative_ = min_alternative;
            min_qual_ = min_qual;
            slocs_.sdist = min_splice;
            min_barcodes_ = min_barcodes;
            min_edge_ = min_edge;
            min_af_ = min_af;
            if(refs != nullptr){
                targets_.set_refs(*refs);
                tloaded_ = true;
            }
        }

        std::vector<PositionCount>    positions;
        std::vector<PositionCoverage> coverage;
        std::vector<unsigned int>     bcoverage;
        std::vector<unsigned int>     bbases;
        //std::vector<BarcodeRate>      barcode_rates;
        //BarcodeRate                   total_rates;


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
        std::vector<uint16_t>                                         btypes_;
        TargetFinder                                                  targets_;
        SpliceLocations                                               slocs_;

        const TXIndex                                               & txidx_;
        const Fastas                                                & genome_;
        BamBuffer                                                   * buffer_;
        std::vector<PileupOut>                                        out_;
        std::vector<uint32_t>                                         gids_;
        std::string                                                   tmp_;
        phmap::flat_hash_map<uint32_t, std::array<uint16_t, 8>>       bcounter_;

        unsigned int                                                  min_alternative_ = 10;
        unsigned int                                                  min_qual_ = 20;
        unsigned int                                                  min_barcodes_ = 15;
        double                                                        min_af_ = 0.01;
        unsigned int                                                  min_edge_ = 5;
        bool                                                          tloaded_ = false;
        bool                                                          debug_ = false;
};

}
