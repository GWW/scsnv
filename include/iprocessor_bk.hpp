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

#include "gzstream.hpp"
#include "htslib/htslib/hts_endian.h"
#include "htslib/htslib/sam.h"
#include "index.hpp"
#include "reader.hpp"
#include "bam_genes_aux.hpp"
#include "fasta.hpp"
#include "pileup_worker.hpp"
#include <thread>
#include <mutex>
#include <condition_variable>

namespace gwsc{

class BamGenomeProcessor : public ProcessorBase {
    public: 
        BamGenomeProcessor(const TXIndex & index, StrandMode ST){
            index_ = &index;
            ST_ = ST;
        }

        bool operator()(BamDetail & read, unsigned int fno){
            if(read.b->core.qual >= qual_) {
                read.filenum = fno;
                return true;
            }
            return false;
        }

        void set_mapq(int mapq){
            qual_ = mapq;
        }

    private:
        unsigned int    qual_ = 30;
};


struct AccSNV{
    std::string chrom;
    std::string line;
    unsigned int pos;
    uint32_t     founds = 0;
    char ref;
    char alt;
    char strand;
};

class IndexProcessor{
    public:
        IndexProcessor(const std::string & name, ProcessorBase * processor, const TXIndex & index, Fastas & genome, bool use_dups, 
                std::string & bam_file, const spp::sparse_hash_map<std::string, unsigned int> & bchash, bool filter_barcodes) 
            : bchash_(bchash), name_(name), processor_(processor), pileup_(bchash.size(), index, genome, use_dups), filter_barcodes_(filter_barcodes)

        {
            //code 0  for genome 1 for scsnv 2 for cellranger
            bf_ = sam_open(bam_file.c_str(), "r");
            bh_ = sam_hdr_read(bf_);
            bi_ = sam_index_load(bf_, bam_file.c_str());
            targets_.resize(genome.size());

            if(bi_ == NULL){
                std::cerr << "Could not open bam index for " << bam_file << "\n";
                exit(1);
            }
        }

        void make_targets(const std::vector<AccSNV> & snvs){
            for(auto & s : snvs){
                int tid = bam_name2id(bh_, s.chrom.c_str());
                targets_[tid].push_back(TargetFinder::Target(s.pos, s.alt));
            }
        }

       ~IndexProcessor(){
            if(processor_ != nullptr) delete processor_;
            if(bi_ != nullptr){
                hts_idx_destroy(bi_);
            }

            if(bh_ != nullptr){
                bam_hdr_destroy(bh_);
            }

            if(bf_ != nullptr){
                sam_close(bf_);
            }
            for(auto * d : recs_) delete d;
            recs_.clear();
            if (thread_.joinable()) thread_.join();
        }

        void run_position(const std::string & chrom, unsigned int pos);
        //void pileup_position();
        //void pileup_position_debug(const std::string & chrom, unsigned int pos);

        const std::string & name() const {
            return name_;
        }

        bool valid() const {
            return valid_;
        }

        PileupWorker & worker() {
            return pileup_;
        }

        void start(){
            thread_ = std::thread(std::ref(*this));
        }

        void wait();

        void done();

        size_t N() const {
            return curr_;
        }

        TargetFinder::trefs & targets() {
            return targets_;
        }

        std::vector<PositionCount>::const_iterator entry() const {
            return entry_;
        }

        void header(gzofstream & zout) const;
        void write(gzofstream & zout, unsigned int refi, unsigned int alti) const;
        void operator()();

    private:
        std::string             chrom_;
        std::mutex              mtx_wait_;
        std::condition_variable cv_;
        std::thread             thread_;
        const spp::sparse_hash_map<std::string, unsigned int> & bchash_;
        std::string             name_;
        ProcessorBase         * processor_ = nullptr;
        PileupWorker            pileup_;
        std::vector<PositionCount>::const_iterator entry_;
        
        std::vector<BamDetail*> recs_;
        TargetFinder::trefs     targets_;
        std::string             btmp_;
        hts_idx_t             * bi_ = nullptr;
        bam_hdr_t             * bh_ = nullptr;
        samFile               * bf_ = nullptr;
        unsigned int            pos_ = 0;
        unsigned int            curr_ = 0;
        unsigned int            read_ = 0;
        int                     tid_ = -1;
        bool                    filter_barcodes_ = false;
        bool                    valid_ = false;
        bool                    ready_ = false;
        bool                    done_ = false;
        bool                    processed_ = false;
};

}
