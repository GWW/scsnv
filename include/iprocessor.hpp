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
#include "parallel-hashmap/parallel_hashmap/phmap.h"
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
    unsigned int index;
    char ref;
    char alt;
    //char strand;
    //bool found = false;
    uint32_t pref= 0;
    uint32_t palt= 0;
    uint32_t mref= 0;
    uint32_t malt= 0;
    uint32_t r_barcodes = 0;
    uint32_t a_barcodes = 0;
    uint32_t pr_barcodes = 0;
    uint32_t pa_barcodes = 0;
    uint32_t par_barcodes = 0;
    uint32_t mr_barcodes = 0;
    uint32_t ma_barcodes = 0;
    uint32_t mar_barcodes = 0;
    uint32_t tar_barcodes = 0;
    uint32_t total_barcodes = 0;
};

struct AccSNVGroup{
    AccSNVGroup(std::vector<AccSNV>::iterator start, std::vector<AccSNV>::iterator end, size_t end_pos)
        : start(start), end(end), end_pos(end_pos)
    {

    }
    std::vector<AccSNV>::iterator start;
    std::vector<AccSNV>::iterator end;
    size_t                        end_pos;
};

struct BarcodeCall{
    BarcodeCall(){
    }

    BarcodeCall(uint32_t barcode, uint32_t sindex, uint32_t pref, uint32_t mref, uint32_t palt, uint32_t malt)
        :barcode(barcode), sindex(sindex), pref(pref), mref(mref), palt(palt), malt(malt)
    {
    }

    bool operator<(const BarcodeCall & b) const {
        return std::tie(sindex, barcode) < std::tie(b.sindex, b.barcode);
    }

    uint32_t barcode;
    uint32_t sindex;
    uint32_t pref = 0;
    uint32_t mref = 0;
    uint32_t palt = 0;
    uint32_t malt = 0;
};

class IndexProcessor{
    public:
        using getrange = std::function<std::pair<size_t, size_t>(size_t)>;
        IndexProcessor(const std::string & name, ProcessorBase * processor, const TXIndex & index, Fastas & genome, bool use_dups, 
                std::string & bam_file, const phmap::flat_hash_map<std::string, unsigned int> & bchash, size_t BC, bool filter_barcodes) 
            : bchash_(bchash), name_(name), processor_(processor), pileup_(bchash.size(), index, genome, use_dups), filter_barcodes_(filter_barcodes)

        {
            //code 0  for genome 1 for scsnv 2 for cellranger
            bf_ = sam_open(bam_file.c_str(), "r");
            bh_ = sam_hdr_read(bf_);
            bi_ = sam_index_load(bf_, bam_file.c_str());
            barcodes_ = BC;
            targets_.resize(genome.size());

            if(bi_ == NULL){
                std::cerr << "Could not open bam index for " << bam_file << "\n";
                exit(1);
            }
        }

        void make_targets(const std::vector<AccSNV> & snvs){
            total_targets_ = 0;
            for(auto & s : snvs){
                int tid = bam_name2id(bh_, s.chrom.c_str());
                targets_[tid].push_back(TargetFinder::Target(s.pos, s.alt, s.index));
                total_targets_ = std::max(s.index, total_targets_);
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

        const std::string & name() const {
            return name_;
        }

        PileupWorker & worker() {
            return pileup_;
        }

        size_t N() const {
            return curr_;
        }

        TargetFinder::trefs & targets() {
            return targets_;
        }

        void start(getrange rng, std::vector<AccSNVGroup> & snvs){
            getrange_ = rng;
            snvs_ = &snvs;
            thread_ = std::thread(std::ref(*this));
        }

        bool fbarcodes() const {
            return filter_barcodes_;
        }

        void join(){
            thread_.join();
        }

        void operator()();

        std::vector<BarcodeCall>     calls;

    private:
        std::thread                  thread_;
        const phmap::flat_hash_map<std::string, unsigned int> & bchash_;
        getrange                     getrange_;
        std::string                  name_;
        ProcessorBase              * processor_ = nullptr;
        std::vector<AccSNVGroup>   * snvs_;
        PileupWorker                 pileup_;
        std::vector<BamDetail*>      recs_;
        TargetFinder::trefs          targets_;
        std::string                  btmp_;

        hts_idx_t                  * bi_ = nullptr;
        bam_hdr_t                  * bh_ = nullptr;
        samFile                    * bf_ = nullptr;
        unsigned int                 curr_ = 0;
        unsigned int                 read_ = 0;
        unsigned int                 total_targets_ = 0;
        unsigned int                 barcodes_ = 0;
        int                          tid_ = -1;
        bool                         filter_barcodes_ = false;
};


struct DataManager{
    std::vector<AccSNV>           snvs;
    std::vector<AccSNVGroup>      groups;
    std::vector<IndexProcessor*>  processors;
    std::string                   name;
    size_t                        curr_ = 0;
    size_t                        last_ = 0;
    std::mutex                    mtx_;

    std::pair<size_t, size_t> get_range(size_t N) {
        std::lock_guard<std::mutex> lock(mtx_);
        size_t start = curr_;
        curr_ = std::min(curr_ + N, groups.size());
        if((curr_ - last_) >= 10000){
            tout << "Processing " << name << " " << curr_ << " group out of " << groups.size() << "\n";
            last_ = curr_;
        }
        return {start, curr_};
    }

    ~DataManager(){
        for(auto & p : processors) delete p;
        processors.clear();
    }

    void run();
    void write(const std::string & out, const std::string & header, const std::vector<std::string> & barcode_ids);
    std::vector<BarcodeCall>     calls;
};

}
