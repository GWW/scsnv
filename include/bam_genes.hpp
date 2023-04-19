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


#pragma once

#include <string>
#include <cstring>
#include <iostream>
#include <functional>
#include <vector>
#include "misc.hpp"
#include "index.hpp"
#include "sbam_merge.hpp"
#include "sequence.hpp"
#include <mutex>
#include <random>
#include "bam_genes_aux.hpp"
namespace gwsc{

template <typename T, typename R, typename P>
class BamGeneReader{
    public:

        BamGeneReader() : rf_(index, T::LibraryStrand){
        }

        virtual ~BamGeneReader(){

        }

        template <typename IT>
        void add_bams(IT start, IT end) {
            bm_.add_bams(start, end);
        }

        void set_bam(const std::string & bf){
            bm_.set_bam(bf);
        }

        void set_threads(unsigned int threads){
            bm_.set_threads(threads);
        }

        void prepare(){
            if(!get_()) {
                std::cerr << "bam shouldn't be empty\n";
                exit(1);
            }
            ltid_ = next_.b->core.tid;
            max_rgt_ = std::max(static_cast<int>(bam_endpos(next_.b)), static_cast<int>(index.gene(next_.gid).rgt));
        }

        const bam_hdr_t * header() const {
            return bm_.header();
        }

        size_t total() {
            return total_;
        }

        size_t mtotal() {
            return bm_.total();
        }

        unsigned int read_genes(BamBuffer & buffer, size_t N);

        unsigned int read_genes_safe(BamBuffer & buffer, size_t N){
            std::lock_guard<std::mutex> lock(mutex_);
            return read_genes(buffer, N);
        }

        const std::vector<unsigned int> & file_counts() const {
            return bm_.file_totals();
        }

        void set_max_reads(uint64_t max_reads){
            max_reads_ = max_reads;
        }

        TXIndex                    index;

    protected:
        virtual bool get_();
        unsigned int read_(BamBuffer & buffer);

        uint64_t                  max_reads_ = 5000000;
        P                         rf_;
        R                         bm_;
        BamDetail                 next_;
        std::mutex                mutex_;
        unsigned int              total_ = 0;
        unsigned int              idx_ = 0;
        int                       max_rgt_ = 0;
        int                       ltid_ = -1;
        bool                      done_ = false;
};

template <typename T, typename R, typename P>
inline unsigned int BamGeneReader<T, R, P>::read_genes(BamBuffer & buffer, size_t N){
    buffer.reset();
    unsigned int read = 0;
    {
        buffer.idx_ = idx_;
        for(size_t i = 0; i < N; i++){
            if(done_) {
                buffer.starts_.push_back(buffer.count_);
                return read;
            }
            buffer.starts_.push_back(buffer.count_);
            read += read_(buffer);
        }
        buffer.starts_.push_back(buffer.count_);
    }
    return read;
}

template <typename T, typename R, typename P>
inline bool BamGeneReader<T,R,P>::get_(){
    unsigned int fno = 0;
    while(!done_ && bm_.next(next_.b, &fno) != nullptr){
        if(rf_(next_, fno)) {
            total_++;
            return true;
        }
    }
    done_ = true;
    return false;
}

template <typename T, typename R, typename P>
inline unsigned int BamGeneReader<T, R, P>::read_(BamBuffer & buffer){
    size_t r = 0;
    int rstart = 0;
    size_t start_count = buffer.count();
    bool warned = false;
    //std::cout << "Reading gene group max_rgt = " << max_rgt_ << "\n";
    //int lpos = next_.b->core.pos;
    while(next_.b->core.tid == ltid_ && next_.b->core.pos < max_rgt_){
        int32_t gid = next_.gid;
        rstart = next_.b->core.pos;
        max_rgt_ = std::max(static_cast<int>(index.gene(gid).rgt), max_rgt_);
        max_rgt_ = std::max(static_cast<int>(bam_endpos(next_.b)), max_rgt_);
        if(r < max_reads_){
            warned = true;
            buffer.add(next_);
        }
        r++;
        if(!get_()) return r;
    }

    if(warned){
        std::cerr << "[Warning] exceeded " << max_reads_ << " skipping this gene region " << r << " reads. Region " 
            << index.tname(ltid_) << ":" << rstart << " - " << max_rgt_ << "\n";
        buffer.revert(start_count);
    }
    //std::cout << "  Done " << r << " total\n";
    idx_++;
    int32_t gid = next_.gid;
    ltid_ = index.gene(gid).tid;
    max_rgt_ = std::max(static_cast<int>(bam_endpos(next_.b)), static_cast<int>(index.gene(next_.gid).rgt));
    return r;
}

template <typename T, typename R, typename P>
class BamGeneReaderFiltered: public BamGeneReader<T, R, P> {
    public:
        using bcfilter_func = std::function<unsigned int(const std::string & bc, unsigned int fno)>;
        BamGeneReaderFiltered(bcfilter_func & filter, std::vector<double> downsamples = std::vector<double>()) 

            : BamGeneReader<T,R,P>(), ds_(downsamples), filter_(filter), dist_(0, 1)
        {
            std::random_device rd;
            seed_ = rd();
            gen_ = std::mt19937(seed_);
            ds_kept.resize(ds_.size());
            ds_skipped.resize(ds_.size());
        }

        void add_rnd(){
            add_rand_ = true;
        }

        void set_seed(unsigned int seed){
            gen_ = std::mt19937(seed);
            seed_ = seed;
        }

        unsigned int get_seed() const {
            return seed_;
        }

        virtual ~BamGeneReaderFiltered() {

        }

        std::vector<unsigned int>               ds_skipped;
        std::vector<unsigned int>               ds_kept;

    protected:
        std::string btmp_;

        virtual bool get_(){
            unsigned int fno = 0;

            while(!BamGeneReader<T,R,P>::done_ && BamGeneReader<T,R,P>::bm_.next(BamGeneReader<T,R,P>::next_.b, &fno) != nullptr){
                if(!BamGeneReader<T,R,P>::rf_(BamGeneReader<T,R,P>::next_, fno)){
                    continue;
                }
                btmp_ = bam_aux2Z(bam_aux_get(BamGeneReader<T,R,P>::next_.b, "CB"));
                size_t dash = btmp_.find('-');
                if(dash != std::string::npos){
                    btmp_ = btmp_.erase(dash);
                }
                unsigned int bid = filter_(btmp_, fno);
                if(bid == std::numeric_limits<unsigned int>::max()) continue;
                if(!ds_.empty()){
                    if(dist_(gen_) >= ds_[fno]){
                        ds_skipped[fno]++;
                        continue;
                    }
                    ds_kept[fno]++;
                }else if(add_rand_){
                    BamGeneReader<T,R,P>::next_.rnd = dist_(gen_);
                }
                BamGeneReader<T,R,P>::next_.barcode = bid;
                BamGeneReader<T,R,P>::next_.filenum = fno;
                BamGeneReader<T,R,P>::total_++;
                return true;
            }

            BamGeneReader<T,R,P>::done_ = true;
            return false;
        }

    private:
        std::vector<double>                     ds_;
        bcfilter_func                           filter_;
        std::mt19937                            gen_;
        std::uniform_real_distribution<double>  dist_;
        unsigned int                            seed_;
        bool                                    add_rand_;
};

}
