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

#include <string>
#include <vector>
#include <cstdint>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <htslib/hts.h>
#include <htslib/thread_pool.h>
#include <htslib/sam.h>
#include <iostream>
#include <sstream>
#include <thread>
#include "../scmap/align_aux.hpp"
#include "../scmap/reader.hpp"

namespace gwsc{

struct SortBamTidPos{
    uint32_t max_tid = 0;
    bool operator()(const bam1_t * lhs, const bam1_t * rhs) const {

        int32_t tid = lhs->core.tid == -1 ? max_tid : lhs->core.tid;
        uint64_t h1 = (static_cast<uint64_t>(tid) << 32) | (lhs->core.pos+1)<<1 | bam_is_rev(lhs);
        tid = rhs->core.tid == -1 ? max_tid : rhs->core.tid;
        uint64_t h2 = (static_cast<uint64_t>(tid) << 32) | (rhs->core.pos+1)<<1 | bam_is_rev(rhs);
        return h1 < h2;
    }
};

class BamHeader{
    BamHeader( const BamHeader& ) = delete;
    BamHeader& operator=(const BamHeader&) = delete;

    public:
        BamHeader() {

        }

        struct RefData{
            RefData() {

            }

            RefData(const std::string & name, uint32_t len, int32_t tid) :
                name(name), tid(tid), len(len) {

            }

            std::string name;
            int32_t     tid;
            uint32_t    len;

            bool operator <(const RefData & rhs) const {
                return tid < rhs.tid;
            }
        };

        ~BamHeader() {
            if(hdr_ != NULL) {
                bam_hdr_destroy(hdr_);
                hdr_ = NULL;
            }
        }

        void add_name(const std::string & name, uint32_t len, int32_t tid){
            refs_.push_back(RefData(name, len, tid));
        }

        const bam_hdr_t * bam_hdr() const {
            return hdr_;
        }

        bam_hdr_t * bam_hdr(){
            return hdr_;
        }

        void build_header(){
            if(hdr_ != NULL){
                std::cout << "Header already built\n";
                exit(1);
            }
            std::sort(refs_.begin(), refs_.end());
            hdr_ = bam_hdr_init();
            hdr_->n_targets     = refs_.size();
            hdr_->target_name   = (char**)malloc(hdr_->n_targets * sizeof(char*));
            hdr_->target_len    = (uint32_t*)malloc(hdr_->n_targets * sizeof(uint32_t));

            for(size_t i = 0; i < refs_.size(); i++){
                auto const & f = refs_[i];
                hdr_->target_len[i]  = f.len;
                hdr_->target_name[i] = (char*)malloc((f.name.size() + 1) * sizeof(char));
                strcpy(hdr_->target_name[i], f.name.c_str());
            }
        }

        void set_text(std::string text){
            if(hdr_ == NULL){
                std::cout << "Must call build_header before adding text\n";
                exit(1);
            }
            hdr_->text = (char*)malloc((text.size() + 1) * sizeof(char));
            hdr_->l_text = text.size() + 1;
            strcpy(hdr_->text, text.c_str());
        }

    private:
        std::vector<RefData>   refs_;
        bam_hdr_t            * hdr_ = NULL;

};

class SortedBamWriter{
    SortedBamWriter( const SortedBamWriter& ) = delete;
    SortedBamWriter& operator=(const SortedBamWriter&) = delete;
    public:
        using read_buffer = std::vector<bam1_t*>;

        SortedBamWriter() : bh(), pool_{NULL, 0} {

        }

        ~SortedBamWriter(){
            for(auto p : out_){
                bam_destroy1(p);
            }
            out_.clear();
            if(pool_.pool != NULL){
                hts_tpool_destroy(pool_.pool);
            }
        }

        void make_pool(unsigned int threads){
            if(threads > 1){
                pool_.pool = hts_tpool_init(threads);
            }
            pool_threads_ = threads;
        }

        void set_thread_buffer(unsigned int out_sz, unsigned int thread_sz){
            out_sz_ = out_sz;
            thread_sz_ = thread_sz;
            
            for(size_t i = 0; i < out_sz_; i++){
                out_.push_back(bam_init1());
            }
        }

        unsigned int total_reads() const {
            return total_reads_;
        }

        unsigned int file_reads() const {
            return file_reads_;
        }

        unsigned int total_files() const {
            return file_number_;
        }

        unsigned int thread_size() const {
            return thread_sz_;
        }

        void prepare_thread_buffer(read_buffer & buff) {
            for(size_t i = 0; i < thread_sz_; i++){
                buff.push_back(bam_init1());
            }
        }

        void set_prefix(const std::string & prefix){
            prefix_ = prefix;
        }

        // Add a read to the output buffer, write the file if necessary **THREAD SAFE**
        void write(const AlignGroup & g, const Read & r, read_buffer & buffer, unsigned int & rcount, const std::string & gid);
        // To write any left over reads directly **NOT THREAD SAFE**
        void merge_buffer(read_buffer & buffer, unsigned int & rcount);
        // Write everything in the buffer **NOT THREAD SAFE** 
        void force_write();
        
        BamHeader bh;

    private:
        void _align2bam(bam1_t * bam, const AlignGroup & g, const Read & r, const std::string & gid);
        std::mutex            mtx_write_;
        unsigned int          thread_sz_ = 10000;
        unsigned int          out_sz_    = 2000000;
        read_buffer           out_;
        unsigned int          file_number_ = 0;
        unsigned int          total_reads_ = 0;
        unsigned int          file_reads_ = 0;
        unsigned int          pool_threads_ = 1;
        std::string           prefix_;
        htsThreadPool         pool_;
};

}
