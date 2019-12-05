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


#pragma once

#include <string>
#include <cstring>
#include <iostream>
#include <vector>
#include "../util/misc.hpp"
#include "../scmap/index.hpp"
#include "sbam_merge.hpp"
#include <htslib/hts.h>
#include <htslib/sam.h>
#include "../util/sequence.hpp"
#include <mutex>

namespace gwsc{

struct BamDetail {
    typedef uint64_t hash_type;
    BamDetail(){
        b = bam_init1();
    }

    ~BamDetail(){
        if(b != nullptr) bam_destroy1(b);
        b = nullptr;
    }

    void copy(BamDetail & src){
        //bam_copy1(b, src.b);
        std::swap(b, src.b);
        umi = src.umi;
        gid = src.gid;
        rnum = src.rnum;
        barcode = src.barcode;
        dup = src.dup;
        strand = src.strand;
    }

    void make_hash(){
        hash = barcode;
        hash <<= 32;
        hash |= umi;
    }

    hash_type hash;
    bam1_t * b = nullptr;
    AlignSummary::bint barcode = 0;
    uint32_t           umi = 0;
    uint32_t           gid = 0;
    uint32_t           rnum = 0;
    int                pos = 0;
    int                lft = 0;
    int                rgt = 0;
    char               xt = '?';
    char               strand = '?';
    bool               corrected = false;
    bool               dup = false;
};

class BamBuffer {
    template <typename T>
    friend class BamGeneReader;

    public: 
        using read_vec = std::vector<BamDetail*>;
        using rit = read_vec::iterator;
        using rpair = std::pair<rit, rit>;

        void reset(){
            count_ = 0;
            idx_ = 0;
            current_ = 0;
            starts_.clear();
        }

        rit begin() {
            return reads_.begin();
        }

        rit end() {
            return reads_.end();
        }
    
        void add(BamDetail & d){
            assure_size();
            reads_[count_++]->copy(d);
        }

        size_t count() const {
            return count_;
        }

        void assure_size(){
            if((count_ + 1) >= reads_.size()){
                for(size_t i = 0; i < 20; i++){
                    reads_.push_back(new BamDetail());
                }
            }
        }

        rpair get_all_unsafe(){
            return {reads_.begin(), reads_.begin() + starts_.back()};
        }

        bool get_next(rpair & ret){
            std::lock_guard<std::mutex> lock(mutex_);
            if(current_ < (starts_.size() - 1)){
                rit start  = reads_.begin() + starts_[current_++];
                rit end = reads_.begin() + starts_[current_];
                ret = {start, end};
                return true;
            }
            ret = {reads_.end(), reads_.end()};
            return false;
        }

        ~BamBuffer(){
            //std::cout << "Destroying buffer with " << reads_.size() << " reads\n";
            for(size_t i = 0; i < reads_.size(); i++){
                delete reads_[i];
            }
            reads_.clear();
            count_ = 0;
        }

    private:
        read_vec                  reads_;
        std::vector<unsigned int> starts_;
        std::mutex                mutex_;
        unsigned int              count_ = 0;
        unsigned int              idx_ = 0;
        unsigned int              current_ = 0;
};

template <typename T>
class BamGeneReader{
    public:

        BamGeneReader(bool cellranger) : cellranger_(cellranger){
        }

        ~BamGeneReader(){
        }

        template <typename IT>
        void add_bams(IT start, IT end) {
            bm_.add_bams(start, end);
        }

        void prepare(){
            if(!get_()) {
                std::cerr << "bam shouldn't be empty\n";
                exit(1);
            }
            ltid_ = index.gene(next_.gid).tid;
            max_rgt_ = index.gene(next_.gid).rgt;
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

        TXIndex                    index;

    private:
        unsigned int read_(BamBuffer & buffer);
        bool         get_();

        TXIndex::Ref::itree::intervalVector overlaps_;
        BamMerger                 bm_;
        BamDetail                 next_;
        std::mutex                mutex_;
        unsigned int              total_ = 0;
        unsigned int              idx_ = 0;
        unsigned int              umi_len_ = T::UMI_LEN;
        int                       max_rgt_ = 0;
        int                       ltid_ = -1;
        bool                      done_ = false;
        bool                      cellranger_ = false;
};

template <typename T>
inline unsigned int BamGeneReader<T>::read_genes(BamBuffer & buffer, size_t N){
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

template <typename T>
inline bool BamGeneReader<T>::get_(){
    if(cellranger_){
        while(!done_ && bm_.next(next_.b) != nullptr){
            if(next_.b->core.qual != 255) {
                continue;
            }

            if(bam_aux_get(next_.b, "CB") == NULL || bam_aux_get(next_.b, "UB") == NULL){
                continue;
            }

            char A = bam_aux2A(bam_aux_get(next_.b, "RE"));
            auto ptr = bam_aux_get(next_.b, "MM");
            if(A == 'E' && ptr != NULL){
                unsigned int MM = bam_aux2i(ptr);
                if(MM == 1) {
                    continue;
                }
            }
            uint32_t gid = std::numeric_limits<uint32_t>::max();
            bool rev = bam_is_rev(next_.b);
            char xs = ((T::LibraryStrand == StrandMode::TAG_FWD && !rev) || (T::LibraryStrand == StrandMode::TAG_REV && rev)) ? '+' : '-';
            if(A == 'N'){
                auto const & ref = index.ref(next_.b->core.tid);
                if(!ref.empty()){
                    overlaps_.clear();
                    ref.tree.findOverlapping(next_.b->core.pos, bam_endpos(next_.b) - 1, overlaps_);
                    unsigned int count = 0;
                    for(size_t i = 0; i < overlaps_.size(); i++){
                        auto tgid = overlaps_[i].value;
                        const GeneEntry & gene = index.gene(tgid);
                        if(gene.strand == xs){
                            gid = tgid;
                            count++;
                        }
                    }
                    if(count == 1){
                        const GeneEntry & gene = index.gene(gid);
                        const char * gstr = gene.gene_id.c_str();
                        bam_aux_append(next_.b, "XG", 'Z', strlen(gstr) + 1, reinterpret_cast<const uint8_t *>(gstr));
                    }else{
                        gid = std::numeric_limits<uint32_t>::max();
                    }
                }
            }else if(A == 'E'){
                ptr = bam_aux_get(next_.b, "GX");
                char * XG = ptr == NULL ? NULL : bam_aux2Z(ptr);
                if(XG != NULL && strchr(XG, ';') == NULL){
                    gid = index.gid_from_str(XG);
                    bam_aux_append(next_.b, "XG", 'Z', strlen(XG) + 1, reinterpret_cast<const uint8_t *>(XG));
                }
            }else{
                continue;
            }

            if(gid != std::numeric_limits<uint32_t>::max()){
                bam_aux_append(next_.b, "XS", 'A', 1, reinterpret_cast<const uint8_t*>(&xs));
                next_.rnum = total_;
                next_.b->core.flag &= ~BAM_FDUP;
                next_.gid = gid;
                total_++;
                return true;
            }
        }
        done_ = true;
        return false;
    }
    if(!done_ && bm_.next(next_.b) != nullptr){
        next_.gid = index.gene_from_id(bam_aux2Z(bam_aux_get(next_.b, "XG"))).gid;
        next_.rnum = total_;
        total_++;
        return true;
    }
    done_ = true;
    return false;
}

template <typename T>
inline unsigned int BamGeneReader<T>::read_(BamBuffer & buffer){
    size_t r = 0;
    while(next_.b->core.tid == ltid_ && next_.b->core.pos <= max_rgt_){
        int32_t gid = next_.gid;
        max_rgt_ = std::max(static_cast<int>(index.gene(gid).rgt), max_rgt_);
        buffer.assure_size();
        buffer.add(next_);
        r++;
        if(!get_()) return r;
    }
    idx_++;
    int32_t gid = next_.gid;;
    ltid_ = index.gene(gid).tid;
    max_rgt_ = index.gene(gid).rgt;
    return r;
}

}
