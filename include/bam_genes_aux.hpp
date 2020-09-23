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
#include "htslib/htslib/hts.h"
#include "htslib/htslib/sam.h"
#include "align_aux.hpp"
#include "index.hpp"

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
        std::swap(b, src.b);
        umi = src.umi;
        gid = src.gid;
        barcode = src.barcode;
        dup = src.dup;
        strand = src.strand;
        filenum = src.filenum;
        rnd = src.rnd;
    }

    void make_hash(){
        hash = barcode;
        hash <<= 32;
        hash |= umi;
    }

    hash_type hash;
    double rnd = 0.0;
    bam1_t * b = nullptr;
    AlignSummary::bint barcode = 0;
    unsigned int       filenum = 0;
    uint32_t           umi = 0;
    uint32_t           gid = 0;
    int                pos = 0;
    int                lft = 0;
    int                rgt = 0;
    char               xt = '?';
    char               strand = '?';
    bool               corrected = false;
    bool               dup = false;
    bool               processed = false;
};

inline std::string debug_read(const BamDetail & d) {
    std::stringstream ss;
    const bam1_t * b = d.b;
    const bam1_core_t & c = b->core;
    std::string qn(bam_get_qname(b));
    ss << qn << " " << d.barcode << " umi = " << d.umi << " pos = " << d.pos
       << " lft = " << c.pos << " cigar = ";
    uint32_t *cigar = bam_get_cigar(b);
    for (unsigned int i = 0; i < c.n_cigar; i++) {
        ss << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    char * flag = bam_flag2str(c.flag);
    ss << " flag = " << flag;
    free(flag);
    ss << " RE = " <<  d.xt << " seq ";
    auto s = bam_get_seq(b);
    for(int i = 0; i < b->core.l_qseq; i++){
        ss << "=ACMGRSVTWYHKDBN"[bam_seqi(s, i)];
    }

    return ss.str();
}

inline std::string debug_read(const bam1_t * b) {
    std::stringstream ss;
    const bam1_core_t & c = b->core;
    std::string qn(bam_get_qname(b));
    ss << qn << " " << "lft = " << c.pos << " cigar = ";
    uint32_t *cigar = bam_get_cigar(b);
    for (unsigned int i = 0; i < c.n_cigar; i++) {
        ss << bam_cigar_oplen(cigar[i]) << bam_cigar_opchr(cigar[i]);
    }
    char * flag = bam_flag2str(c.flag);
    ss << " flag = " << flag;
    free(flag);
    ss << " seq ";
    auto s = bam_get_seq(b);
    for(int i = 0; i < b->core.l_qseq; i++){
        ss << "=ACMGRSVTWYHKDBN"[bam_seqi(s, i)];
    }

    return ss.str();
}

class BamBuffer {
    template <typename T, typename R, typename P>
    friend class BamGeneReader;

    template <typename T, typename R, typename P>
    friend class BamGeneReaderFiltered;

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

        void set_starts(std::vector<unsigned int>::const_iterator start, std::vector<unsigned int>::const_iterator end){
            starts_.assign(start, end);
        }
    
        void add(BamDetail & d){
            assure_size();
            reads_[count_++]->copy(d);
        }

        size_t count() const {
            return count_;
        }

        size_t reserved() const {
            return reads_.size();
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
            if(!starts_.empty() && current_ < (starts_.size() - 1)){
                //std::cout << "Giving reads N = " << starts_[current_] << " - " << starts_[current_ + 1] 
                //    << " N = " << (starts_[current_ + 1] - starts_[current_]) 
                //    << "\n";
                rit start  = reads_.begin() + starts_[current_++];
                rit end = reads_.begin() + starts_[current_];
                ret = {start, end};
                return true;
            }
            ret = {reads_.end(), reads_.end()};
            return false;
        }

        unsigned int get_next_ranges(std::vector<rpair> & ranges, unsigned int N){
            ranges.clear();
            std::lock_guard<std::mutex> lock(mutex_);
            unsigned int count = 0;
            while(count < N && (current_ < starts_.size() - 1)){
                rit start  = reads_.begin() + starts_[current_++];
                rit end = reads_.begin() + starts_[current_];
                ranges.push_back({start, end});
                count++;
            }
            return count;
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

class ProcessorBase {
    public:

        virtual ~ProcessorBase(){

        }

        virtual bool operator()(BamDetail & read, unsigned int fno) = 0;

    protected:
        const TXIndex * index_;
        StrandMode      ST_;
};

class BamScSNVProcessor : public ProcessorBase {
    public: 
        BamScSNVProcessor(const TXIndex & index, StrandMode ST){
            index_ = &index;
            ST_ = ST;
        }

        bool operator()(BamDetail & read, unsigned int fno){
            char A = bam_aux2A(bam_aux_get(read.b, "RA"));
            if((A == 'E' || A == 'N') && read.b->core.qual == 255){
                read.gid = index_->gene_from_id(bam_aux2Z(bam_aux_get(read.b, "XG"))).gid;
                read.filenum = fno;
                return true;
            }
            return false;
        }

    private:
        const TXIndex * index_;
        StrandMode      ST_;
};


class BamCellRangerProcessor : public ProcessorBase{
    public: 
        BamCellRangerProcessor(const TXIndex & index, StrandMode ST){
            index_ = &index;
            ST_ = ST;
        }

        bool operator()(BamDetail & read, unsigned int fno){
            if(read.b->core.qual != 255) {
                return false;
            }

            if(bam_aux_get(read.b, "CB") == NULL || bam_aux_get(read.b, "UB") == NULL){
                return false;
            }
            auto ptr = bam_aux_get(read.b, "RE");
            char A = '?';
            if(ptr != NULL){
                A = bam_aux2A(ptr);
            }
            ptr = bam_aux_get(read.b, "MM");
            if(A == 'E' && ptr != NULL){
                unsigned int MM = bam_aux2i(ptr);
                if(MM == 1) {
                    return false;
                }
            }
            uint32_t gid = std::numeric_limits<uint32_t>::max();
            bool rev = bam_is_rev(read.b);
            char xs = ((ST_ == StrandMode::TAG_FWD && !rev) || (ST_ == StrandMode::TAG_REV && rev)) ? '+' : '-';
            read.filenum = fno;
            ptr = bam_aux_get(read.b, "GX");
            if(A == 'N' || (A == '?' && ptr == NULL)){
                auto const & ref = index_->ref(read.b->core.tid);
                if(!ref.empty()){
                    overlaps_.clear();
                    ref.tree.findOverlapping(read.b->core.pos, bam_endpos(read.b) - 1, overlaps_);
                    unsigned int count = 0;
                    for(size_t i = 0; i < overlaps_.size(); i++){
                        auto tgid = overlaps_[i].value;
                        const GeneEntry & gene = index_->gene(tgid);
                        if(gene.strand == xs){
                            gid = tgid;
                            count++;
                        }
                    }
                    if(count == 1){
                        const GeneEntry & gene = index_->gene(gid);
                        const char * gstr = gene.gene_id.c_str();
                        bam_aux_append(read.b, "XG", 'Z', strlen(gstr) + 1, reinterpret_cast<const uint8_t *>(gstr));
                    }else{
                        gid = std::numeric_limits<uint32_t>::max();
                    }
                }
            }else if(A == 'E' || (A == '?' && ptr != NULL)){
                char * XG = ptr == NULL ? NULL : bam_aux2Z(ptr);
                if(XG != NULL && strchr(XG, ';') == NULL){
                    gid = index_->gid_from_str(XG);
                    bam_aux_append(read.b, "XG", 'Z', strlen(XG) + 1, reinterpret_cast<const uint8_t *>(XG));
                }
            }else{
                return false;
            }

            if(gid != std::numeric_limits<uint32_t>::max()){
                bam_aux_append(read.b, "XS", 'A', 1, reinterpret_cast<const uint8_t*>(&xs));
                //read.b->core.flag &= ~BAM_FDUP;
                read.gid = gid;
                return true;
            }
            return false;
        }

    private:
        TXIndex::Ref::itree::intervalVector overlaps_;
};

}
