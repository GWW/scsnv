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
#include <functional>
#include <iostream>
#include <queue>
#include <tuple>
#include "../util/misc.hpp"
#include "../htslib/htslib/hts.h"
#include "../htslib/htslib/sam.h"

namespace gwsc{

class BamMerger {

    struct MergeNode{
        MergeNode() {
            read = bam_init1();
        }

        ~MergeNode(){
            if(read != NULL)
                bam_destroy1(read);
            read = NULL;
        }

        bam1_t *     read = NULL;
        uint64_t     hash;
        unsigned int fno;
    };

    struct MergePtrCmp{
        bool operator()(MergeNode *p1, MergeNode * p2) const {
            return p1->hash > p2->hash;
        }
    };

    public:
        BamMerger(){
        }

        ~BamMerger(){

            for(auto & p : files_){
                if(p.second != nullptr){
                    bam_hdr_destroy(p.second);
                    p.second = nullptr;
                }
            }
        }

        template <typename IT>
        void add_bams(IT start, IT end);

        bam1_t * next(bam1_t * read = nullptr);

        const bam_hdr_t * header() const {
            if(files_.empty()) return nullptr;
            return files_.front().second;
        }

        size_t total() const {
            return total_;
        }

    private:
        std::vector<std::pair<samFile*, bam_hdr_t*>>                          files_;
        std::priority_queue<MergeNode*, std::vector<MergeNode*>, MergePtrCmp> pq_;
        size_t                                                                total_ = 0;


};

inline bam1_t * BamMerger::next(bam1_t * read) {
    if(pq_.empty()) return nullptr;
    auto n = pq_.top();
    pq_.pop();
    auto & p = files_[n->fno];
    if(read == nullptr) read = bam_init1();
    if(bam_copy1(read, n->read) == NULL) {
        std::cerr << "Error copying bam record\n"; exit(1);
    }

    if(p.first != NULL && sam_read1(p.first, p.second, n->read) >= 0){
        total_++;
        int32_t tid = n->read->core.tid == -1 ? p.second->n_targets : n->read->core.tid;
        n->hash = (static_cast<uint64_t>(tid) <<32) | (n->read->core.pos+1)<<1 | bam_is_rev(n->read);
        pq_.push(n);
    }else{
        //std::cout << "Done file " << n->fno << " done\n";
        sam_close(p.first);
        p.first = NULL;
        delete n;
    }
    return read;
}

template <typename IT>
inline void BamMerger::add_bams(IT start, IT end){
    while(start != end){
        samFile * sf = sam_open(start->c_str(), "r");
        bam_hdr_t * bh = sam_hdr_read(sf);
        files_.push_back({sf, bh});
        auto node = new MergeNode();
        node->fno = files_.size() - 1;
        if(sam_read1(sf, bh, node->read) < 0){
            std::cerr << "Error reading file " << files_.back().first->fn << "\n";
        }
        total_++;
        int32_t tid = node->read->core.tid == -1 ? bh->n_targets : node->read->core.tid;
        node->hash = (static_cast<uint64_t>(tid) << 32) | (node->read->core.pos+1)<<1 | bam_is_rev(node->read);
        //std::cout << "  Read " << node->fno << " " << node->read->core.tid << " pos = " << node->read->core.pos << "\n";
        pq_.push(node);
        start++;
    }
}

}
