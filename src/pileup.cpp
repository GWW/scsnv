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
#include "pileup.hpp"
#include <iostream>
#include <limits>

using namespace gwsc;
using namespace std;


void Pileup::build_reads(BamBuffer::rpair range) {
    count_ = 0;
    idx_ = 0;
    for(auto it = range.first; it != range.second; it++){
        if((count_ + 1) >= reads_.size()){
            for(size_t i = 0; i < 16; i++){
                reads_.push_back(new PileupRead());
            }
        }

        if((!dups_ && (*it)->dup)) continue;
        reads_[count_++]->init(*(*it));
    }
}

//template <typename T>
bool Pileup::next(std::vector<PileupOut> & out){
    if(buffer_.empty() && (reads_.empty() || idx_ >= count_)){
        return false;
    }
    get_reads_();
    pileup_reads_(out);
    return true;
}

void Pileup::pileup_reads_(std::vector<PileupOut> & out){
    auto it = buffer_.begin();
    unsigned int nextp = (idx_ < count_ && reads_[idx_]->tid() == ctid_) 
                              ? reads_[idx_]->rpos : std::numeric_limits<unsigned int>::max();
    out.clear();
    while(it != buffer_.end()){
        auto & r = *(*it);
        assert(r.rpos >= cpos_);
        if(r.rpos == cpos_){
            out.push_back(PileupOut(r));
            if(!r.next()){
                it = buffer_.erase(it);
            }else{
                nextp = std::min(r.rpos, nextp);
                it++;
            }
        }else{
            nextp = std::min(r.rpos, nextp);
            it++;
        }

    }

    pos = cpos_;
    tid = ctid_;

    assert(nextp > cpos_);
    if(nextp <= cpos_) {
        std::cout << "   There seems to be an issue with the sorting of your bam file tid = " << tid << " ctid_ = " << ctid_ << " nextp = " << nextp << " cpos = " << cpos_ << "\n";
        exit(0);
    }
    cpos_ = nextp;
}

void Pileup::get_reads_(){
    // No new reads to fetch
    if(idx_ >= count_){
        return;
    }
    if(buffer_.empty()){
        ctid_ = reads_[idx_]->tid();
        cpos_ = reads_[idx_]->rpos;
    }
    while(idx_ < count_ && reads_[idx_]->tid() == ctid_ && reads_[idx_]->rpos == cpos_){
        buffer_.push_back(reads_[idx_++]);
    }
}


void PileupRead::init(BamDetail & d){
    cigar.clear();
    ibases = 0;
    dbases = 0;
    abases = 0;
    this->d = &d;
    auto * ptr = bam_aux_get(d.b, "NM");
    if(ptr == NULL){
        ptr = bam_aux_get(d.b, "nM");
    }
    if(ptr == NULL){
        NM = 0;
    }else{
        NM = bam_aux2i(ptr);
    }

    qpositions.clear();
    auto xc = bam_aux_get(d.b, "XC");
    char * mstr = (xc == NULL ? NULL : bam_aux2Z(xc));

    strand = bam_aux2A(bam_aux_get(d.b, "XS"));
    uint32_t * cig = bam_get_cigar(d.b);
    unsigned int cs = 0;
    unsigned int qlft = 0;
    int qstart = 0;
    if(bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP){
        auto op = bam_cigar_oplen(cig[0]);
        qlft = op;
        qstart = op;
        cs++;
    }

    //std::cout << "Need to remove qpositions that are internal to a read block\n";

    qpositions.push_back(qstart);

    unsigned int rlft = d.lft;
    size_t sidx = 0;
    for(size_t j = cs; j < d.b->core.n_cigar; j++){
        cigar.push_back(CigarElement(cig[j]));
        auto op = cigar.back().op;
        auto len = cigar.back().len;
        switch(op){
            case Cigar::INS:
                ibases += len;
                qlft += len;
                break;
            case Cigar::DEL:
                rlft += len;
                dbases += len;
                break;
            case Cigar::SOFT_CLIP:
                cigar.pop_back();
                break;
            case Cigar::MATCH:
                rlft += len;
                qlft += len;
                abases += len;
                break;
            case Cigar::REF_SKIP:
                if(mstr != NULL && mstr[sidx] == '1'){
                    qpositions.push_back(qlft - 1);
                    qpositions.push_back(qlft);
                }
                sidx++;
                break;
            default:
                break;
        }
    }

    qpositions.push_back(qlft  - 1);

    qdist = std::numeric_limits<unsigned int>::max();
    for(auto & q : qpositions){
        unsigned int d = q < qpos ? qpos - q : q - qpos;
        qdist = std::min(qdist, d);
    }

    NM -= (ibases + dbases);
    cpos = 0;
    rpos = d.lft;
    qpos = qstart;
    rend = d.rgt + 1;
    it = cigar.begin();
}

bool PileupRead::next() {
    cpos++;
    qpos++;
    if(it != cigar.end() && cpos >= it->len){
        it++;
        while(it != cigar.end() && it->op != Cigar::MATCH){
            switch(it->op){
                case Cigar::DEL: case Cigar::REF_SKIP:
                    rpos += it->len;
                    break;
                case Cigar::INS:
                    qpos += it->len;
                    break;
                default:
                    break;
            }
            it++;
        }
        cpos = 0;
    }
    rpos++;
    qdist = std::numeric_limits<unsigned int>::max();
    for(auto & q : qpositions){
        unsigned int d = q < qpos ? qpos - q : q - qpos;
        //std::cout << " q = " << q << " d = " << d << ", ";
        qdist = std::min(qdist, d);
    }
    //std::cout << "  Read qpos = " << qpos << " qdist = " << qdist << "\n";
    return rpos < rend;
}

