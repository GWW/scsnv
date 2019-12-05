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
            // TODO: **This assuption may be wrong for RNA-seq reads?**
            //if(debug_) cout << "    Skip: " << *it << "\n";
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
        NM = bam_aux2i(bam_aux_get(d.b, "nM"));
    }else{
        NM = bam_aux2i(ptr);
    }
    strand = bam_aux2A(bam_aux_get(d.b, "XS"));
    uint32_t * cig = bam_get_cigar(d.b);
    unsigned int cs = 0;
    unsigned int qlft = 0;
    if(bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP){
        qlft += bam_cigar_oplen(cig[0]);
        cs++;
    }
    for(size_t j = cs; j < d.b->core.n_cigar; j++){
        cigar.push_back(CigarElement(cig[j]));
        if(cigar.back().op == Cigar::INS){
            ibases += cigar.back().len;
        }else if(cigar.back().op == Cigar::DEL){
            dbases += cigar.back().len;
        }else if(cigar.back().op == Cigar::SOFT_CLIP){
            cigar.pop_back();
        }else if(cigar.back().op == Cigar::MATCH){
            abases += cigar.back().len;
        }
    }

    NM -= (ibases + dbases);
    cpos = 0;
    rpos = d.lft;
    qpos = qlft;
    rend = d.rgt + 1;
    it = cigar.begin();
    //cout << "Init: " << read.qname << " " << read.tname << ": " << read.lft << " - " << read.rgt << " cigar: " << read.cigar << " rpos = " << rpos << " qpos = " << qpos << "\n";
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
    return rpos < rend;
}

