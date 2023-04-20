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
#include "transcript_align.hpp"
#include <cctype>
using namespace gwsc;

void TranscriptAlign::load(const std::string & prefix, unsigned int min_overhang){
    args_ = mem_opt_init();
    args_->b = ascore.mismatch;
    args_->a = ascore.match;
    bwa_fill_scmat(args_->a, args_->b, args_->mat);
    min_overhang_ = min_overhang;
    bidx_ = bwa_idx_load((prefix + "_bwa").c_str(), BWA_IDX_ALL);
    if(bidx_ == nullptr) {
        std::cout << "Index load failed\n";
        exit(1);
    }
}

void TranscriptAlign::align(AlignGroup & ad, const std::string & seq, unsigned int len) const {
    ad.transcript_score = 0;
    ad.transcript_ar = mem_align1(args_, bidx_->bwt, bidx_->bns, bidx_->pac, len, seq.c_str());
    if(ad.transcript_ar.n <= 0 || ad.transcript_ar.a[0].score < static_cast<int>(ascore.min_score)){
        free(ad.transcript_ar.a);
        ad.transcript_ar.a = NULL;
        ad.transcript_ar.n = 0;
        ad.transcript_ar.m = 0;
        return;
    }

    ad.transcript_score = ad.transcript_ar.a[0].score;
    /*
    for(size_t i = 1; i < ar.n; i++){
        if(ar.a[i].score > ar.a[i - 1].score){
            std::cout << "not sorted\n";
        }

    }
    */
}

void TranscriptAlign::get(AlignGroup & ad, size_t i, const std::string & seq, unsigned int len) const{
    auto & ar = ad.transcript_ar;
    //std::cout << "Transcript get " << i << " " << ar.n << " " << ad.acount << " " << " " << ad.alns.size() << "\n";
    mem_aln_t a;
    assert(i < ar.n);
    a = mem_reg2aln(args_, bidx_->bns, bidx_->pac, len, seq.c_str(), &ar.a[i]);
    //std::cout << "  a.score = " << ar.a[i].score << " a.truesc = " << ar.a[i].truesc << " reg score = " << a.score << "\n";
    if((ad.acount + 1) >= ad.alns.size()){
        ad.alns.resize(ad.acount + 5);
    }

    assert(ad.acount < ad.alns.size());

    make_align_(ad, a);
    auto & aln = ad.alns[ad.acount - 1];
    if(len < seq.size()){
        if(aln.rev){
            if(aln.cigar.front().op == Cigar::SOFT_CLIP){
                aln.cigar.front().len += seq.size() - len;
            }else{
                aln.cigar.insert(aln.cigar.begin(), CigarElement(seq.size() - len, Cigar::SOFT_CLIP));
                //aln.cigar.push_front(CigarElement(seq.size() - len, Cigar::SOFT_CLIP));
            }
        }else{
            if(aln.cigar.back().op == Cigar::SOFT_CLIP){
                aln.cigar.back().len += seq.size() - len;
            }else{
                aln.cigar.push_back(CigarElement(seq.size() - len, Cigar::SOFT_CLIP));
            }
        }
    }

    project(aln, ad.tcigar);
    free(a.cigar);
    gidx_->rescore(aln, seq);
}

void TranscriptAlign::make_align_(AlignGroup & ad, mem_aln_t & a) const {

    auto & aln = ad.alns[ad.acount++];
    /*
    printf("  Align: %c\t%s\t%d\t%ld\t%d\t", "+-"[a.is_rev], idx_->bns->anns[a.rid].name, a.rid, (long)a.pos, a.mapq);
    for (int k = 0; k < a.n_cigar; ++k) // print CIGAR
            printf("%d%c", a.cigar[k]>>4, "MIDSH"[a.cigar[k]&0xf]);
    printf("\t%d\t%d\t%s\n", p.score, a.NM, (char*)(a.cigar + a.n_cigar)); // print edit distance
    */


    int matches = 0, dels = 0;
    int rgt = a.pos;
    //size_t aidx = ad.acount++;
    //auto & aln = ad.alns[aidx];
    aln.reset();
    auto & cigar = aln.cigar;
    for(int k = 0; k < a.n_cigar; k++){
        uint32_t op  = a.cigar[k] & 0xF;
        if(op >= 3) a.cigar[k]++;
        cigar.push_back(CigarElement(a.cigar[k]));
        if(cigar.back().op == Cigar::DEL) {
            rgt += cigar.back().len;
            dels += cigar.back().len;
        }else if(cigar.back().op == Cigar::MATCH){
            matches += cigar.back().len;
            rgt += cigar.back().len;
        }
    }
    if(cigar.front().op == Cigar::HARD_CLIP) cigar.front().op = Cigar::SOFT_CLIP;
    if(cigar.back().op == Cigar::HARD_CLIP) cigar.back().op = Cigar::SOFT_CLIP;

    aln.NM = a.NM; aln.score = a.score;
    aln.total_bases = matches; aln.dels = dels; aln.mapq = a.mapq;
    aln.lft = a.pos; aln.rgt = rgt - 1;
    aln.txid = a.rid; aln.rev = a.is_rev;
    auto & t = idx_.transcript(aln.txid);
    aln.tid = t.tid;
    aln.gid = t.gid;
    aln.exonic = aln.total_bases + aln.dels;
    aln.valid_transcript = true;


}

void TranscriptAlign::project(AlignData & a, CigarString & tc) const {
    // Project a transcriptome only alignment to the genome
    auto const & tx = idx_.transcript(a.txid);
    bool tplus = tx.strand == '+';
    if(!tplus){
        a.rev = !a.rev;
    }
    if((smode_ == StrandMode::TAG_FWD && tplus == a.rev) || (smode_ == StrandMode::TAG_REV && tplus != a.rev)){
        a.atype=AlignType::ANTISENSE;
        a.xs = "+-"[tplus];
    }else{
        a.atype=AlignType::CDNA;
        a.xs = "-+"[tplus];
    }
    //std::cout << "  " << a.print(tx.transcript_id) << "\n";
    auto it = std::lower_bound(tx.rexons.begin(), tx.rexons.end(), a.lft, [](const PosEntry & p1, const uint32_t & p2){ return p1.rgt < p2; });
    if(it == tx.rexons.end()){
        std::cout << "Odd Genome Projection error?";
        exit(0);
    }
    unsigned int istart = it - tx.rexons.begin();
    tc.clear();
    auto cit = a.cigar.begin();
    unsigned int pos = a.lft;
    unsigned int clen = cit->len;
    if(cit->op == Cigar::SOFT_CLIP){
        tc.push_back(*cit);
        cit++;
        clen = cit->len;
    }
    unsigned int pcount = 0;
    auto cend = a.cigar.end();
    if(a.cigar.back().op == Cigar::SOFT_CLIP) cend--;

    while(it != tx.rexons.end() && it->lft <= a.rgt && cit != cend){
        unsigned int needed = std::min(it->rgt, a.rgt) - pos + 1;
        pcount++;
        //std::cout << "  pos = " << pos << " exon = " << it->lft << " - " << it->rgt << " needed = " << needed << "\n";
        while(cit != cend && needed > 0){
            //std::cout << "     Element = " << *cit << " clen = " << clen << " needed = " << needed << "\n";
            if(cit->op == Cigar::MATCH || cit->op == Cigar::DEL){
                if(clen > needed){
                    tc.push_back(CigarElement(needed, cit->op));
                    clen -= needed;
                    needed = 0;
                }else{
                    tc.push_back(CigarElement(clen, cit->op));
                    needed -= clen;
                    cit++;
                    clen = cit != cend ? cit->len : 0;
                }
            }else{
                tc.push_back(*cit);
                cit++;
                clen = cit != cend ? cit->len : 0;
            }
            //std::cout << "     After   = " << *cit << " clen = " << clen << " needed = " << needed << "\n";
        }
        // Create intron here
        pos = it->rgt + 1;
        it++;
        if(it != tx.rexons.end() && it->lft <= a.rgt){
            unsigned int isize = tx.isizes[(it - tx.rexons.begin()) - 1];
            tc.push_back(CigarElement(isize, Cigar::REF_SKIP));
        }
    }
    if(a.cigar.back().op == Cigar::SOFT_CLIP) tc.push_back(a.cigar.back());

    unsigned int nlft = 0;
    if(tx.strand == '+'){
        /*
        cout << tx.exons[istart].lft << " - " << tx.exons[istart].rgt << " | " 
             << tx.rexons[istart].lft << " - " << tx.rexons[istart].rgt << "\n";
         */
        nlft = tx.exons[istart].lft + (a.lft - tx.rexons[istart].lft);
    }else{
        /*
        cout << tx.exons[istart + pcount - 1].lft << " - " << tx.exons[istart + pcount - 1].rgt << " | " 
             << tx.rexons[istart + pcount - 1].lft << " - " << tx.rexons[istart + pcount - 1].rgt << "\n";
        */
        nlft = tx.exons[istart + pcount - 1].lft + (tx.rexons[istart + pcount - 1].rgt - a.rgt + 1) - 1;
        std::reverse(tc.begin(), tc.end());
        // Need to flip the strand too
    }

    auto l = nlft;
    //std::cout << "  New Cigar " << tc << " for transcript on " << tx.strand << " strand\n";
    cit = tc.begin();
    while(cit != tc.end()){
        if(cit->op == Cigar::MATCH || cit->op == Cigar::DEL){
            /*
            if(cit->op == Cigar::MATCH || cit->op == Cigar::DEL){
                std::cout << "    Alignment Region " << l << " - " << (l + cit->len - 1) << " " << *cit << "\n";
            }
            */
            l += cit->len;
        }else if(cit->op == Cigar::REF_SKIP){
            //std::cout << "    Ref Skip Region " << l << " - " << (l + cit->len - 1) << " " << *cit << "\n";
            l += cit->len;
        }
        cit++;
    }

    /*
    for(size_t i = istart; i < tx.rexons.size(); i++){
        cout << "       Exon " << tx.exons[i].lft << " - " << tx.exons[i].rgt
             << " | " << tx.rexons[i].lft << " - " << tx.rexons[i].rgt << "\n";
    }
    */



    //std::cout << "spliced = " << spliced << " fcount = " << fcount << " lcount = " << lcount << "\n";
    std::swap(a.cigar, tc);
    a.lft = nlft;
    a.rgt = l - 1;
    //a.rescore = true;
    trim_splices_(a);
}

template <typename T>
std::tuple<int, int, T> trim_to_splice(AlignData & a, T cit, T end, unsigned int mo){
    unsigned int bases = 0, del = 0, clip = 0;
    while(cit != end && cit->op != Cigar::REF_SKIP){
        if(cit->op == Cigar::MATCH){
            bases += cit->len;
        }else if(cit->op == Cigar::INS || cit->op == Cigar::SOFT_CLIP){
            clip += cit->len;
        }else if(cit->op == Cigar::DEL){
            del += cit->len;
        }
        cit++;
    }
    // No splice to trim off
    if(cit != end && bases <= mo) {
        a.total_bases -= bases;
        a.exonic -= (bases + del);
        return {bases + del + cit->len, clip + bases, ++cit};
    }
    return {0, 0, end};
}

void TranscriptAlign::trim_splices_(AlignData & a) const {
    {
        auto res = trim_to_splice(a, a.cigar.begin(), a.cigar.end(), min_overhang_);
        int bases = std::get<0>(res);
        int clipped = std::get<1>(res);
        auto it = std::get<2>(res);
        if(bases > 0){
            a.rescore = true;
            a.lft += bases;
            a.cigar.erase(a.cigar.begin() + 1, it);
            a.cigar.front().op = Cigar::SOFT_CLIP;
            a.cigar.front().len = clipped;
        }
    }
    {
        auto res = trim_to_splice(a, a.cigar.rbegin(), a.cigar.rend(), min_overhang_);
        int bases = std::get<0>(res);
        int clipped = std::get<1>(res);
        auto it = std::get<2>(res);
        if(bases > 0){
            a.rescore = true;
            a.rgt -= bases;
            a.cigar.erase(it.base() + 1, a.cigar.end());
            a.cigar.back().op = Cigar::SOFT_CLIP;
            a.cigar.back().len = clipped;
        }
    }
}
