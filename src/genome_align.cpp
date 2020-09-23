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
#include "genome_align.hpp"
#include <iostream>
#include <string>
#include <stdlib.h>
using namespace gwsc;

void GenomeAlign::load(const std::string & prefix, unsigned int min_overhang) {
    args_ = mem_opt_init();
    args_->b = ascore.mismatch;
    args_->a = ascore.match;
    bwa_fill_scmat(args_->a, args_->b, args_->mat);
    bidx_ = bwa_idx_load(prefix.c_str(), BWA_IDX_ALL);
    min_overhang_ = min_overhang;
    if(bidx_ == nullptr) {
        std::cout << "Index load failed\n";
        exit(1);
    }

    //Check to make sure the genome and the transcript index make sense
    for(int32_t i = 0; i < bidx_->bns->n_seqs; i++){
        auto const & r = idx_.ref(i);
        if(bidx_->bns->anns[i].name != r.name){
            std::cout << "Reference mismatch for tid = " << i << " it appears the genome was built with a different reference than the transcript index (" << bidx_->bns->anns[i].name << " vs " << r.name << ")\n";
            exit(EXIT_FAILURE);
        }
    }
}


void GenomeAlign::align(AlignGroup & ad, const std::string & seq, unsigned int len) const {
    ad.genome_score = 0;
    ad.genome_ar = mem_align1(args_, bidx_->bwt, bidx_->bns, bidx_->pac, len, seq.c_str());
    if(ad.genome_ar.n <= 0 || ad.genome_ar.a[0].score < static_cast<int>(ascore.min_score)){
        free(ad.genome_ar.a);
        ad.genome_ar.a = NULL;
        ad.genome_ar.n = 0;
        ad.genome_ar.m = 0;
        return;
    }
    
    ad.genome_score = ad.genome_ar.a[0].score;
}

void GenomeAlign::get(AlignGroup & ad, size_t i, const std::string & seq, unsigned int len) const {
    auto & ar = ad.genome_ar;
    assert(i < ar.n);
    mem_aln_t a;
    a = mem_reg2aln(args_, bidx_->bns, bidx_->pac, len, seq.c_str(), &ar.a[i]);
    //std::cout << "  a.score = " << ar.a[i].score << " a.truesc = " << ar.a[i].truesc << " reg score = " << a.score << "\n";

    if((ad.acount + 1) >= ad.alns.size()){
        ad.alns.resize(ad.acount + 5);
    }

    //std::cout << "len = " << len << " score = " << ar.a[i].score << " " << ar.a[i].truesc << "\n";
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
    ad.toverlaps.clear();
    //std::cout << "Annotate cigar " << aln.cigar << " rescore = " << aln.rescore << "\n";
    //std::cout << "Pre = " << aln.print() << "\n";
    annotate_(aln, ad.toverlaps, ad.tbases);
    //std::cout << "After = " << aln.print() << "\n";
    //verify(aln, seq, "genome after verification");
    free(a.cigar);
    rescore(aln, seq);
}



void GenomeAlign::make_align_(AlignGroup & ad, mem_aln_t & a) const {
    /*
    printf("  Align: %c\t%s\t%d\t%ld\t%d\t", "+-"[a.is_rev], idx_->bns->anns[a.rid].name, a.rid, (long)a.pos, a.mapq);
    for (int k = 0; k < a.n_cigar; ++k) // print CIGAR
            printf("%d%c", a.cigar[k]>>4, "MIDSH"[a.cigar[k]&0xf]);
    printf("\t%d\t%d\t%s\n", p.score, a.NM, (char*)(a.cigar + a.n_cigar)); // print edit distance
    */


    int matches = 0, dels = 0;
    int rgt = a.pos;
    size_t aidx = ad.acount++;
    auto & aln = ad.alns[aidx];
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
    //std::cout << "cigar done\n";
    aln.total_bases = matches; aln.dels = dels; aln.mapq = a.mapq;
    aln.lft = a.pos; aln.rgt = rgt - 1;
    aln.tid = a.rid; aln.rev = a.is_rev;
    aln.score = a.score; aln.NM = a.NM;
    aln.xs = ((smode_ == StrandMode::TAG_FWD && !aln.rev) || (smode_ == StrandMode::TAG_REV && aln.rev)) ? '+' : '-';

    auto const & ref = idx_.ref(aln.tid);
    //std::cout << "Original cigar " << aln.cigar << " score before trim = " << a.score << "\n";
    trim_splices_(aln, ref, ad.toverlaps);
}

template <typename T>
std::tuple<int, int, T> trim_(AlignData & a, T cit, T end, unsigned int needed){
    unsigned int dels = 0, bases = 0;
    unsigned int clipped = 0;
    //std::cout << "    Needed = " << needed << "\n";
    while(cit != end){
        if(cit->op == Cigar::INS || cit->op == Cigar::SOFT_CLIP){
            clipped += cit->len;
            //std::cout << "    Clipped += " << (*cit) << "\n";
        }else if(cit->op == Cigar::MATCH){
            if(cit->len > needed){
                //std::cout << "    Clipped += " << (*cit) << " needed = " << needed << "\n";
                clipped += needed;
                cit->len -= needed;
                bases += needed;
                needed = 0;
                break;
            }else{
                //std::cout << "    Clipped += " << (*cit) << " needed = " << needed << "\n";
                clipped += cit->len;
                needed -= cit->len;
                bases += cit->len;
            }
        }else if(cit->op == Cigar::DEL){
            // Take the whole deletion
            dels += cit->len;
            //std::cout << "    Clipped += " << (*cit) << " needed = " << needed << "\n";
            if(cit->len > needed){
                needed = 0;
                cit++;
                break;
            }else{
                needed -= cit->len;
            }
        }
        cit++;
    }

    while(cit->op == Cigar::INS || cit->op == Cigar::DEL){
        if(cit->op == Cigar::INS){
            clipped += cit->len;
        }else{
            dels += cit->len;
        }
        cit++;
    }

    a.total_bases -= bases;
    a.dels -= dels;
    return {bases + dels, clipped, cit};
}

void GenomeAlign::trim_splices_(AlignData & a, const TXIndex::Ref & ref, TXIndex::Ref::itree::intervalVector & overlaps) const {

    unsigned int trim_lft = std::numeric_limits<unsigned int>::max();
    unsigned int trim_rgt = std::numeric_limits<unsigned int>::max();
    const TXIndex::Ref::itree * lsites = nullptr, *rsites = nullptr;
    if(a.xs == '+'){
        lsites = &ref.lp_splices;
        rsites = &ref.rp_splices;
    }else{
        lsites = &ref.lm_splices;
        rsites = &ref.rm_splices;
    }

    overlaps.clear();
    lsites->findOverlapping(a.lft, a.lft, overlaps);
    if(!overlaps.empty()){
        trim_lft = 0;
        for(auto & o:overlaps){
            trim_lft = std::max(trim_lft, o.stop);
        }
    }

    overlaps.clear();
    rsites->findOverlapping(a.rgt, a.rgt, overlaps);
    if(!overlaps.empty()){
        for(auto & o:overlaps){
            trim_rgt = std::min(trim_rgt, o.start);
        }
    }

    if(trim_lft < std::numeric_limits<unsigned int>::max()){
        unsigned int needed = trim_lft - a.lft + 1;
        //std::cout << "    Trim Left: " << a.tid << ": " << a.lft << " - " << a.rgt << " Needed = " << needed << " cigar = " << a.cigar << " trim_lft = " << trim_lft << "\n";
        auto res = trim_(a, a.cigar.begin(), a.cigar.end(), needed);
        auto it = std::get<2>(res);
        int clipped = std::get<1>(res);
        a.rescore = true;
        a.lft += std::get<0>(res);
        //std::cout << "  Lft Cigar = " << a.cigar << " last = " << (*it) << " clipped = " << clipped << "\n";;
        if(it == a.cigar.begin()){
            a.cigar.insert(a.cigar.begin(), CigarElement(clipped, Cigar::SOFT_CLIP));
        }else{
            a.cigar.erase(std::next(a.cigar.begin()), it);
            a.cigar.front().op = Cigar::SOFT_CLIP;
            a.cigar.front().len = clipped;
        }
        //std::cout << "  After erase: " << a.cigar << "\n";
    }

    if(trim_rgt < std::numeric_limits<unsigned int>::max()){
        unsigned int needed = a.rgt - trim_rgt + 1;
        //std::cout << "    Trim Right: " << a.tid << ": " << a.lft << " - " << a.rgt << " Needed = " << needed << " cigar = " << a.cigar << " trim_rgt = " << trim_rgt << "\n";
        auto res = trim_(a, a.cigar.rbegin(), a.cigar.rend(), needed);
        auto it = std::get<2>(res);
        int clipped = std::get<1>(res);
        a.rescore = true;
        a.rgt -= std::get<0>(res);
        //std::cout << "Rgt Cigar = " << a.cigar << " last = " << (*it) << " clipped = " << clipped << "\n";;
        auto t = it.base();
        a.cigar.erase(t, a.cigar.end());
        a.cigar.push_back(CigarElement(clipped, Cigar::SOFT_CLIP));
        //std::cout << "  After erase: " << a.cigar << "\n";
    }
}

void GenomeAlign::calculate_bases_(AlignData & a, const GeneEntry & g, AlignBases & b) const {
    unsigned int overlap = (std::min(a.rgt, g.rgt) - std::max(a.lft, g.lft) + 1);
    unsigned int tb = a.total_bases + a.dels;
    b.intergenic = tb - overlap;
    auto it = std::lower_bound(g.introns.begin(), g.introns.end(), g.lft, [](const PosEntry & p, unsigned int lft) { return p.rgt < lft; });
    b.intronic = 0;
    while(it != g.introns.end() && it->lft < a.rgt){
        if(it->lft <= a.rgt && a.lft <= it->rgt){
            b.intronic += (std::min(a.rgt, it->rgt) - std::max(a.lft, it->lft) + 1);
            //cout << "    Intron index = " << (it - g.introns.begin()) << " " << it->lft << " - " << it->rgt << " N = " << it->length() << " overlap = " << overlap << "\n";
        }
        it++;
    }
    b.exonic = overlap - b.intronic;
}

void GenomeAlign::annotate_(AlignData & a, TXIndex::Ref::itree::intervalVector & overlaps, 
        std::vector<AlignBases> & bases) const {
    auto const & ref = idx_.ref(a.tid);
    //std::cout << "annotate " << a.print() << "\n";
    if(ref.empty()){
        a.gid = std::numeric_limits<uint32_t>::max();
        a.intergenic = a.total_bases + a.dels;
        a.atype = AlignType::INTERGENIC;
        return;
    }
    overlaps.clear();
    ref.tree.findOverlapping(a.lft, a.rgt, overlaps);
    if(overlaps.empty()){
        a.gid = std::numeric_limits<uint32_t>::max();
        a.intergenic = a.total_bases + a.dels;
        a.atype = AlignType::INTERGENIC;
        return;
    }
    
    bases.clear();
    bases.resize(overlaps.size());

    int mexonic = -1;
    int mindex = -1;

    int aexonic = -1;
    int aindex = -1;
    for(size_t i = 0; i < overlaps.size(); i++){
        auto ov = overlaps[i];
        auto const & g = idx_.gene(ov.value);
        calculate_bases_(a, g, bases[i]);
        if(g.strand == a.xs && mexonic < (int)bases[i].exonic){
            mexonic = bases[i].exonic;
            mindex = i;
        }else if(g.strand != a.xs && aexonic < (int)bases[i].exonic){
            aexonic = bases[i].exonic;
            aindex = i;
        }
    }

    unsigned int ovc = 0, ova = 0;
    for(int i = 0; i < (int)overlaps.size(); i++){
        auto ov = overlaps[i];
        auto const & g = idx_.gene(ov.value);
        auto const & b = bases[i];
        //std::cout << "  " << g.gene_name << " bases e = " << b.exonic << " i = " << b.intronic << " n = " << b.intergenic << " diff = " << (mexonic - b.exonic) << "\n";
        if(g.strand == a.xs && (mexonic - b.exonic) < 5){
            a.gids.push_back(ov.value);
            ovc++;
        }else if(g.strand != a.xs && (aexonic - b.exonic) < 5){
            ova++;
        }
    }
    //std::cout << "  ova = " << ova << " ovc = " << ovc << "\n";

    a.exonic = 0;
    a.intronic = 0;
    a.intergenic = 0;
    if(ova == 0 && ovc == 0){
        //std::cout << "    double\n";
        // No alignments that overlap a gene that isn't already mapped to a transcript
        a.atype = AlignType::UNKNOWN;
    }else if(ovc == 0 && ova > 1){
        //Antisense alignment to more than one gene
        a.atype = AlignType::ANTISENSE;
    }else if(ovc > 1){
        a.atype = AlignType::AMBIGUOUS;
    }else if(ovc == 0 && ova == 1){
        a.atype = AlignType::ANTISENSE;
        a.gid = idx_.gene(overlaps[aindex].value).gid;
        a.exonic = bases[aindex].exonic;
        a.intergenic = bases[aindex].intergenic;
        a.intronic = bases[aindex].intronic;
        a.txid = std::numeric_limits<uint32_t>::max();
    }else{
        a.gid = idx_.gene(overlaps[mindex].value).gid;
        a.exonic = bases[mindex].exonic;
        a.intergenic = bases[mindex].intergenic;
        a.intronic = bases[mindex].intronic;
        a.txid = std::numeric_limits<uint32_t>::max();
        if(a.exonic >= 25 && a.intronic < 5){
            a.atype = AlignType::CDNA;
        }else if(a.intronic >= 5){
            a.atype = AlignType::INTRON;
        }else{
            a.atype = AlignType::INTERGENIC;
            a.gid = std::numeric_limits<uint32_t>::max();
        }
    }
}

struct BWAit{
    BWAit(const BWAit & c) : pac(c.pac), astart(c.astart), aend(c.aend), pos(c.pos)
    {

    }

    BWAit & operator=(const BWAit & c){
        pac = c.pac;
        astart = c.astart;
        aend = c.aend;
        pos = c.pos;
        return *this;
    }

    BWAit(const bntseq_t * bnt, const uint8_t * pac, unsigned int lft, int tid) : pac(pac){
        pos = bnt->anns[tid].offset + lft;
        aend = bnt->ambs + bnt->n_holes;
        astart = std::lower_bound(bnt->ambs, aend, pos, [](const bntamb1_t & a, unsigned int pos){
            return (a.offset + a.len) < pos;    
        });

    }

    char operator*() const{
        /*
        if(astart != aend){
            std::cout << "ofset = " << lft << " ambig area " << astart->offset << " - " << astart->len << "\n";
        }else{
            std::cout << "astart at end offset = " << lft << " last = " << (aend - 1)->offset << " " << (aend - 1)->len << "\n";
        }
        */
        if(astart != aend && pos >= astart->offset && pos < (astart->offset + astart->len)){
            //std::cout << "Inside ambiguous\n";
            return astart->amb;
        }
        char b = "ACGT"[(pac[pos>>2] >> ((~(pos)&3)<<1)) & 3];
        return b;
    }


    BWAit & operator+=(int i){
        pos += i;
        while(astart != aend && (astart->offset + astart->len) < pos) astart++;
        return *this;
    }

    BWAit operator++(int){
        BWAit tmp(*this);
        ++(*this);
        return tmp;
    }

    BWAit & operator++(){
        pos++;
        while(astart != aend && (astart->offset + astart->len) < pos) astart++;
        return *this;
    }

    const uint8_t  * pac;
    bntamb1_t      * astart;
    bntamb1_t      * aend;
    unsigned int     pos;
};

template<typename T>
void print_(BWAit rit, AlignData & a, T qlft, mem_opt_t * args_, const std::string & msg){
    auto lft = a.lft;
    std::string qy, ar, re;
    int score = 0, NM = 0;
    std::cout << a.print() << "\n";
    for(auto & c : a.cigar){
        if(c.op == Cigar::MATCH){
            for(size_t i = 0; i < c.len; i++, lft++, qlft++, rit++){
                char b = *rit;
                char q = *qlft;
                if(a.rev) {
                    q = reverse_cmpl_[static_cast<int>(q)];
                }
                qy.push_back(q);
                re.push_back(b);
                ar.push_back(b == q ? '|' : '*');
                if(q == 'N' || b == 'N'){
                    score--;
                    NM++;
                }else{
                    score += (b == q ? args_->a : -args_->b);
                    NM += (b != q);
                }

            }
        }else if(c.op == Cigar::DEL){
            score -= (args_->o_del + c.len * args_->e_del);
            NM += c.len;
            for(size_t i = 0; i < c.len; i++, lft++, rit++){
                char b = *rit;
                re.push_back(b);
                ar.push_back('D');
                qy.push_back('-');
            }
        }else if(c.op == Cigar::SOFT_CLIP){
            for(size_t i = 0; i < c.len; i++, qlft++){
                char q = *qlft;
                if(a.rev) {
                    q = reverse_cmpl_[static_cast<int>(q)];
                }
                qy.push_back(q);
                re.push_back(' ');
                ar.push_back(' ');
            }
        }else if(c.op == Cigar::INS){
            score -= (args_->o_ins + c.len * args_->e_ins);
            NM += c.len;
            for(size_t i = 0; i < c.len; i++, qlft++){
                char q = *qlft;
                if(a.rev) {
                    q = reverse_cmpl_[static_cast<int>(q)];
                }
                re.push_back(' ');
                qy.push_back(q);
                ar.push_back('I');
            }
        }else if(c.op == Cigar::REF_SKIP){
            lft += c.len;
            rit += c.len;
            qy.push_back('^');
            ar.push_back('^');
            re.push_back('^');
        }
    }
    std::cout << msg << " NM Error original = " << a.NM << " vs " << NM << " rescored = " << a.rescore << "\n";
    std::cout << "  qry: " << qy << "\n";
    std::cout << "  aln: " << ar << "\n";
    std::cout << "  ref: " << re << "\n";
    std::cout << "Score before = " << a.score << " score after = " << score << "\n\n";
}

template<typename T>
void rescore_(BWAit rit, AlignData & a, T qlft, mem_opt_t * args_, bool check, const std::string & msg = ""){
    auto lft = a.lft;
    int score = 0, NM = 0;
    T qorig = qlft;
    auto rit_cp = rit;
    for(auto & c : a.cigar){
        if(c.op == Cigar::MATCH){
            for(size_t i = 0; i < c.len; i++, lft++, qlft++, rit++){
                char b = *rit;
                char q = *qlft;
                if(a.rev) {
                    q = reverse_cmpl_[static_cast<int>(q)];
                }
                if(b == q && q != 'N' && b != 'N'){
                    score += args_->a;
                }else{
                    score-=args_->b;
                    NM++;
                }
            }
        }else if(c.op == Cigar::DEL){
            NM += c.len;
            score -= (args_->o_del + c.len * args_->e_del);
            lft += c.len;
            rit += c.len;
        }else if(c.op == Cigar::SOFT_CLIP){
            qlft += c.len;
        }else if(c.op == Cigar::INS){
            NM += c.len;
            score -= (args_->o_ins + c.len * args_->e_ins);
            qlft += c.len;
        }else if(c.op == Cigar::REF_SKIP){
            lft += c.len;
            rit += c.len;
        }
    }
    if((check && (int)a.NM != NM)){
        print_(rit_cp, a, qorig, args_, msg);
    }
    a.score = score;
    a.NM = NM;
}

void GenomeAlign::verify(AlignData & a, const std::string & seq, const std::string & msg) const {
    //auto ref = bidx_->bns->anns[a.tid];
    //auto offset = ref.offset;
    //auto pacseq = bidx_->pac;
    BWAit rit(bidx_->bns, bidx_->pac, a.lft, a.tid);
    if(a.rev){
        rescore_(rit, a, seq.rbegin(), args_, true, msg);
        //if(a.score < 20 || std::abs(oscore - a.score) > 20){
        //    rescore_(rit, a, seq.rbegin(), args_, true);
        //}
    }else{
        rescore_(rit, a, seq.begin(), args_, true, msg);
        //if(a.score < 20 || std::abs(oscore - a.score) > 20){
        //    rescore_(rit, a, seq.begin(), args_, true);
        //}
    }

}
void GenomeAlign::rescore(AlignData & a, const std::string & seq) const {
    //auto ref = bidx_->bns->anns[a.tid];
    //auto offset = ref.offset;
    //auto pacseq = bidx_->pac;
    BWAit rit(bidx_->bns, bidx_->pac, a.lft, a.tid);
    if(a.rev){
        rescore_(rit, a, seq.rbegin(), args_, false);
        //if(a.score < 20 || std::abs(oscore - a.score) > 20){
        //    rescore_(rit, a, seq.rbegin(), args_, true);
        //}
    }else{
        rescore_(rit, a, seq.begin(), args_, false);
        //if(a.score < 20 || std::abs(oscore - a.score) > 20){
        //    rescore_(rit, a, seq.begin(), args_, true);
        //}
    }
}
