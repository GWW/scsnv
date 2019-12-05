#include "collapse_worker.hpp"
#include <vector>
#include <numeric>
#include <algorithm>
#include <sstream>

#include "../scmap/align_aux.hpp"
#include "../scmap/aux.hpp"
#include "../util/sequence.hpp"

using namespace gwsc;

void CollapseWorker::operator()() {
    BamBuffer::rpair range;
    while(buffer_->get_next(range)){
        for(auto it = range.first; it != range.second; it++){
            BamDetail & d = *(*it);
            d.corrected = false;
            d.xt = bam_aux2A(bam_aux_get(d.b, "RE"));
            tmp_ = bam_aux2Z(bam_aux_get(d.b, "CB"));
            if(cellranger_){
                std::size_t pos = tmp_.find("-"); 
                tmp_ = tmp_.substr(0, pos);
            }
            d.barcode = bhash_(tmp_);
            tmp_ = bam_aux2Z(bam_aux_get(d.b, "UB"));
            seq2int<gwsc::ADNA4, uint32_t>(tmp_, d.umi);
            d.make_hash();
            d.lft = d.b->core.pos;
            d.rgt = bam_endpos(d.b) - 1;

            if(bam_is_rev(d.b)){
                d.pos = d.rgt;
            }else{
                d.pos = d.lft;
            }
        }
        sort_reads_index(range.first, range.second, sidx_);
        size_t istart = 0;
        uint32_t lgid = (*(sidx_.front() + range.first))->gid;
        for(size_t i = 0; i < sidx_.size(); i++){
            auto it = range.first + sidx_[i];
            BamDetail & d = *(*it);
            if(d.gid != lgid){
                process_gene_(istart, i, range.first);
                istart = i;
                lgid = d.gid;
            }
        }
        if(istart < sidx_.size()){
            process_gene_(istart, sidx_.size(), range.first);
        }
    }
}

void CollapseWorker::process_gene_(size_t start, size_t end, BamBuffer::rit it){
    uint32_t gid = (*(it + sidx_[start]))->gid;
    auto ret = ghash_.find(gid);
    uhash_.clear();
    size_t tot = 0, cor = 0;
    //std::cout << " Processing gid = " << gid << " range = " << start << " - " << end << "\n";
    if(ret != ghash_.end()){
        size_t gstart = ret->second;
        while(gstart < umap_.size() && umap_[gstart].gene_id == gid){
            //std::cout << "    uhash gid = " << umap_[gstart].gene_id << " barcode = " << umap_[gstart].barcode << " umi = " << umap_[gstart].umi_from << "\n";
            uint64_t k = (static_cast<uint64_t>(umap_[gstart].barcode) << 32) | umap_[gstart].umi_from;
            uhash_[k] = umap_[gstart].umi_to;
            gstart++;
            tot++;
        }
    }

    // Correct the UMIs
    for(size_t i = start; i < end; i++){
        BamDetail & d = *(*(it + sidx_[i]));
        uint64_t k = (static_cast<uint64_t>(d.barcode) << 32) | d.umi;
        auto res = uhash_.find(k);
        if(res != uhash_.end()){
            d.umi = res->second;
            d.make_hash();
            d.corrected = true;
            corrected++;
            uint32_t mask = (1 << ADNA4::size_) - 1;
            uint8_t * cumi = bam_aux_get(d.b, "UB") + 1;
            //std::cout << "Uncorrected umi = " << reinterpret_cast<char*>(cumi) << "\n";
            uint32_t val = d.umi;
            for(size_t i = 0; i < umi_len_; i++){
                cumi[umi_len_ - i - 1] = ADNA4::alphabet_str_[val & mask];
                val >>= ADNA4::size_;
            }
            //std::cout << "Corrected umi   = " << reinterpret_cast<char*>(cumi) << "\n";
            cor++;
        }
    }
    //std::cout << " umap size = " << uhash_.size() << " tot = " << tot << " cor = " << cor << "\n";

    AlignSummary::bint lbarcode = (*(it + sidx_[start]))->barcode;
    barcodes_.clear();
    bool corrected = false;
    size_t i = start;
    while(i < end){
        BamDetail & d = *(*(it + sidx_[i]));
        if(d.barcode != lbarcode){
            process_barcode_(corrected);
            barcodes_.clear();
            corrected = false;
            lbarcode = d.barcode;
        }
        barcodes_.push_back(&d);
        corrected |= d.corrected;
        i++;
    }
    if(!barcodes_.empty()){
        process_barcode_(corrected);
    }
}

void CollapseWorker::process_barcode_(bool resort){
    if(resort){
        sort(barcodes_.begin(), barcodes_.end(),
            [](const BamDetail * p1, const BamDetail * p2) {
                return std::tie(p1->hash, p1->pos) < std::tie(p2->hash, p2->pos);
            });
    }

    //std::cout << "  Barcode " << barcodes_.front()->barcode << "\n";
    umis_.clear();
    uint32_t lumi = barcodes_.front()->umi;
    for(auto d : barcodes_){
        //std::cout << "    " << d->barcode << " " << d->umi << " " << d->pos << "\n";
        if(d->umi != lumi){
            collapse_umi_();
            umis_.clear();
            lumi = d->umi;
        }
        umis_.push_back(d);
    }
    if(!umis_.empty()) collapse_umi_();
}

std::string debug_read(const BamDetail & d){
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

void CollapseWorker::collapse_umi_(){
    total++;
    rreads += umis_.size();
    build_contigs_();
    fdups_ = 0;
    for(auto & u : umis_){
        if((u->b->core.flag & BAM_FDUP) > 0){
            fdups_++;
        }
    }

    if(build_islands_()){
        //std::cout << "Island group\n";
        fcigar_.clear();
        fquals_.clear();
        fbases_.clear();
        fsplices_.clear();
        fcoverage_.clear();
        freads_ = 0;
        for(size_t i = 0; i < icount_; i++){
            //std::cout << "  Processing group " << i << "\n";
            if(i > 0){
                int d = islands_[i]->lft - islands_[i - 1]->rgt - 1;
                if(d > 0) fcigar_.push_back(CigarElement(d, Cigar::REF_SKIP));
            }
            islands_[i]->merge(umis_, fbases_, fquals_, fcigar_, fsplices_, fcoverage_);
            //bam1_t * bf = umis_[islands_[0]->contigs[0]->index]->b;
            //islands_[i]->debug(genome_[bf->core.tid]);
        }

        freads_ = umis_.size();
        rdups += fdups_;
        rcollapsed += freads_;
        //std::cout << "freads_ = " << freads_ << " umis size = " << umis_.size() << "\n";
        creads++;

        std::sort(fsplices_.begin(), fsplices_.end());
        fsplices_.erase(std::unique(fsplices_.begin(), fsplices_.end()), fsplices_.end());

        if((ccount + 1) >= collapsed.size()){
            cbuffer_.get(collapsed, 100);
        }

        bam1_t * bf = umis_[islands_[0]->contigs[0]->index]->b;
        auto & ref = genome_[bf->core.tid];
        unsigned int lft = islands_[0]->lft;
        unsigned int qlft = 0;
        //std::string rstr, qstr, astr;
        unsigned int NM = 0;
        unsigned int bases = 0;
        for(auto c : fcigar_){
            switch(c.op){
                case Cigar::MATCH:
                    bases += c.len;
                    for(size_t i = 0; i < c.len; i++, lft++, qlft++){
                        NM += (ref.seq[lft] != fbases_[qlft]);
                        //rstr.push_back(ref.seq[lft]);
                        //qstr.push_back(fbases_[qlft]);
                        //astr.push_back(fbases_[qlft] == ref.seq[lft] ? '|' : '*');
                    }
                    break;
                case Cigar::DEL:
                    /*
                    for(size_t i = 0; i < c.len; i++, lft++){
                        rstr.push_back(ref.seq[lft]);
                        qstr.push_back('_');
                        astr.push_back(' ');
                    }
                    */
                    NM += c.len;
                    lft += c.len;
                    break;
                case Cigar::INS:
                    /*
                    for(size_t i = 0; i < c.len; i++, qlft++){
                        rstr.push_back('_');
                        qstr.push_back(fbases_[qlft]);
                        astr.push_back(' ');
                    }
                    */
                    bases += c.len;
                    NM += c.len;
                    qlft += c.len;
                    break;
                case Cigar::REF_SKIP:
                    /*
                    rstr += "   ";
                    qstr += "<->";
                    astr += "   ";
                    */
                    lft += c.len;
                    break;
                default:
                    break;
            }
        }

        ldata.push_back({freads_, bases});

        /*
        //std::cout << "  Fin:  " << fbases_ << "\n";
        std::cout << "  Query: " << qstr << "\n";
        std::cout << "  Align: " << astr << "\n";
        std::cout << "  Ref:   " << rstr << "\n";
        std::cout << "  Qual: ";
        for(size_t i = 0; i < fquals_.size(); i++){
            std::cout << static_cast<char>(fquals_[i] + 33);
        }
        std::cout <<  "\n";
        std::cout << "  Cig : " << fcigar_ << " Reads = " << freads_ << " dups = " << fdups_ << "\n";
        std::cout << "  Qual size = " << fquals_.size() << " seq size = " << fbases_.size() << " cov size = " << fcoverage_.size() << "\n";
        std::cout << "  Cov : ";
        for(auto c : fcoverage_) std::cout << c << ", ";
        std::cout << "\n\n";
        */


        bam1_t * out = collapsed[ccount++];
        make_read_(out, NM);
    }else{
        rdups += fdups_;
        ambig++;
        rlost += umis_.size();
    }
}

void CollapseWorker::make_read_(bam1_t * bam, uint32_t NM) {

    bam1_t * bf = umis_[islands_[0]->contigs[0]->index]->b;

    bam->core.flag = 0;
    if(bam_is_rev(bf)){
        bam->core.flag |= BAM_FREVERSE;
    }

    bam->core.qual = 255; //a.mapq;  All "uniquely" mapped reads
    bam->core.n_cigar = fcigar_.size();
    bam->core.tid = bf->core.tid;
    bam->core.mtid = -1;
    bam->core.mpos = -1;
    bam->core.isize = -1;
    bam->core.pos = islands_[0]->lft;

    std::stringstream ss;

    fqname_.clear();
    fqname_ += bam_aux2Z(bam_aux_get(bf, "CB"));
    fqname_ += '_';
    fqname_ += std::to_string(umis_[islands_[0]->contigs[0]->index]->gid);
    fqname_ += '_';
    fqname_ += bam_aux2Z(bam_aux_get(bf, "UB"));

    auto qlen = fqname_.size() + 1;
    bam->core.l_extranul = (4 - (qlen & 3)) & 3;
    bam->core.l_qname = qlen + bam->core.l_extranul; // Because of the null terminator
    bam->core.l_qseq  = fbases_.size();

    uint32_t dl = bam->core.n_cigar*4 + bam->core.l_qname + bam->core.l_qseq + (bam->core.l_qseq + 1)/2;

    if(bam->m_data < dl) {
        bam->m_data = dl;
        kroundup32(bam->m_data);
        bam->data = (uint8_t*)realloc(bam->data, bam->m_data);
    }

    bam->l_data = dl;
    // Copy the qname
    uint8_t * p = bam->data;
    memcpy(p, fqname_.c_str(), qlen);
    for(size_t i = 0; i < bam->core.l_extranul; i++){
        p[qlen + i] = 0;
    }
    p += bam->core.l_qname;

    for(auto c : fcigar_){
        uint32_t v = c.packed();
        memcpy(p, &v, 4);
        p+=4;
    }

    for(size_t i = 0; i < fbases_.size(); i++){
        uint8_t * pp = p + (i / 2);
        uint8_t mask = 0xFU << (((~i) & 1) << 2);
        *pp = ((*pp) & ~mask) | seq_nt16_table[static_cast<size_t>(fbases_[i])] << (((~i) & 1) << 2);
    }

    p += (bam->core.l_qseq + 1)/2;
    for(size_t i = 0; i < fquals_.size(); i++){
        *p++ = fquals_[i];
    }

    bam->core.bin = hts_reg2bin(bam->core.pos, bam_endpos(bam), 14, 5);
    if((p - bam->data) != dl){
        std::cout << "Data length wrong dl = " << dl << " compared to " << (p - bam->data) << "\n";
        exit(1);
    }

    /*
    auto ptr = bam_aux_get(bf, "XG");
    if(ptr == NULL) {
        std::cout << "Missing XG tag for read " << bf->data << " freads_ = " << freads_ << "\n";
        //print_tags(bf);
    }
    */
    char * x = bam_aux2Z(bam_aux_get(bf, "XG"));

    bam_aux_append(bam, "XG", 'Z', strlen(x) + 1, reinterpret_cast<uint8_t*>(x));

    x = bam_aux2Z(bam_aux_get(bf, "CB"));
    bam_aux_append(bam, "CB", 'Z', strlen(x) + 1, reinterpret_cast<uint8_t*>(x));

    x = bam_aux2Z(bam_aux_get(bf, "UB"));
    bam_aux_append(bam, "UB", 'Z', strlen(x) + 1, reinterpret_cast<uint8_t*>(x));

    bam_aux_append(bam, "NM", 'i', 4, reinterpret_cast<uint8_t*>(&NM));
    bam_aux_append(bam, "ND", 'i', 4, reinterpret_cast<uint8_t*>(&fdups_));
    bam_aux_append(bam, "NR", 'i', 4, reinterpret_cast<uint8_t*>(&freads_));

    char tmpc = bam_aux2A(bam_aux_get(bf, "RE"));
    bam_aux_append(bam, "RE", 'A', 1, reinterpret_cast<uint8_t*>(&tmpc));

    tmpc = bam_aux2A(bam_aux_get(bf, "XS"));
    bam_aux_append(bam, "XS", 'A', 1, reinterpret_cast<uint8_t*>(const_cast<char*>(&tmpc)));

    // The base coverage
    bcov_.clear();
    // 1 for the array type, 4 for the size, and 4 * fcoverage size for 32 bit ints
    bcov_.resize(5 + fcoverage_.size() * 4);
    bcov_[0] = 'I';
    uint8_t * bp = &bcov_[0];
    (*bp++) = 'I';
    uint32_t l = fcoverage_.size();
    memcpy(bp, &l, sizeof(uint32_t));
    bp+=4;
    memcpy(bp, &fcoverage_[0], fcoverage_.size() * 4);
    bam_aux_append(bam, "XC", 'B', bcov_.size(), reinterpret_cast<uint8_t*>(&bcov_[0]));


    fqname_.clear();
    for(size_t i = 0; i < fsplices_.size(); i++){
        if(i > 0) fqname_ += '|';
        fqname_ += std::to_string(fsplices_[i].first);
        fqname_ += ',';
        fqname_ += std::to_string(fsplices_[i].second);
    }

    // Splice junctions
    bam_aux_append(bam, "XJ", 'Z', fqname_.size() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(fqname_.c_str())));
}

void CollapseWorker::build_contigs_() {
    ccount_ = 0;
    for(size_t i = 0; i < umis_.size(); i++){
        if(i > 0 && umis_[i - 1]->pos == umis_[i]->pos){
            umis_[i]->b->core.flag |= BAM_FDUP;
        }

        //unsigned int fs = ccount_;
        ccount_++;
        if(ccount_ >= contigs_.size()) {
            for(size_t i = 0; i < 10; i++) contigs_.push_back(new ReadContig());
        }

        ReadContig * f = contigs_[ccount_-1];
        f->cig.clear();
        f->ins.clear();
        f->lft = umis_[i]->b->core.pos;
        f->rgt = umis_[i]->b->core.pos;
        f->qlft = 0;
        f->index = i;
        f->next = nullptr;
        f->prev = nullptr;

        uint32_t * cig = bam_get_cigar(umis_[i]->b);
        unsigned int cs = 0;
        if(bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP){
            f->qlft += bam_cigar_oplen(cig[0]);
            cs++;
        }
        f->qrgt = f->qlft;
        auto s = bam_get_seq(umis_[i]->b);
        auto q = bam_get_qual(umis_[i]->b);
        for(size_t j = cs; j < umis_[i]->b->core.n_cigar; j++){
            int op =  bam_cigar_op(cig[j]);
            int len = bam_cigar_oplen(cig[j]);
            if(op == BAM_CREF_SKIP){
                // Want one based
                ccount_++;
                if(ccount_ >= contigs_.size()) {
                    for(size_t i = 0; i < 10; i++) contigs_.push_back(new ReadContig());
                }

                auto fn = contigs_[ccount_ - 1];
                fn->cig.clear();
                fn->ins.clear();
                fn->index = i;
                fn->lft = f->rgt + len;
                fn->rgt = fn->lft;
                fn->qlft = f->qrgt;
                fn->qrgt = f->qrgt;
                fn->next = nullptr;
                fn->prev = f;
                f->next = fn;
                f = fn;
            }else if(op == BAM_CMATCH){
                f->qrgt += len;
                f->rgt  += len;
                f->cig.push_back(CigarElement(cig[j]));
            }else if(op == BAM_CDEL){
                f->rgt  += len;
                f->cig.push_back(CigarElement(cig[j]));
            }else if(op == BAM_CINS){
                f->ins.push_back({f->rgt});
                for(int x = 0; x < len; x++, f->qrgt++){
                    f->ins.back().bases.push_back("=ACMGRSVTWYHKDBN"[bam_seqi(s, f->qrgt)]);
                    f->ins.back().quals.push_back(q[f->qrgt]);
                }
                f->cig.push_back(CigarElement(cig[j]));
            }
        }
    }

    // Need to subtract one from each rgt
    for(size_t i = 0; i < ccount_; i++){
        auto & f = *contigs_[i];
        f.rgt--;
        f.qrgt--;
    }

    sort(contigs_.begin(), contigs_.begin() + ccount_,
        [](const ReadContig * p1, const ReadContig * p2) {
            return p1->lft < p2->lft;
    });
}

bool CollapseWorker::build_islands_(){
    icount_ = 1;
    if(icount_ >= islands_.size()) {
        for(size_t i = 0; i < 10; i++) islands_.push_back(new ReadIsland());
    }

    ReadIsland * isl = islands_[icount_ - 1];
    isl->contigs.clear();
    isl->splices.clear();
    isl->starts.clear();
    isl->ends.clear();
    isl->conns.clear();
    isl->lft = contigs_.front()->lft;
    isl->rgt = contigs_.front()->rgt;
    for(size_t i = 0; i < ccount_; i++){
        auto & f = *contigs_[i];
        //std::cout << "  Check " << isl->lft << " - " << isl->rgt << " vs " << f.lft << " - " << f.rgt
        //    << " c1 = " << (isl->lft <= f.rgt) << " c2 = " << (f.lft <= isl->rgt)
        //    << "\n";
        if(isl->lft <= f.rgt && f.lft <= isl->rgt){
            isl->rgt = std::max(isl->rgt, f.rgt);
            isl->contigs.push_back(contigs_[i]);
            contigs_[i]->island = icount_ - 1;
            if(f.prev != nullptr){
                islands_[f.prev->island]->conns.push_back(f.island);
                islands_[f.prev->island]->splices.push_back({f.prev->rgt, f.lft});
                islands_[f.prev->island]->ends.push_back(f.prev->rgt);
                islands_[f.island]->starts.push_back(f.lft);
            }
        }else{
            icount_++;
            if(icount_ >= islands_.size()) {
                for(size_t j = 0; j < 10; j++) islands_.push_back(new ReadIsland());
            }
            isl = islands_[icount_ - 1];
            isl->contigs.clear();
            isl->conns.clear();
            isl->splices.clear();
            isl->starts.clear();
            isl->ends.clear();
            islands_[icount_ - 2]->conns.push_back(icount_ - 1);
            isl->contigs.push_back(contigs_[i]);
            contigs_[i]->island = icount_ - 1;
            isl->lft = f.lft;
            isl->rgt = f.rgt;
            if(f.prev != nullptr){
                //isl->CollapseSplices.push_back({f.prev->lft, f.rgt});
                islands_[f.prev->island]->splices.push_back({f.prev->rgt, f.lft});
                islands_[f.prev->island]->conns.push_back(f.island);
                islands_[f.prev->island]->ends.push_back(f.prev->rgt);
                islands_[f.island]->starts.push_back(f.lft);
            }
        }
    }

    for(auto i : islands_){
        std::sort(i->conns.begin(), i->conns.end());
        i->conns.erase(std::unique(i->conns.begin(), i->conns.end()), i->conns.end());
        std::sort(i->starts.begin(), i->starts.end());
        i->starts.erase(std::unique(i->starts.begin(), i->starts.end()), i->starts.end());
        std::sort(i->ends.begin(), i->ends.end());
        i->ends.erase(std::unique(i->ends.begin(), i->ends.end()), i->ends.end());
    }
    /*
    for(size_t i = 0; i < umis_.size(); i++){
        std::cout << "  Read: "  << i << " " << debug_read(*umis_[i]) << "\n";
    }
    */
    bool has_skip = false;
    for(size_t i = 0; i < icount_; i++){
        auto & isr = *islands_[i];
        //std::cout << "  Island " << i << " lft " << isr.lft << " - " << isr.rgt << " conns: [";
        bool skip = false;
        if(!isr.conns.empty()){
            //std::cout << isr.conns.front();
            skip |= isr.conns.front() != (i + 1);
            for(size_t j = 1; j < isr.conns.size(); j++){
                //std::cout << ", " << isr.conns[j];
                skip |= isr.conns[i] != (i + 1);
            }
        }

        unsigned int max_sdist = 0, max_edist = 0;

        for(size_t j = 0; j < isr.starts.size(); j++) max_sdist = std::max(isr.starts[j] - isr.lft, max_sdist);
        for(size_t j = 0; j < isr.ends.size(); j++) max_sdist = std::max(isr.rgt - isr.ends[j], max_sdist);
        if(max_sdist > 5 || max_edist > 5){
            skip = true;
        }
        has_skip |= skip;
    }

    return !has_skip;
}
