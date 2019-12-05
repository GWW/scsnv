#include "collapse_aux.hpp"
#include "../util/sbam_writer.hpp"
#include <algorithm>
#include <numeric>

using namespace gwsc;

void ReadIsland::merge(std::vector<BamDetail *> & umis, std::string & fbases, std::string & fquals, 
        CigarString & fcigar, std::vector<CollapseSplice> & fsplices, std::vector<unsigned int> & coverage){

    N = contigs.size();
    merge_ins_();

    for(auto f : contigs){
        if(f->next){
            splices.push_back({f->rgt, f->next->lft});
        }
    }

    std::sort(splices.begin(), splices.end());
    splices.erase(std::unique(splices.begin(), splices.end()), splices.end());
    fsplices.insert(fsplices.end(), splices.begin(), splices.end());

    L = rgt - lft + 1;

    merge_bases_(umis);

    size_t istart = 0;
    for(size_t i = 0; i < L; i++){
        if(istart < ins.size() && i == ins[istart].lft){
            unsigned int cov = std::accumulate(icounts.begin() + istart * N, icounts.begin() + (istart + 1) * N, 0);
            for(size_t x = 0; x < ins[istart].len(); x++){
                fquals.push_back(ins[istart].quals[x]);
                fbases.push_back(ins[istart].bases[x]);
                coverage.push_back(cov);
            }
            fcigar.push_back({ins[istart].len(), Cigar::INS});
        }

        if(cbases[i] == '_'){
            if(fcigar.empty() || fcigar.back().op != Cigar::DEL){
                fcigar.push_back(CigarElement(1, Cigar::DEL));
            }else{
                fcigar.back().len++;
            }
            coverage.push_back(bcounts[i * 6 + 5]);
        }else{
            fbases.push_back(cbases[i]);
            fquals.push_back(cquals[i]);
            int bc = ADNA5::ltable_[static_cast<unsigned int>(cbases[i])];
                coverage.push_back(bcounts[i * 6 + bc]);

            if(fcigar.empty() || fcigar.back().op != Cigar::MATCH){
                fcigar.push_back(CigarElement(1, Cigar::MATCH));
            }else{
                fcigar.back().len++;
            }
        }
    }

    //debug();
}

void ReadIsland::debug(){
    std::cout << "  Island " << " lft " << lft << " - " << rgt << "\n";
    std::cout << "    Contigs\n";
    for(auto f : contigs){
        std::cout << "      " << *f << "\n";
    }

    if(!conns.empty()){
        std::cout << "    Connections = ";
        for(size_t j = 0; j < conns.size(); j++){
            std::cout << " " << conns[j];
        }
        std::cout << "\n";
    }

    if(!splices.empty()){
        std::cout << "    Splices = ";
        for(size_t j = 0; j < splices.size(); j++){
            std::cout << " " << splices[j].first << " - " << splices[j].second;
        }
        std::cout << "\n";
    }

    if(!ins.empty()){
        std::cout << "    Insertions = ";
        for(size_t j = 0; j < ins.size(); j++){
            std::cout << " " << ins[j].lft << " L = " << ins[j].len() 
                << " N = " << std::accumulate(icounts.begin() + j * N, icounts.begin() + (j + 1) * N, 0)
                << " C = " << icoverage[j];
        }
        std::cout << "\n";
    }

    std::cout << "          ";
    for(size_t j = 0; j < L; j++){
        if(j % 10 == 0){
            std::cout << ((j / 10) % 10);
        }else{
            std::cout << " ";
        }
    }
    std::cout << "\n";
    std::cout << "    Cons: " << cbases << "\n";
    for(size_t i = 0; i < N; i++){
        std::cout << "          " << bases.substr(i * L, L) << "\n";
    }
    std::cout << "\n";

}

void ReadIsland::merge_bases_(std::vector<BamDetail *> & umis){
    bcounts.clear();
    cquals.clear();
    cbases.clear();
    qmax.clear();
    bcounts.resize(L * 6, 0);
    cquals.resize(L, 0);
    cbases.resize(L, 0);
    qmax.resize(L * 6, 0);
    bases.clear();
    bases.resize(L * N, ' ');

    quals.clear();
    quals.resize(L * N);

    unsigned int ilft = lft;

    for(size_t j = 0; j < N; j++){
        auto & f = *contigs[j];
        unsigned int qlft = f.qlft;
        unsigned int cs = (j * L) + (f.lft - ilft);
        auto b = umis[f.index]->b;
        auto s = bam_get_seq(b);
        auto q = bam_get_qual(b);
        for(auto c : f.cig){
            if(c.op == Cigar::MATCH){
                for(size_t i = 0; i < c.len; i++, qlft++, cs++){
                    bases[cs] = "=ACMGRSVTWYHKDBN"[bam_seqi(s, qlft)];
                    quals[cs] = q[qlft];
                }
            }else if(c.op == Cigar::DEL){
                for(size_t i = 0; i < c.len; i++, cs++){
                    bases[cs] = '_';
                }
            }else if(c.op == Cigar::INS){
                qlft += c.len;
            }
        }
    }

    for(size_t i = 0; i < N; i++){
        unsigned int cs = (i * L);
        for(size_t j = 0; j < L; j++, cs++){
            auto b = bases[cs];
            if(b != ' '){
                int bc = 5;
                if(b != '_'){
                    bc = ADNA5::ltable_[static_cast<unsigned int>(bases[cs])];
                    if(bc == -1) bc = 4;
                }
                bcounts[j * 6 + bc]++;
                qmax[j * 6 + bc] = std::max(quals[cs], qmax[j * 6 + bc]);
            }
        }
    }

    for(size_t i = 0; i < L; i++){
        unsigned int tot = 0;
        unsigned int mj = 0;
        unsigned int mv = 0;
        for(size_t j = 0; j < 6; j++){
            auto v = bcounts[i * 6 + j];
            if(v > mv){
                mj = j;
                mv = v;
            }
            tot += v;
        }

        auto c = bcounts[i * 6 + mj];
        double m = 1.0 * c / tot;
        if(m >= 0.6){
            cbases[i] = "ACGTN_"[mj];
            cquals[i] = qmax[i * 6 + mj];
        }else{
            cbases[i] = 'N';
            cquals[i] = 0;
        }

    }

}

void ReadIsland::merge_ins_(){
    ins.clear();
    for(auto f : contigs){
        ins.insert(ins.end(), f->ins.begin(), f->ins.end());
    }

    if(ins.empty()) return;
    std::sort(ins.begin(), ins.end());

    auto start = ins.begin();
    for(auto it = std::next(ins.begin()); it != ins.end(); it++){
        if(*start == *it){
            //Merge the quals
            for(size_t i = 0; i < start->len(); i++){
                start->quals[i] = std::max(start->quals[i], it->quals[i]);
            }
        }else{
            if(++start != it){
                *start = std::move(*it);
            }
        }
    }
    ins.erase(++start, ins.end());

    icounts.clear();
    icounts.resize(ins.size() * N);
    icoverage.clear();
    icoverage.resize(ins.size());

    for(size_t i = 0; i < N; i++){
        auto & f = *contigs[i];
        // Check how often each CollapseInsertion is found
        for(size_t j = 0; j < ins.size(); j++){
            auto & ii = ins[j];
            for(auto & fi : f.ins)
                if(ii == fi) icounts[j * N + i]++;
            if(f.lft <= ii.lft && ii.lft <= f.rgt){
                icoverage[j]++;
            }
        }
    }

    size_t i = 0; 
    for(size_t j = 0; j < ins.size(); j++){
        unsigned int tot = std::accumulate(icounts.begin() + j * N, icounts.begin() + (j + 1) * N, 0);
        double d = 1.0 * tot / icoverage[j];
        if(d >= 0.6){
            if(i != j) {
                std::swap(ins[i], ins[j]);
                std::swap(icoverage[i], icoverage[j]);
                for(size_t x = 0; x < N; x++){
                    std::swap(icounts[j * N + x], icounts[i * N + x]);
                }
            }
            ibases += ins[i].len();
            i++;
        }
    }
    ins.resize(i);
    icounts.resize(i * N);
    icoverage.resize(i);

    for(auto & i : ins){
        i.lft -= lft;
    }
}

CollapsedBamWriter::CollapsedBamWriter(const std::string & out, unsigned int bam_write_threads, const bam_hdr_t * bh){
    bh_ = bam_hdr_dup(bh);
    std::string sf = out + ".bam";
    bam_out_ = sam_open(sf.c_str(), "wb");
    if(bam_write_threads > 1) hts_set_threads(bam_out_, bam_write_threads);
    if(sam_hdr_write(bam_out_, bh_) < 0) {
        std::cerr << "Error writing header\n";
        exit(1);
    }
    sf = out + "_collapsed.bam";
    merged_out_ = sam_open(sf.c_str(), "wb");
    if(bam_write_threads > 1) hts_set_threads(merged_out_, bam_write_threads);
    if(sam_hdr_write(merged_out_, bh_) < 0) {
        std::cerr << "Error writing header\n";
        exit(1);
    }
}

CollapsedBamWriter::~CollapsedBamWriter(){
    sam_close(bam_out_);
    sam_close(merged_out_);
    bam_hdr_destroy(bh_);

    for(auto b : collapsed){
        bam_destroy1(b);
    }
    collapsed.clear();
}

void CollapsedBamWriter::start(BamBuffer * buffer){
    buffer_ = buffer;
    thread_ = std::thread(std::ref(*this));
}

void CollapsedBamWriter::join() {
    thread_.join();
}

void CollapsedBamWriter::operator()(){
    auto ret = buffer_->get_all_unsafe();
    for(auto it = ret.first; it != ret.second; it++){
        if(sam_write1(bam_out_, bh_, (*it)->b) < 0){
            std::cerr << "Error writing bam record\n";
            exit(1);
        }
        wrote++;
    }

    SortBamTidPos psort;
    std::sort(collapsed.begin(), collapsed.end(), psort);

    for(auto c : collapsed){
        if(sam_write1(merged_out_, bh_, c) < 0){
            std::cerr << "Error writing bam record\n";
            exit(1);
        }
        cwrote++;
    }
}
