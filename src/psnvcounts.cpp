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
#include "psnvcounts.hpp"
#include "collapse_hist.hpp"
#include <unordered_map>
#include "parallel-hashmap/parallel_hashmap/phmap.h"
#include <zlib.h>
#include <fstream>
#include "bam_genes_aux.hpp"
#include "bam_genes.hpp"
#include "tokenizer.hpp"
#include "sbam_merge.hpp"
#include "gzstream.hpp"
#include "htslib/htslib/hts_endian.h"

using namespace gwsc;

/*
struct CTagOut{
    using bint = uint32_t;
    CTagOut() {

    }

    CTagOut(bint barcode, bool intronic, uint32_t umi, uint32_t gene_id, uint32_t snv_id)
        : barcode(barcode), intronic(intronic), umi(umi), gene_id(gene_id), snv_id(snv_id)
    {

    }

    bool operator<(const CTagOut & rhs) const {
        return std::tie(barcode, umi, gene_id, snv_id) < std::tie(rhs.barcode, rhs.umi, rhs.gene_id, snv_id);
    }

    bint     barcode:31;
    bint     intronic:1;
    uint32_t umi;
    uint32_t gene_id;
    uint32_t snv_id;
    //uint32_t bases;
};
*/

argagg::parser ProgSNVCounts::parser() const {
    argagg::parser argparser {{
        { "txidx", {"-i", "--index"},
          "Transcript Index", 1},
        { "barcodes", {"-b", "--barcodes"},
          "Barcode count prefix", 1},
        { "output", {"-o", "--output"},
          "Output prefix", 1},
        { "snvs", {"-s", "--snvs"},
          "Tab separated list of strand specific SNVs must have the following columns chrom, pos, ref, alt, strand", 1},
        { "cellranger", {"-c", "--cellranger"},
          "Indicates the merged bam file is from cell ranger", 0},
        { "tags", {"-t", "--tags"},
          "If this bam file is collapsed write tags that can be used to look for poorly mapping barcode UMI gene combinations", 0},
        { "library", {"-l", "--library"},
          "libary type (V2)", 1},
        { "help", {"-h", "--help"},
          "shows this help message", 1},
      }};
    return argparser;
}


void ProgSNVCounts::load() {
    iprefix_ = args_["txidx"].as<std::string>();
    isnvs_ = args_["snvs"].as<std::string>();
    bcin_ = args_["barcodes"].as<std::string>();
    cellranger_ = args_["cellranger"];
    //tags_ = args_["tags"];
    lib_ = args_["library"].as<std::string>("V2");
    outp_ = args_["output"].as<std::string>();
    if(args_.pos.size() != 1){
        throw std::runtime_error("Missing the output option");
    }
    bamin_ = args_.as<std::string>(0);
}


void ProgSNVCounts::parse_snvs_(){
    std::unordered_map<std::string, int> tidmap;
    for(auto & r : idx_.refs()){
        tidmap[r.name] = r.tid;
    }
    snvs_.resize(tidmap.size());

    map_idx_ = 0;

    FileWrapper fb(isnvs_);
    {
        ParserTokens toks;
        std::vector<unsigned int> tmp;
        unsigned int chromi, posi, refi, alti, strandi;
        chromi = posi = refi = alti = strandi = std::numeric_limits<unsigned int>::max();
        fb.tokenize_line(toks);
        for(size_t i = 0; i < toks.size(); i++){
            if(toks[i] == "chrom")       chromi = i;
            else if(toks[i] == "pos")    posi = i;
            else if(toks[i] == "ref")    refi = i;
            else if(toks[i] == "alt")    alti = i;
            else if(toks[i] == "strand") strandi = i;
        }

        bool failed = false;
        if(chromi == std::numeric_limits<unsigned int>::max()){
            failed = true;
            std::cout << "Could not find a chrom column\n";
        }
        if(posi == std::numeric_limits<unsigned int>::max()){
            failed = true;
            std::cout << "Could not find a pos column\n";
        }
        if(refi == std::numeric_limits<unsigned int>::max()){
            failed = true;
            std::cout << "Could not find a ref column\n";
        }
        if(alti == std::numeric_limits<unsigned int>::max()){
            failed = true;
            std::cout << "Could not find a alt column\n";
        }
        if(strandi == std::numeric_limits<unsigned int>::max()){
            failed = true;
            std::cout << "Could not find a strand column\n";
        }
        if(failed) exit(EXIT_FAILURE);
        while(fb.tokenize_line(toks) >= 0){
            if(toks.empty() || toks[chromi].empty()) continue;
            int tid = tidmap[toks[chromi]];
            unsigned int pos = std::stoi(toks[posi]);
            char ref = toks[refi][0];
            char alt = toks[alti][0];
            char strand = toks[strandi][0];
            if((size_t)tid >= snvs_.size()){
                snvs_.resize(tid + 1);
            }
            snvs_[tid].push_back(SNV(tid, pos, ref, alt, strand, snv_count_));
            snv_count_++;
        }
        size_t tid = 0;
        for(auto & snvs : snvs_){
            if(snvs.empty()) continue;
            //std::cout << " " << tid << " " << idx_.ref(tid).name << " snvs = " << snvs.size() << "\n";
            std::sort(snvs.begin(), snvs.end());
            tid++;
        }
    }
    map_idx_ = snv_count_;
}

void ProgSNVCounts::read_passed_(unsigned int blength){
    FileWrapper in(bcin_);
    std::string line;
    in.get_line(line);
    size_t index = 0;
    while(in.get_line(line) > -1){
        barcodes_.push_back(line.substr(0, blength));
        bchash_[barcodes_.back()] = index++;
    }
    tout << "Read " << barcodes_.size() << " passed barcodes\n";
}


unsigned int ProgSNVCounts::count_dups_(std::string & xr) {
    positions_.clear();
    //positions_.push_back({lft, rgt});
    Tokenizer::get(xr, ';', ttoks_);
    for(auto & t : ttoks_){
        Tokenizer::get(t, ',', ptoks_);
        unsigned int lft = std::stoi(ptoks_[0]);
        unsigned int rgt = std::stoi(ptoks_[1]);
        positions_.push_back({lft, rgt});
    }
    std::sort(positions_.begin(), positions_.end());
    size_t osize = positions_.size();
    positions_.erase(std::unique(positions_.begin(), positions_.end()), positions_.end());
    return osize - positions_.size();
}

void print_aux(bam1_t * b){
    auto s = bam_get_aux(b); // aux
    auto end = b->data + b->l_data;
    while (end - s >= 4) {
        uint8_t type, key[2];
        key[0] = s[0]; key[1] = s[1];
        s += 2; type = *s++;
        std::cout << "  " << (char)key[0] << (char)key[1] << " x = " << (s) << " :";
        if (type == 'A') {
            std::cout << "A:" << *s;
            ++s;
        } else if (type == 'C') {
            std::cout << "i:" << *s;
            ++s;
        } else if (type == 'c') {
            std::cout << "i:" << *(int8_t*)s;
            ++s;
        } else if (type == 'S') {
            if (end - s >= 2) {
                std::cout << "i:" << le_to_u16(s);
                s += 2;
            } else { std::cout << "  Aux problem end = " << end << " s = " << (unsigned long int)s << " diff = " << (end - s) << "\n"; return; };
        } else if (type == 's') {
            if (end - s >= 2) {
                std::cout << "i:" << le_to_i16(s);
                s += 2;
            } else { std::cout << "  Aux problem end = " << end << " s = " << (unsigned long int)s << " diff = " << (end - s) << "\n"; return; };
        } else if (type == 'I') {
            if (end - s >= 4) {
                std::cout << "i:" << le_to_u32(s);
                s += 4;
            } else { std::cout << "  Aux problem end = " << end << " s = " << (unsigned long int)s << " diff = " << (end - s) << "\n"; return; };
        } else if (type == 'i') {
            if (end - s >= 4) {
                std::cout << "i:" << le_to_i32(s);
                s += 4;
            } else { std::cout << "  Aux problem end = " << end << " s = " << (unsigned long int)s << " diff = " << (end - s) << "\n"; return; };
        } else if (type == 'f') {
            if (end - s >= 4) {
                std::cout << "f:" << le_to_float(s);
                s += 4;
            } else { std::cout << "  Aux problem end = " << end << " s = " << (unsigned long int)s << " diff = " << (end - s) << "\n"; return; };

        } else if (type == 'd') {
            if (end - s >= 8) {
                std::cout << "d:" << le_to_double(s);
                s += 8;
            } else { std::cout << "  Aux problem end = " << end << " s = " << (unsigned long int)s << " diff = " << (end - s) << "\n"; return; };
        } else if (type == 'Z' || type == 'H') {
            std::cout << type << ":";
            while (s < end && *s) std::cout << *s++;
            if (s >= end)
                { std::cout << "  Aux problem end = " << end << " s = " << (unsigned long int)s << " diff = " << (end - s) << "\n"; return; };
            ++s;
        } else { // Unknown type
            std::cout << "Unknown type\n";
            return;
        }
        std::cout << "\n";
    }
}

template <typename T, typename B>
inline int ProgSNVCounts::run_()  {
    tout << "Loading the index\n";
    idx_.load(iprefix_);
    read_passed_(B::BARCODE_LEN);
    tout << "Parsing the SNVs\n";
    parse_snvs_();

    ReadBaseHist hists(bchash_.size());

    tout << "Processing " << bamin_ << "\n";
    T proc(idx_, B::LibraryStrand);
    BamReader bin;
    BamDetail read;
    bin.set_bam(bamin_);
    int ltid = -1;
    std::vector<SNV>::iterator start, end;
    size_t total = 0, written = 0;
    std::string btmp;
    std::string umi_tmp;

    //std::vector<CTagOut> ctags;
    gzofstream zout(outp_ + "_reads.txt.gz");
    zout << "snv_idx\tbarcode\tbases\treads\tstrand\n";
    while(bin.next(read.b) != nullptr){
        auto bam = read.b;
        if(!proc(read, 0)) continue;

        if(ltid != bam->core.tid){
            ltid = bam->core.tid;
            start = snvs_[ltid].begin();
            end = snvs_[ltid].end();
        }

        //std::cout << bam->data << " tags: ";
        //print_aux(bam);
        //std::cout << "\n";
        int pos = bam->core.pos;
        int qpos = 0;
        while(start != end && start->pos < pos){
            start++;
        }
        auto it = start;
        uint32_t * cig = bam_get_cigar(bam);
        char strand = bam_aux2A(bam_aux_get(bam, "XS"));
        uint32_t bases = 0;
        auto seq = bam_get_seq(bam);
        overlaps_.clear();
        refalt_.clear();

        total++;
        if(total % 10000000 == 0){
            tout << "Processed current = " << idx_.ref(ltid).name << " " << pos << " " << total << " used = " << written << " map size: " << snvmap_.size() 
                << " SNV pair size = " << snvpairs_.size() << "\n"; 
        }

        btmp = bam_aux2Z(bam_aux_get(bam, "CB"));

        size_t dash = btmp.find('-');
        if(dash != std::string::npos){
            btmp = btmp.erase(dash);
        }
        auto ptr = bam_aux_get(bam, "NR");
        int NR = 1;
        if(ptr != NULL){
            NR = bam_aux2i(ptr);
        }       
        auto bit = bchash_.find(btmp);
        if(bit == bchash_.end()) {
            continue;
        }
        auto barcode = bit->second;
        for(size_t j = 0; j < bam->core.n_cigar; j++){
            CigarElement elem(cig[j]);
            auto op = elem.op;
            auto len = elem.len;
            switch(op){
                case Cigar::INS:
                    qpos += len;
                    break;
                case Cigar::DEL:
                    pos += len;
                    break;
                case Cigar::SOFT_CLIP:
                    qpos += len;
                    break;
                case Cigar::MATCH:
                    while(it != end && it->pos < pos) it++;
                    if(it != end){
                        for(size_t i = 0; i < len; i++){
                            if(it != end){
                                if(it->pos < pos) it++;
                                if(pos == it->pos && strand == it->strand){
                                    char base = "NACNGNNNTNNNNNNN"[bam_seqi(seq, qpos)];
                                    if(base == it->alt){
                                        overlaps_.push_back(it->index);
                                        refalt_.push_back(RefAlt(it->index, true));
                                    }else if(base == it->ref){
                                        overlaps_.push_back(it->index);
                                        refalt_.push_back(RefAlt(it->index, false));
                                    }
                                }
                            }
                            pos++;
                            qpos++;
                        }
                    }
                    bases += len;
                    break;
                case Cigar::REF_SKIP:
                    pos += len;
                    break;
                default:
                    break;
            }
        }
        ptr = bam_aux_get(bam, "XR");
        unsigned int dups = 0;
        if(ptr != NULL){
            btmp = bam_aux2Z(ptr);
            dups = count_dups_(btmp);
        }
        hists.add_count(barcode, NR, dups, bases);
        uint32_t snvidx = std::numeric_limits<uint32_t>::max();
        if(!overlaps_.empty()){
            if(overlaps_.size() == 1){
                snvidx = overlaps_.front();
            }else{
                std::sort(overlaps_.begin(), overlaps_.end());
                auto it = snvmap_.insert({overlaps_, map_idx_});
                if(it.second) map_idx_++;
                snvidx = it.first->second;
            }
            zout << snvidx << "\t" << barcodes_[barcode] << "\t" << bases << "\t" << NR << "\t" << strand << "\n";
            written++;
        }

        if(!refalt_.empty()){
            for(size_t i = 0; i < refalt_.size() - 1; i++){
                for(size_t j  = i + 1; j < refalt_.size(); j++){
                    uint64_t key = (static_cast<uint64_t>(refalt_[i].snv_idx) << 32) | refalt_[j].snv_idx;
                    auto it = snvpairs_.insert(std::make_pair(key, std::array<uint32_t, 4>{}));
                    uint32_t ckey = (refalt_[i].alt << 1) | refalt_[j].alt;
                    it.first->second[ckey]++;
                }
            }
        }



        /*
        if(tags_){
            char A = bam_aux2A(bam_aux_get(bam, "RA"));
            ptr = bam_aux_get(bam, "UB");
            umi_tmp = bam_aux2Z(ptr);
            uint32_t umi = seq2int<gwsc::ADNA4, uint32_t>(umi_tmp);
            ctags.push_back(CTagOut(barcode, A == 'N', umi, read.gid, snvidx));
        }
        */
    }
    tout << "Writing the SNV occurence map\n";
    write_map_(outp_ + "_map.txt.gz");
    hists.write(outp_ + "_", barcodes_);

    std::vector<EdgeOut> out;
    std::vector<SNV*> smap;
    for(auto & v : snvs_){
        for(auto & s : v){
            if(s.index >= smap.size()){
                smap.resize(s.index + 1);
            }
            smap[s.index] = &s;
        }
    }

    std::cout << "Building SNV Edges\n";
    for(auto & k : snvpairs_){
        uint32_t i1 = (k.first >> 32);
        uint32_t i2 = k.first & 0xFFFFFFFF;
        if(i1 >= smap.size() || i2 >= smap.size()){
            std::cout << "i1 = " << i1 << " i2 = "<< i2 << " k = " << k.first << " size = " << smap.size() << "\n";
        }
        auto & s1 = *smap[i1];
        auto & s2 = *smap[i2];
        out.push_back(EdgeOut(s1.index, s2.index, s1.tid, s1.pos, s2.pos, s1.ref, s1.alt, s2.ref, s2.alt, s1.strand, k.second[0], k.second[3], k.second[2], k.second[1]));
    }

    std::cout << "Built " << out.size() << " edges\n";
    std::sort(out.begin(), out.end());

    {
        gzofstream zeout(outp_ + "_edges.txt.gz");
        zeout << "snv_idx_1\tsnv_idx_2\tchrom\tpos_1\tpos_2\tref_1\talt_1\tref_2\talt_2\tstrand\tRR\tAA\tRA\tAR\n";
        for(auto & o : out){
            zeout << o.i1 << "\t" << o.i2 << "\t" << idx_.refs()[o.tid].name << "\t" << o.pos1 << "\t" << o.pos2 << "\t" << o.ref1 << "\t" << o.alt1 << "\t" 
                << o.ref2 << "\t" << o.alt2 << "\t" << o.strand << "\t" << o.rr << "\t" << o.aa << "\t" << o.ra << "\t" << o.ar << "\n";
        }
    }

    /*
    if(tags_){
        tout << "Sorting and writing the read tags\n";
        std::string outf = outp_ + "tags.gz";
        std::sort(ctags.begin(), ctags.end());
        gzFile zp = gzopen(outf.c_str(), "w");
        for(auto & tag : ctags){
            gzwrite(zp, reinterpret_cast<char*>(&tag), sizeof(CTagOut));
        }
        tout << "Wrote " << ctags.size() << " tags\n";
    }
    */
    
    tout << "Done total = " << total << " used = " << written << "  map size: " << snvmap_.size() << "\n"; 
    return EXIT_SUCCESS;
}

void ProgSNVCounts::write_map_(const std::string & out){
    gzofstream zo(out);
    zo << "snv_id\tsnvs\n";
    std::vector< const SNVKey * > tout;
    tout.resize(map_idx_ - snv_count_, nullptr);
    for(auto const & m : snvmap_){
        if(m.second >= snv_count_)
            tout[m.second - snv_count_] = &m.first;
    }
    for(size_t i = 0; i < tout.size(); i++){
        zo << (i + snv_count_) << "\t" << tout[i]->front();
        for(size_t j = 1; j < tout[i]->size(); j++){
            zo << "," << tout[i]->at(j);
        }
        zo << "\n";
    }
}

int ProgSNVCounts::run(){
    if(cellranger_){
        if(lib_ == "V2"){
            return run_<BamCellRangerProcessor, Reader10X_V2>();
        }else if(lib_ == "V3"){
            return run_<BamCellRangerProcessor, Reader10X_V3>();
        }else if(lib_ == "V2_5P"){
            return run_<BamCellRangerProcessor, Reader10X_V2_5P>();
        }else if(lib_ == "V3_5P"){
            return run_<BamCellRangerProcessor, Reader10X_V3_5P>();
        }
    }

    if(lib_ == "V2"){
        return run_<BamScSNVProcessor, Reader10X_V2>();
    }else if(lib_ == "V3"){
        return run_<BamScSNVProcessor, Reader10X_V3>();
    }else if(lib_ == "V2_5P"){
        return run_<BamScSNVProcessor, Reader10X_V2_5P>();
    }else if(lib_ == "V3_5P"){
        return run_<BamScSNVProcessor, Reader10X_V3_5P>();
    }
    return EXIT_FAILURE;
}


