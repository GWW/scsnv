#pragma once
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

#include "align_aux.hpp"
#include "dfix.hpp"
#include "sequence.hpp"
#include "sparsepp/sparsepp/spp.h"
#include <cassert>
#include <stack>
#include <iomanip>
#include <vector>
#include <unordered_set>

namespace gwsc{

// For each barcode I need to keep a list of UMI's that are associated 
// with each gene so I can look at the repeats later

struct DupNode {
    void reset(){
        umi = 0;
        count = 0;
        //intronic = 0;
        //positions.clear();
    }

    void init(const AlignSummary & a){
        umi = a.umi;
        count = 1;
        //intronic += a.intronic;
        //positions.clear();
        //positions.push_back(a.pos);
    }

    bool merge(const AlignSummary & a){
        if(a.umi != umi) return false;
        count++;
        //intronic += a.intronic;
        //positions.push_back(a.pos);
        return true;
    }

    /*
    void prepare(){
        if(!positions.empty()){
            std::sort(positions.begin(), positions.end());
            auto last = std::unique(positions.begin(), positions.end());
            positions.erase(last, positions.end());
        }
    }
    */

    uint32_t hash() const {
        return umi;
    }

    //std::vector<uint64_t> positions;

    uint32_t umi;
    uint32_t count;
    //uint32_t intronic;
};

template <typename N>
class DupGraph {
    public:
        void reset(){
            for(auto & e :edges_) e.clear();
            for(auto & n :nodes_) n.reset();
            N_ = 0;
            E_ = 0;
        }

        void add_edge(size_t i, size_t j){
            edges_[i].push_back(j);
            edges_[j].push_back(i);
            E_++;
        }

        bool empty() const {
            return N_ == 0;
        }

        size_t nodes() const {
            return N_;
        }

        size_t edges() const {
            return E_;
        }

        const N & node(size_t i) const {
            return nodes_[i];
        }

        N & node(size_t i) {
            return nodes_[i];
        }

        const std::vector<uint32_t> & nedges(size_t i) const {
            return edges_[i];
        }

        std::vector<uint32_t> & nedges(size_t i) {
            return edges_[i];
        }

        N & add_node(){
            N_++;
            if(nodes_.size() < N_) {
                nodes_.push_back(N());
                edges_.push_back(std::vector<uint32_t>());
            }
            return nodes_[N_ - 1];
        }

        N & back() {
            return nodes_[N_ - 1];
        }

        const N & back() const {
            return nodes_[N_ - 1];
        }

    private:
        std::vector<std::vector<uint32_t>>         edges_;
        std::vector<N>                             nodes_;
        unsigned int                               N_;
        unsigned int                               E_;
};

template <typename N, typename A>
class Dedup {
    using umi_map = spp::sparse_hash_map<uint32_t, uint32_t>; 
    using group_hash = spp::sparse_hash_map<uint32_t, std::vector<uint32_t>>;
    Dedup(const Dedup & d) = delete;
    Dedup & operator=(const Dedup & d) = delete;

    public:

        enum CountType{
            MOLECULES = 0,
            GENES = 1,
            PCR_DUPS = 2,
            READS = 3,
            INTRONIC = 4,
            IPCR_DUPS = 5,
            IREADS = 6,
            DREADS = 7,
            ELEM_COUNT = 8
        };

        static std::string ctype2str(CountType t){
            switch(t){
                case MOLECULES: return "molecules";
                case GENES:     return "genes";
                case PCR_DUPS:  return "pcr_dups";
                case READS:     return "reads";
                case INTRONIC:  return "intron_molecules";
                case IPCR_DUPS: return "intron_pcr_dups";
                case IREADS:    return "intron_reads";
                case DREADS:    return "discarded_reads";
                default:        return "";
            }
        }

        using aiterator = typename std::vector<A>::iterator;
        Dedup(unsigned int umi_len, uint32_t total_barcodes, uint32_t num_groups, const group_hash & gh) 
            : gh_(gh), umi_len_(umi_len), num_groups_(num_groups)
        {
            molecule_counts.resize(total_barcodes * ELEM_COUNT);
            group_counts.resize(total_barcodes * num_groups_);
        }

        void process(aiterator start, aiterator end);

        std::vector<uint32_t>        molecule_counts;
        std::vector<uint32_t>        group_counts;
        std::vector<GeneCount>       gene_counts;
        std::vector<UMIMap>          umi_correct;
        std::vector<UMIBad>          umi_bad;

        unsigned int                 total_reads = 0;
        unsigned int                 total_reads2 = 0;
        unsigned int                 total_reads3 = 0;


    private:
        void process_gene_(aiterator start, aiterator end);
        void apply_corrections_(AlignSummary::bint barcode, uint32_t gene_id, aiterator start, aiterator end);
        void count_gene_(AlignSummary::bint barcode, uint32_t gene_id, aiterator gstart, aiterator gend, size_t bstart);

        void make_edges_(size_t nid, const N & node);
        void find_subgraphs_();

        DupGraph<N>                                   G_;
        umi_map                                       H_;
        std::vector<bool>                             V_;
        std::vector<uint32_t>                         SG_;
        std::vector<uint32_t>                         SC_;
        std::stack<uint32_t, std::vector<uint32_t>>   stack_;
        std::vector<uint64_t>                         positions_;
        TagFixer<A>                                   tfix_;
        size_t                                        cstart_ = 0;

        const group_hash                            & gh_;

        unsigned int                                  umi_len_;
        unsigned int                                  dcount_ = 0;
        unsigned int                                  curr_ = 0;
        uint32_t                                      num_groups_;
    
};

template <typename N, typename A>
void Dedup<N, A>::process(aiterator start, aiterator end){
    if(start == end) return;
    //std::cout << "Processing " << curr_ << " with " << (end - start) << " tags\n";
    curr_++;
    dcount_ = 0;
    auto ostart = start;
    cstart_ = umi_correct.size();
    for(auto it = std::next(start); it != end; it++){
        if(it->gene_id != start->gene_id){
            //std::cout << "  gene: " << start->gene_id << " N = " << (it - start) << "\n";
            process_gene_(start, it);
            total_reads += (it - start);
            start = it;
        }
    }
    if(start != end) {
        process_gene_(start, end);
        total_reads += (end - start);
    }

    size_t bstart = umi_bad.size();
    tfix_.process_barcode(ostart, end, umi_bad);
    std::sort(ostart, end);
    std::sort(umi_bad.begin() + bstart, umi_bad.end());

    start = ostart;
    size_t bg = bstart;
    for(auto it = std::next(start); it != end; it++){
        if(it->gene_id != start->gene_id){
            //std::cout << "  gene: " << start->gene_id << " N = " << (it - start) << " bg = " << bg << "\n";
            while(bg < umi_bad.size() && umi_bad[bg].gene_id < start->gene_id) bg++;
            count_gene_(start->barcode, start->gene_id, start, it, bg);
            start = it;
        }
    }
    if(start != end) {
        while(bg < umi_bad.size() && umi_bad[bg].gene_id < start->gene_id) bg++;
        //std::cout << "  gene: " << start->gene_id << " N = " << (end - start) << " bg = " << bg << "\n";
        count_gene_(start->barcode, start->gene_id, start, end, bg);
    }
    //std::cout << "  " << (umi_bad.size() - bstart) << " bad UMI/Gene combos to remove\n";
    if(tfix_.discarded != dcount_){
        std::cout << "  Discarded problem " << tfix_.discarded << " vs " << dcount_ << "\n";
    }

    //count_gene_(ostart->barcode, ostart->gene_id, ostart, end, bstart);
}

template <typename N, typename A>
void Dedup<N, A>::process_gene_(aiterator start, aiterator end){
    if(start == end) return;
    G_.reset();
    H_.clear();

    uint32_t gene_id = start->gene_id;
    auto barcode = start->barcode;

    for(auto it = start; it != end; it++){
        assert(start->barcode == it->barcode);
        assert(start->gene_id == it->gene_id);
        if(G_.empty() || !G_.back().merge(*it)){
            G_.add_node().init(*it);
            H_[G_.back().hash()] = G_.nodes() - 1;
        }
    }

    for(size_t i = 0; i < G_.nodes(); i++) {
        //G_.node(i).prepare();
        total_reads2 += G_.node(i).count;
    }

    for(size_t i = 0; i < G_.nodes(); i++){
        make_edges_(i, G_.node(i));
    }
    find_subgraphs_();
    apply_corrections_(barcode, gene_id, start, end);
}

template <typename N, typename A>
void Dedup<N, A>::apply_corrections_(AlignSummary::bint barcode, uint32_t gene_id, aiterator gstart, aiterator gend){

    size_t ustart = umi_correct.size();
    for(size_t i = 0; i < SC_.size() - 1; i++){
        size_t max_count = 0, total_umis = 0;
        uint32_t max_umi = 0;
        auto start = SG_.begin() + SC_[i];
        auto end = SG_.begin() + SC_[i + 1];
        size_t r = 0;
        positions_.clear();
        while(start != end){
            auto const & n = G_.node(*start);
            total_umis++;
            if(n.count > max_count) {
                max_count = n.count;
                max_umi = n.umi;
            }
            r += n.count;
            start++;
        }

        if(total_umis > 1){
            for(auto it = SG_.begin() + SC_[i]; it != end; it++){
                auto const & n = G_.node(*it);
                if(n.umi != max_umi) umi_correct.push_back(UMIMap(barcode, gene_id, n.umi, max_umi));
            }
        }
    }

    if(ustart == umi_correct.size()){
        //std::cout << "  " << barcode << " corrected 0 out of " << (gend - gstart) << "\n";
        return;
    }
    size_t corrected = 0;
    std::sort(umi_correct.begin() + ustart, umi_correct.end());
    auto rit = umi_correct.begin() + ustart;
    for(auto it = gstart; it != gend; it++){
        while(rit != umi_correct.end() && std::tie(rit->gene_id, rit->umi_from) < std::tie(it->gene_id, it->umi)){
            rit++;
        }
        if(rit != umi_correct.end() && rit->gene_id == it->gene_id && rit->umi_from == it->umi){
            it->umi = rit->umi_to;
            corrected++;
        }
    }
    //std::cout << "  " << barcode << " corrected " << corrected << " out of " << (gend - gstart) << "\n";
}

template <typename N, typename A>
void Dedup<N, A>::count_gene_(AlignSummary::bint barcode, uint32_t gene_id, aiterator start, aiterator end, size_t bstart){
    uint32_t molecules = 0, dups = 0, intronic_dups = 0, 
             intronic = 0, reads = 0, intronic_reads = 0,
             dreads = 0;

    auto lumi = start->umi;
    bool intron = false;

    auto rit = umi_bad.begin() + bstart;
    auto rend = umi_bad.end();
    //std::cout << "Bad tags " << (rend - rit) << "\n";
    //for(auto it = rit; it != rend; it++){
    //    if(it->gene_id == start->gene_id)
    //        std::cout << "    Bad " << it->barcode << " " << it->gene_id << " " << it->umi << "\n";
    //}
    positions_.clear();
    for(auto it = start; it != end; it++){
        //std::cout << "    tag = " << it->barcode << " " << it->gene_id << " " << it->umi << "\n";
        while(rit != rend && (*rit) < (*it)) {
            //std::cout << "        skip " << rit->barcode << " " << rit->gene_id << " " << rit->umi << "\n";
            rit++;
        }
        //if(rit == rend) std::cout << "        no bad\n";
        //else            std::cout << "        bad = " << rit->gene_id << " " << rit->umi << "\n";
        if(rit != rend && (*rit) == (*it)){
            dreads++;
            continue;
        }
        if(it->umi != lumi){
            if(!positions_.empty()){
                std::sort(positions_.begin(), positions_.end());
                size_t r = positions_.size();
                positions_.erase(std::unique(positions_.begin(), positions_.end()), positions_.end());
                size_t d = r - positions_.size();
                //std::cout << "    " << start->gene_id << " Molecule reads = " << r << " dups = " << d << " intronic = " << intron << "\n";
                if(!intron){
                    reads += r;
                    dups += d;
                    molecules++;
                }else{
                    intronic_reads += r;
                    intronic_dups += d;
                    intronic++;
                }
            }
            intron = false;
            positions_.clear();
            lumi = it->umi;
        }
        intron |= it->intronic;
        positions_.push_back(it->pos);
    }
    if(!positions_.empty()){
        std::sort(positions_.begin(), positions_.end());
        size_t r = positions_.size();
        positions_.erase(std::unique(positions_.begin(), positions_.end()), positions_.end());
        size_t d = r - positions_.size();
        //std::cout << "    " << start->gene_id << " Molecule reads = " << r << " dups = " << d << " intronic = " << intron << "\n";
        if(!intron){
            reads += r;
            dups += d;
            molecules++;
        }else{
            intronic_reads += r;
            intronic_dups += d;
            intronic++;
        }
    }
    //std::cout << "  Discarded = " << dreads << " vs " << tfix_.discarded << "\n";
    dcount_ += dreads;
    auto idx = barcode * CountType::ELEM_COUNT;
    molecule_counts[idx + CountType::DREADS] += dreads;
    if(intronic == 0 && molecules == 0) return;
    gene_counts.push_back(GeneCount(barcode, gene_id, molecules, intronic));

    molecule_counts[idx + CountType::MOLECULES] += molecules;
    molecule_counts[idx + CountType::PCR_DUPS] += dups;
    molecule_counts[idx + CountType::READS] += reads;
    molecule_counts[idx + CountType::INTRONIC] += intronic;
    molecule_counts[idx + CountType::IREADS] += intronic_reads;
    molecule_counts[idx + CountType::IPCR_DUPS] += intronic_dups;
    if(molecules > 0){
        molecule_counts[idx + CountType::GENES]++; // Total genes for this barcode
        auto gidx = barcode * num_groups_;
        auto git = gh_.find(gene_id);
        if(git != gh_.end()) {
            //std::cout << " barcode = " << barcode << " num_genes = " << num_groups_ << " gidx = " << gidx << " gid = " << gene_id << " molecules = " << molecules << " groups = ";
            for(auto grp_id : git->second){
                //std::cout << grp_id << " (" << (gidx + grp_id) << ") , ";
                group_counts[gidx + grp_id] += molecules;
            }
            //std::cout << "\n";
        }
    }
}


template <typename N, typename A>
void Dedup<N, A>::make_edges_(size_t nid, const N & node){
    for(size_t i = 0; i < umi_len_; i++){
        uint32_t mask = node.umi & ~getmask<uint32_t>(i, 2);
        for(uint32_t b = 0; b < 4; b++){
            uint32_t m = mask | (b << (i * 2));
            if(m != node.umi){
                auto it = H_.find(m);
                if(it == H_.end()) continue;
                auto c = G_.node(it->second).count;
                if((node.count == 1 && c == 1 && nid < it->second) || ((c > 1 || node.count > 1) && node.count >= (c * 2 - 1))){
                    //std::cout <<  "  Add Edge " << nid << " --> " << it-> second << ", "
                    //    << int2seq<ADNA4>(node.umi, umi_len_) << " --> " << int2seq<ADNA4>(G_.node(it->second).umi, umi_len_) << "\n";
                    G_.add_edge(nid, it->second);
                }
            }
        }
    }
}

template <typename N, typename A>
void Dedup<N, A>::find_subgraphs_(){
    // Each connected component represents a umi group
    //SG_ should be naturally sorted by their gene ID because the nodes in G are sorted
    V_.clear();
    V_.resize(G_.nodes(), false);
    SG_.clear();
    SC_.clear();
    SC_.push_back(SG_.size());
    for(size_t i = 0; i < V_.size(); i++){
        if(V_[i]) continue;
        stack_.push(i);
        //std::cout << "  Subgraph\n";
        while(!stack_.empty()){
            uint32_t t = stack_.top();
            stack_.pop();
            if(V_[t]) continue;
            SG_.push_back(t);
            total_reads3 += G_.node(t).count;
            //std::cout << "    " << t << " umi = " << int2seq<ADNA4>(G_.node(t).umi, umi_len_) << " visited = " << V_[t] << " edges = [";
            V_[t] = true;
            for(auto e : G_.nedges(t)){
                //std::cout << e << ", ";
                if(!V_[e]) stack_.push(e);
            }
            //std::cout << "]\n";
        }
        SC_.push_back(SG_.size());
    }
}

}
