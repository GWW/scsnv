#pragma once
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
IMPLIED, INCLUDInum_groups BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRInum_groupsEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISInum_groups FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALInum_groupsS IN THE
SOFTWARE.
*/

#include "align_aux.hpp"
#include "../util/sequence.hpp"
#include "../sparsepp/sparsepp/spp.h"
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
        intronic = 0;
        positions.clear();
    }

    void init(const AlignSummary & a){
        umi = a.umi;
        count = 1;
        intronic += a.intronic;
        positions.clear();
        positions.push_back(a.pos);
    }

    bool merge(const AlignSummary & a){
        if(a.umi != umi) return false;
        count++;
        intronic += a.intronic;
        positions.push_back(a.pos);
        return true;
    }

    void prepare(){
        if(!positions.empty()){
            std::sort(positions.begin(), positions.end());
            auto last = std::unique(positions.begin(), positions.end());
            positions.erase(last, positions.end());
        }
    }

    uint32_t hash() const {
        return umi;
    }

    std::vector<uint32_t> positions;
    std::vector<uint32_t> read_ids;

    uint32_t umi;
    uint32_t count;
    uint32_t intronic;
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
    using group_hash = spp::sparse_hash_map<uint32_t, uint32_t>;
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
            ELEM_COUNT = 7
        };

        static std::string ctype2str(CountType t){
            switch(t){
                case MOLECULES: return "molecules";
                case GENES:     return "genes";
                case PCR_DUPS:  return "unique_positions";
                case READS:     return "reads";
                case INTRONIC:  return "intron_molecules";
                case IPCR_DUPS: return "intron_unique_positions";
                case IREADS:    return "intron_reads";
                default:        return "";
            }
        }

        using aiterator = typename std::vector<A>::const_iterator;
        Dedup(unsigned int umi_len, uint32_t total_barcodes, uint32_t num_groups, const group_hash & gh,bool bam) 
            : gh_(gh), umi_len_(umi_len), num_groups_(num_groups), bam_(bam)
        {
            molecule_counts.resize(total_barcodes * ELEM_COUNT);
            group_counts.resize(total_barcodes * num_groups_);
        }

        void process(aiterator start, aiterator end);

        std::vector<uint32_t>        molecule_counts;
        std::vector<uint32_t>        group_counts;
        std::vector<GeneCount>       gene_counts;
        std::vector<UMIMap>          umi_correct;

        unsigned int                 total_reads = 0;
        unsigned int                 total_reads2 = 0;
        unsigned int                 total_reads3 = 0;


    private:
        void process_gene_(aiterator start, aiterator end);
        void count_gene_(AlignSummary::bint barcode, uint32_t gene_id);

        void make_edges_(size_t nid, const N & node);
        void find_subgraphs_();

        DupGraph<N>                                   G_;
        umi_map                                       H_;
        std::vector<bool>                             V_;
        std::vector<uint32_t>                         SG_;
        std::vector<uint32_t>                         SC_;
        std::stack<uint32_t, std::vector<uint32_t>>   stack_;
        std::vector<uint32_t>                         positions_;

        const group_hash                            & gh_;

        unsigned int                                  umi_len_;
        uint32_t                                      num_groups_;
        bool                                          bam_;
    
};

template <typename N, typename A>
void Dedup<N, A>::process(aiterator start, aiterator end){
    if(start == end) return;
    for(auto it = std::next(start); it != end; it++){
        if(it->gene_id != start->gene_id){
            process_gene_(start, it);
            total_reads += (it - start);
            start = it;
        }
    }
    if(start != end) {
        process_gene_(start, end);
        total_reads += (end - start);
    }
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
        G_.node(i).prepare();
        total_reads2 += G_.node(i).count;
    }

    for(size_t i = 0; i < G_.nodes(); i++){
        make_edges_(i, G_.node(i));
    }
    find_subgraphs_();
    count_gene_(barcode, gene_id);
}

template <typename N, typename A>
void Dedup<N, A>::count_gene_(AlignSummary::bint barcode, uint32_t gene_id){
    uint32_t molecules = 0, dups = 0, intronic_dups = 0, 
             intronic = 0, reads = 0, intronic_reads = 0;

    for(size_t i = 0; i < SC_.size() - 1; i++){
        size_t max_count = 0, total_umis = 0;
        uint32_t max_umi = 0;
        auto start = SG_.begin() + SC_[i];
        auto end = SG_.begin() + SC_[i + 1];
        bool intron = false;
        size_t r = 0;
        positions_.clear();
        while(start != end){
            auto const & n = G_.node(*start);
            total_umis++;
            if(n.count > max_count) {
                max_count = n.count;
                max_umi = n.umi;
            }
            //Debug stuff
            //std::cout << "    gene_id = " << gene_id << " intronic = " << n.intronic  
            //    << " count = " << n.count << " positions = " << n.positions.front();
            //for(size_t j = 1; j < n.positions.size(); j++) std::cout << ", " << n.positions[j]; 
            //std::cout << "\n";
            if(n.intronic){
                intron = n.intronic;
            }
            r += n.count;
            positions_.insert(positions_.end(), n.positions.begin(), n.positions.end());
            if(n.count < n.positions.size()) {
                std::cout << "Position error\n";
                std::cout << "    gene_id = " << gene_id << " intronic = " << n.intronic  
                    << " count = " << n.count << " positions = " << n.positions.front();
                for(size_t j = 1; j < n.positions.size(); j++) std::cout << ", " << n.positions[j]; 
                std::cout << "\n";
            }
            start++;
        }
        //std::cout << "Total umis = " << total_umis << " total reads = " << r << " intronic = " << intron 
        //    << " dups = " << std::fixed << std::setprecision(2) << (100.0 * (r - d) / r) << "\n";

        std::sort(positions_.begin(), positions_.end());
        unsigned int d = std::unique(positions_.begin(), positions_.end()) - positions_.begin();

        if(bam_ && total_umis > 1){
            for(auto it = SG_.begin() + SC_[i]; it != end; it++){
                auto const & n = G_.node(*it);
                if(n.umi != max_umi) umi_correct.push_back(UMIMap(barcode, gene_id, n.umi, max_umi));
            }
        }

        if(intron){
            intronic_reads += r;
            intronic_dups += (r - d);
            intronic++;
        }else{
            dups += (r - d);
            reads += r;
            molecules++;
        }
    }
   
    if(intronic == 0 && molecules == 0) return;
    gene_counts.push_back(GeneCount(barcode, gene_id, molecules, intronic));
    auto idx = barcode * CountType::ELEM_COUNT;
    auto gidx = barcode * num_groups_;

    molecule_counts[idx + CountType::MOLECULES] += molecules;
    molecule_counts[idx + CountType::PCR_DUPS] += dups;
    molecule_counts[idx + CountType::READS] += reads;
    molecule_counts[idx + CountType::INTRONIC] += intronic;
    molecule_counts[idx + CountType::IREADS] += intronic_reads;
    molecule_counts[idx + CountType::IPCR_DUPS] += intronic_dups;
    if(molecules > 0){
        molecule_counts[idx + CountType::GENES]++; // Total genes for this barcode
        auto git = gh_.find(gene_id);
        if(git != gh_.end()) group_counts[gidx + git->second] += molecules;
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
