#include "index.hpp"
#include "../util/read_buffer.hpp"
#include <iostream>

using namespace gwsc;

void TXIndex::load(const std::string & prefix){
    genes_ = parse_genes(prefix);
    transcripts_ = parse_transcripts(prefix, genes_);
    {
        FileWrapper in(prefix + "_lenghts.txt");
        ParserTokens toks;
        while(in.tokenize_line(toks) >= 0){
            refs_.push_back(Ref(refs_.size(), toks[0], std::stoul(toks[1])));
        }
    }

    for(auto & t : transcripts_){
        for(size_t i = 1; i < t.exons.size(); i++){
            refs_[t.tid].left_splices.push_back({t.exons[i].lft, t.strand});
            refs_[t.tid].right_splices.push_back({t.exons[i - 1].rgt, t.strand});
        }
    }

    unsigned int ltid = 0;
    size_t i = 0;
    for(auto & g : genes_){
        if(g.tid != ltid){
            refs_[ltid].end = i;
            refs_[g.tid].start = i;
            ltid = g.tid;
        }
        i++;
    }
    refs_[ltid].end = i;

    for(auto & r : refs_){
        if(r.start == r.end) continue;
        
        std::sort(r.left_splices.begin(), r.left_splices.end());
        r.left_splices.erase(std::unique(r.left_splices.begin(), r.left_splices.end()), r.left_splices.end());

        std::sort(r.right_splices.begin(), r.right_splices.end());
        r.right_splices.erase(std::unique(r.right_splices.begin(), r.right_splices.end()), r.right_splices.end());


        Ref::itree::intervalVector values;
        for(size_t i = r.start; i < r.end; i++){
            values.push_back(Ref::itree::interval(genes_[i].lft, genes_[i].rgt, i));
        }
        r.tree = Ref::itree(values);
    }
    for(auto & g : genes_) gidmap_[g.gene_id] = g.gid;
}


