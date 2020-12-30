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
#include "index.hpp"
#include "read_buffer.hpp"
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

    size_t ridx = 0;
    for(auto & r : refs_){
        ref_map_[r.name] = ridx++;
        if(r.start == r.end) continue;
        Ref::itree::intervalVector values;
        for(size_t i = r.start; i < r.end; i++){
            values.push_back(Ref::itree::interval(genes_[i].lft, genes_[i].rgt, i));
        }
        r.tree = Ref::itree(values);
    }
    for(auto & g : genes_) gidmap_[g.gene_id] = g.gid;
}

void TXIndex::build_splice_index(unsigned int overhang){
    for(auto & r : refs_){
        if(r.start == r.end) continue;

        Ref::itree::intervalVector plfts, prgts, mlfts, mrgts;
        for(size_t i = r.start; i < r.end; i++){
            for(size_t j = genes_[i].tstart; j < genes_[i].tend; j++){
                auto & t = transcripts_[j];
                if(t.strand == '+'){
                    for(size_t i = 1; i < t.exons.size(); i++){
                        plfts.push_back({t.exons[i].lft - overhang, t.exons[i].lft - 1, t.gid});
                        prgts.push_back({t.exons[i - 1].rgt + 1, t.exons[i - 1].rgt + overhang, t.gid});
                    }
                }else{
                    for(size_t i = 1; i < t.exons.size(); i++){
                        mlfts.push_back({t.exons[i - 1].lft - overhang, t.exons[i - 1].lft - 1, t.gid});
                        mrgts.push_back({t.exons[i].rgt + 1, t.exons[i].rgt + overhang, t.gid});
                    }
                }
            }
        }
        std::sort(plfts.begin(), plfts.end());
        std::sort(prgts.begin(), prgts.end());
        std::sort(mlfts.begin(), mlfts.end());
        std::sort(mrgts.begin(), mrgts.end());

        plfts.erase(std::unique(plfts.begin(), plfts.end()), plfts.end());
        prgts.erase(std::unique(prgts.begin(), prgts.end()), prgts.end());
        mlfts.erase(std::unique(mlfts.begin(), mlfts.end()), mlfts.end());
        mrgts.erase(std::unique(mrgts.begin(), mrgts.end()), mrgts.end());

        r.lm_splices = Ref::itree(mlfts);
        r.rm_splices = Ref::itree(mrgts);
        r.lp_splices = Ref::itree(plfts);
        r.rp_splices = Ref::itree(prgts);
    }
}

void TXIndex::build_splice_site_index(){
    for(auto & r : refs_){
        if(r.start == r.end) continue;

        for(size_t i = r.start; i < r.end; i++){
            for(size_t j = genes_[i].tstart; j < genes_[i].tend; j++){
                auto & t = transcripts_[j];
                if(t.strand == '+'){
                    for(size_t i = 1; i < t.exons.size(); i++){
                        r.plus_lsplices.push_back(t.exons[i].lft - 1);
                        r.plus_rsplices.push_back(t.exons[i - 1].rgt + 1);
                    }
                }else{
                    for(size_t i = 1; i < t.exons.size(); i++){
                        r.minus_lsplices.push_back(t.exons[i - 1].lft - 1);
                        r.plus_rsplices.push_back(t.exons[i].rgt + 1);
                    }
                }
            }
        }

        std::sort(r.plus_lsplices.begin(), r.plus_lsplices.end());
        r.plus_lsplices.erase(std::unique(r.plus_lsplices.begin(), r.plus_lsplices.end()), r.plus_lsplices.end());

        std::sort(r.plus_rsplices.begin(), r.plus_rsplices.end());
        r.plus_rsplices.erase(std::unique(r.plus_rsplices.begin(), r.plus_rsplices.end()), r.plus_rsplices.end());

        std::sort(r.minus_lsplices.begin(), r.minus_lsplices.end());
        r.minus_lsplices.erase(std::unique(r.minus_lsplices.begin(), r.minus_lsplices.end()), r.minus_lsplices.end());

        std::sort(r.minus_rsplices.begin(), r.minus_rsplices.end());
        r.minus_rsplices.erase(std::unique(r.minus_rsplices.begin(), r.minus_rsplices.end()), r.minus_rsplices.end());
    }
}

