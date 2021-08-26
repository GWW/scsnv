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
#include "sequence.hpp"
#include "parallel-hashmap/parallel_hashmap/phmap.h"
#include <cassert>
#include <iomanip>
#include <vector>

namespace gwsc {

template <typename T>
struct TagFixer{
    public:
        using tag_list = std::vector<T>;
        using bad_list = std::vector<UMIBad>;

        void process_barcode(typename tag_list::iterator start, typename tag_list::iterator end, bad_list & bad){
            if(start == end) return;

            discarded = 0;
            dmolecules = 0;
            treads = 0;
            tmols = 0;
            std::sort(start, end, [](const T & a1, const T & a2){
                return std::tie(a1.umi, a1.gene_id) < std::tie(a2.umi, a2.gene_id);
            });
            //greads_.clear();
            uint32_t lgid = start->gene_id;
            //std::cout << "  Tag " << int2seq<ADNA4, uint32_t>(start->umi, 10)
            //    << " gid = " << start->gene_id << " Group = " << utags.size() << " gcounts = " << gcounts.back() << "\n";
            utags_.assign(1, start);
            gcounts_.assign(1, 1);
            tot = 1;
            for(auto it = std::next(start); it != end; it++){
                tot++;
                if(it->umi != utags_.back()->umi){
                    utags_.push_back(it);
                    gcounts_.push_back(1);
                    lgid = it->gene_id;
                    //std::cout << "\n";
                }else if(it->gene_id != lgid){
                    gcounts_.back()++;
                    lgid = it->gene_id;
                }
                //std::cout << "  Tag " << int2seq<ADNA4, uint32_t>(it->umi, 10)
                //    << " gid = " << it->gene_id << " Group = " << utags.size() << " gcounts = " << gcounts.back() << "\n";
            }
            utags_.push_back(end);

            unsigned int multig = 0;
            for(size_t i = 0; i < (utags_.size() - 1); i++){
                if(gcounts_[i] > 1){
                     multig++;
                }
            }
            if(multig == 0) return;
            for(size_t i = 0; i < (utags_.size() - 1); i++){
                if(gcounts_[i]  == 1){
                    tmols++;
                    treads += (utags_[i + 1] - utags_[i]);
                }
            }
            correct_groups_(bad);
        }

        size_t                                         treads;
        size_t                                         tmols;
        size_t                                         tot;
        uint32_t                                       discarded = 0;
        uint32_t                                       dmolecules = 0;
    private:

        void correct_groups_(bad_list & bad){
            discarded = 0;
            dmolecules = 0;
            //std::cout << "Barcode = " << utags_.front()->barcode << "\n";
            for(size_t i = 0; i < (utags_.size() - 1); i++){
                if(gcounts_[i] <= 1) continue;
                //auto & t = (*utags_[i]);
                //std::cout << "  Tag " << int2seq<ADNA4, uint32_t>(t.umi, 10)
                //    << " N = " << (utags_[i + 1] - utags_[i])  << " gcounts = " << gcounts_[i] << " Genes =\n";
                uint32_t lgid = utags_[i]->gene_id;
                uint32_t gc = 0;
                uint32_t max_gc = 0;
                uint32_t max_gid = lgid;
                //gids.assign(1, lgid);
                for(auto it = utags_[i]; it != utags_[i + 1]; it++){
                    if(it->gene_id != lgid) {
                        //std::cout << "    " << lgid << " N = " << gc << "\n";
                        //gids.push_back(lgid);
                        if(gc == max_gc){
                            max_gid = std::numeric_limits<uint32_t>::max();
                        }else if(gc > max_gc) {
                            max_gc = gc;
                            max_gid = lgid;
                        }
                        lgid = it->gene_id;
                        gc = 0;
                    }
                    gc++;
                }

                //std::cout << "    " << lgid << " N = " << gc << "\n";
                if(gc == max_gc){
                    max_gid = std::numeric_limits<uint32_t>::max();
                }else if(gc > max_gc) {
                    max_gc = gc;
                    max_gid = lgid;
                }


                lgid = utags_[i]->gene_id;
                gc = 0;
                uint32_t dg = 0, di = 0, kg = 0, ki = 0;
                for(auto it = utags_[i]; it != utags_[i + 1]; it++){
                    if(it->gene_id != lgid){
                        if(lgid != max_gid){
                            di += gc;
                            dg++;
                            bad.push_back(UMIBad(it->barcode, lgid, it->umi));
                            //std::cout << "  Discard " << lgid << " count = " << gc << "\n";
                        }else{
                            kg++;
                            ki += gc;
                        }
                        gc = 0;
                        lgid = it->gene_id;
                    }
                    gc++;
                }
                if(lgid != max_gid){
                    bad.push_back(UMIBad(utags_[i]->barcode, lgid, utags_[i]->umi));
                    di += gc;
                    dg++;
                }else{
                    kg++;
                    ki += gc;
                }

                //std::cout << "  Total = " << (utags_[i + 1] - utags_[i]) << " ki = " << ki << " kg = " << kg << " di = " << di << " dg = " << dg << "\n";

                //std::cout << "  Discarded " << dg << " molecules from " << di << " reads\n";
                discarded += di;
                dmolecules += dg;
                treads += ki;
                tmols += kg;
                //for(size_t i = 0; i < gids.size(); i++) {
                //    if(max_gid != gids[i]) discarded += (utags_[i + 1] - utags[i]);
                //}

                /*
                std::cout << "    " << lgid << " N = " << gc << "\n";
                if(max_gid == std::numeric_limits<uint32_t>::max()){
                    std::cout << "  Fail = " << max_gc << " bad gids = ";
                    for(size_t i = 0; i < gids.size(); i++) std::cout << gids[i] << ", ";
                }else{
                    std::cout << "  Best = " << max_gid << " " << max_gc << " bad gids = ";
                    for(size_t i = 0; i < gids.size(); i++) 
                        if(max_gid != gids[i]) std::cout << gids[i] << ", ";
                }
                std::cout << "\n";
                */
            }
            /*
            //std::cout << "\n";
            std::cout << "Barcode = " << utags_.front()->barcode
                << " Reads Kept " << treads  << " out of " << tot
                << " [" << std::fixed << std::setprecision(2) << (100.0 * treads / tot) << "%]"
                << " Molecules Kept " << tmols  << " out of " << (tmols + dmolecules)
                << " [" << std::fixed << std::setprecision(2) << (100.0 * tmols / (tmols + dmolecules)) << "%]"
                ;
            */
            if((treads + discarded) != tot){
                std::cout << " Error " << ((int)tot - (int)treads - (int)discarded);
            }
            //std::cout << "\n";
        }

        std::vector<typename tag_list::const_iterator> utags_;
        std::vector<unsigned int>                      gcounts_;
        std::vector<T>                                 bad_tags_;
        //std::vector<unsigned int>                      gids;
        //phmap::flat_hash_map<uint32_t, size_t>         greads_;
};

}
