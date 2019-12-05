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
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "aux.hpp"
#include "align_aux.hpp"
#include "reader.hpp"
#include "../util/interval_tree.hpp"
#include "../sparsepp/sparsepp/spp.h"
#include <unordered_map>
#include <vector>
namespace gwsc {

class TXIndex {
    TXIndex( const TXIndex& ) = delete;
    TXIndex& operator=(const TXIndex&) = delete;
    public:
        struct Ref{
            Ref() {

            }

            Ref(int tid, const std::string & name, unsigned int len) : name(name), tid(tid), len(len){

            }

            bool empty() const {
                return end == start;
            }

            std::string  name;
            using itree = IntervalTree<unsigned int, unsigned int>;
            using Splice = std::pair<unsigned int, char>;
            itree                        tree;
            std::vector<Splice>          left_splices;
            std::vector<Splice>          right_splices;
            unsigned int                 tid = 0;
            unsigned int                 len = 0;
            unsigned int                 start = 0;
            unsigned int                 end = 0;
        };

        TXIndex(){

        }

        uint32_t max_gid() const {
            return genes_.back().gid + 1;
        }

        uint32_t max_txid() const {
            return transcripts_.back().txid + 1;
        }

        void load(const std::string & prefix);
        //void load_genome(const std::string & prefix);
        //void align(const std::string & seq, AlignGroup & d) const;

        std::string tname(int tid) const {
            return refs_[tid].name;
        }

        const GeneEntry & gene(size_t i) const {
            return genes_[i];
        }

        const GeneEntry & gene_from_id(const std::string & gid) const {
            auto it = gidmap_.find(gid);
            return genes_[it->second];
        }

        uint32_t gid_from_str(const std::string & gid) const {
            auto it = gidmap_.find(gid);
            if(it == gidmap_.end()) return std::numeric_limits<uint32_t>::max();
            return genes_[it->second].gid;
        }

        const TranscriptEntry & transcript(size_t i) const {
            return transcripts_[i];
        }

        const Ref & ref(size_t i) const {
            return refs_[i];
        }

        const std::vector<Ref> & refs() const {
            return refs_;
        }

    private:
        // Project the transcript alignment to genome coordinates
        void merge_(AlignGroup & d) const;

        std::vector<GeneEntry>                          genes_;
        spp::sparse_hash_map<std::string, unsigned int> gidmap_;
        std::vector<TranscriptEntry>                    transcripts_;
        std::vector<Ref>                                refs_;
        bool                                            genome_;
};

}
