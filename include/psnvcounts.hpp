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

#include "pbase.hpp"
#include "index.hpp"
#include "tokenizer.hpp"
#include "parallel-hashmap/parallel_hashmap/phmap.h"
#include <exception>

namespace gwsc{

struct VectorHasher {
    std::size_t operator()(const std::vector<uint32_t>& v) const {
        std::size_t seed = v.size();
        for(auto & i : v){
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct SNV{
    SNV(){

    }
    SNV(int tid, int pos, char ref, char alt, char strand, unsigned int index) :
        index(index), tid(tid), pos(pos), ref(ref), alt(alt), strand(strand){

    }

    bool operator<(const SNV & s) const {
        return std::tie(pos, alt) < std::tie(s.pos, s.alt);
    }

    unsigned int index;
    int tid;
    int pos;
    char ref;
    char alt;
    char strand;
};

struct RefAlt{
    RefAlt(){
    }

    RefAlt(uint32_t snv_id, bool alt) : snv_idx(snv_id), alt(alt){
    }

    uint32_t snv_idx;
    uint32_t alt;
};

struct EdgeOut{
    EdgeOut(uint32_t i1, uint32_t i2, int32_t tid, uint32_t pos1, uint32_t pos2, 
            char ref1, char alt1, char ref2, char alt2, char strand,
            uint32_t rr, uint32_t aa, uint32_t ra, uint32_t ar)
        : i1(i1), i2(i2), tid(tid), pos1(pos1), pos2(pos2), rr(rr), aa(aa), ra(ra), ar(ar), 
        ref1(ref1), alt1(alt1), ref2(ref2), alt2(alt2), strand(strand)
    {

    }

    EdgeOut(){

    }

    bool operator<(const EdgeOut & rhs) const {
        return std::tie(tid, pos1, pos2) < std::tie(rhs.tid, rhs.pos1, rhs.pos2);
    }

    uint32_t i1;
    uint32_t i2;
    int32_t  tid;
    uint32_t pos1;
    uint32_t pos2;

    uint32_t rr = 0;
    uint32_t aa = 0;
    uint32_t ra = 0;
    uint32_t ar = 0;

    char    ref1;
    char    alt1;
    char    ref2;
    char    alt2;
    char    strand;
};

using SNVSet = std::vector<SNV>;
using SNVKey = std::vector<uint32_t>;
using SNVMap = std::unordered_map<SNVKey, uint32_t, VectorHasher >;


class ProgSNVCounts : public ProgBase {
    public:
        argagg::parser parser() const;
        std::string usage() const {
            return "scsnv snvcount -i index_prefix -s snvs.tsv -b barcode_counts.txt.gz -o out_prefix bam_in";
        }

        int run();
        void load();

    private:
        template <typename T, typename B>
        int run_();        
        void parse_snvs_();
        void write_map_(const std::string & out);
        void read_passed_(unsigned int blength);
        unsigned int count_dups_(std::string & xr);
        phmap::flat_hash_map<std::string, unsigned int>           bchash_;
        phmap::flat_hash_map<uint64_t, std::array<uint32_t, 4>>   snvpairs_;
        std::vector<std::pair<unsigned int, unsigned int>>        positions_;
        std::string              iprefix_;
        std::string              lib_;
        std::string              bcin_;
        std::string              isnvs_;
        std::string              bamin_;
        std::string              outp_;
        std::vector<std::string> barcodes_;
        std::vector<SNVSet>      snvs_;
        SNVMap                   snvmap_;
        SNVKey                   overlaps_;
        std::vector<RefAlt>      refalt_;
        Tokenizer::tokens        ttoks_;
        Tokenizer::tokens        ptoks_;
        TXIndex                  idx_;
        uint32_t                 map_idx_ = 0;
        uint32_t                 snv_count_ = 0;
        bool                     cellranger_ = false;
        //bool                     tags_ = false;
};


}
