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

#pragma once

#include <vector>
#include <string>
#include <unordered_map>
#include <limits>
#include <iostream>
#include <map>
#include <tuple>
#include <functional>
#include "misc.hpp"

namespace gwsc {

struct Block {
    Block(unsigned int lft = 0, unsigned int rgt = 0, char strand = '?') : lft(lft), rgt(rgt), strand(strand) {

    }

    bool operator==(const Block &rhs) const {
        return lft == rhs.lft && rgt == rhs.rgt && strand == rhs.strand;
    }

    bool operator!=(const Block &rhs) const {
        return !(*this == rhs);
    }

    bool operator<(const Block &rhs) const {
        return std::tie(lft, rgt, strand) < std::tie(rhs.lft, rhs.rgt, rhs.strand);
    }

    bool operator>(const Block & rhs) const {
        return std::tie(lft, rgt, strand) > std::tie(rhs.lft, rhs.rgt, rhs.strand);
    }

    bool overlaps(const Block &rhs) const {
        return lft <= rhs.rgt && rhs.lft <= rgt;
    }

    bool overlaps(int rlft, int rrgt) const {
        return lft <= rrgt && rlft <= rgt;
    }

    bool contains(const Block &rhs) const {
        return lft <= rhs.lft && rgt >= rhs.rgt;
    }

    int size() const {
        return rgt - lft + 1;
    }

    int  lft;
    int  rgt;
    char strand;
};

struct ReadBlock{
    ReadBlock() {

    }

    ReadBlock(const Block & ref, const Block & query) : lft(ref.lft), rgt(ref.rgt), qlft(query.lft), qrgt(query.rgt){

    }

    bool operator<(const ReadBlock & rhs) const {
        return std::tie(lft, rgt) < std::tie(rhs.lft, rhs.rgt);
    }

    int qsize() const {
        return (qrgt - qlft + 1);
    }

    int rsize() const {
        return (rgt - lft + 1);
    }

    int lft;
    int rgt;
    int qlft;
    int qrgt;
};

struct Exon : public Block {
    Exon(unsigned int lft = 0, unsigned int rgt = 0, char strand = '?') : Block(lft, rgt, strand) {

    }
    int tlft;
    int trgt;
};

typedef std::vector<ReadBlock> ReadBlocks;


template <typename T>
struct ParentBlock : public Block{
    std::string    ref;
    std::vector<T> children;
};

struct Transcript : ParentBlock<Exon> {
    size_t tlen() const {
        size_t s = 0;
        for(auto & e : children){
            s += e.size();
        }
        return s;
    }

    std::string id;
    std::string biotype;
    std::string name;

    int cds_start = -1;
    int cds_end   = -1;
};

struct Gene : ParentBlock<Transcript> {
    std::string       id;
    std::string       biotype;
    std::string       name;
    std::vector<Exon> introns;
};

inline std::ostream & operator<<(std::ostream & os, const Block & b) {
    return os << b.lft << " - " << b.rgt << " l = " << b.size() << " [" << b.strand << "]";
}

inline std::ostream & operator<<(std::ostream & os, const Transcript & t) {
    return os << t.name << " transcript_id " << t.id << " biotype " << t.biotype << " "
           << t.ref << ": " << t.lft << " - " << t.rgt << " [" << t.strand << "]";
}

inline std::ostream & operator<<(std::ostream & os, const Gene & g) {
    return os << g.name << " gene_id " << g.id << " biotype " << g.biotype << " "
           << g.ref << ": " << g.lft << " - " << g.rgt << " [" << g.strand << "]";
}



template<typename T>
class GenericModel {
    public:
        typedef std::map<std::string, T*, ReadStringCmp>    gid_map;
        typedef T                                           value_type;
        typedef std::vector<T>                              genes;
        typedef std::unordered_map<std::string, genes>      chrom_map;
        chrom_map                                           chroms;
        std::vector<std::string>                            chrom_order;
        gid_map                                             gene_map;

    void make_map() {
        for(auto & c : chroms){
            for(auto & g : c.second){
                gene_map.emplace(g.id, &g);
            }
        }
    }
};

class GeneModel : public GenericModel<Gene> {
    public:

};

using TXMap = std::vector<Transcript*>;
using GeneMap = std::vector<Gene*>;
void sort_model(GeneModel & m);

}
