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

#include <string>
#include "../util/read_buffer.hpp"
#include "../util/annotation.hpp"
#include "../util/task_log.hpp"

namespace gwsc {

struct PosEntry {
    unsigned int length() const {
        return rgt - lft + 1;
    }
    unsigned int lft;
    unsigned int rgt;
};

struct GeneEntry {
    std::string               ref;
    std::string               gene_id;
    std::string               gene_name;
    std::vector<PosEntry>     introns;
    
    unsigned int          tid;
    unsigned int          gid;
    unsigned int          lft;
    unsigned int          rgt;
    unsigned int          tstart;
    unsigned int          tend;
    char                  strand;
};

struct TranscriptEntry {
    std::string               transcript_id;
    std::string               transcript_name;
    std::vector<PosEntry>     exons;
    std::vector<PosEntry>     rexons;
    std::vector<int>          isizes;
    
    unsigned int          tid;
    unsigned int          gid;
    unsigned int          txid;
    unsigned int          lft;
    unsigned int          rgt;
    int                   coding_start;
    int                   coding_end;
    char                  strand;
};

struct FastqPair{
    FastqPair() : total(0) {

    }

    FastqPair(const std::string & f, const std::string & s, const std::string & dir, uint32_t total = 0) : 
        first(f), second(s), dir(dir), total(total) {
        prefix =  first.substr(0, first.find('.'));
    }

    bool operator<(const FastqPair & f) const{
        ReadStringCmp cmp;
        return cmp(dir, f.dir) ? true : cmp(first, f.first);
    }

    std::string first;
    std::string second;
    std::string dir;
    std::string i1;
    std::string i2;
    std::string prefix;
    uint32_t    total;
};
using FastqPairs = std::vector<FastqPair>;

std::vector<GeneEntry> parse_genes(const std::string & prefix);
std::vector<TranscriptEntry> parse_transcripts(const std::string & prefix, std::vector<GeneEntry> & genes);


extern gwsc::TaskLog tout;

}
