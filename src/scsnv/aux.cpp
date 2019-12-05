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
#include <cstdlib>
#include <iostream>
#include "../util/read_buffer.hpp"

gwsc::GeneEntry parse_gene(gwsc::ParserTokens & toks) {
    gwsc::GeneEntry e;
    e.ref = toks[0];
    e.tid = std::stoi(toks[1]);
    e.gid = std::stoi(toks[2]);
    e.gene_id = toks[3];
    e.gene_name = toks[4];
    e.lft = std::stoi(toks[5]);
    e.rgt = std::stoi(toks[6]);
    e.strand = toks[7][0];
    unsigned int N = std::stoi(toks[8]);
    e.introns.resize(N);
    size_t i1 = 0;
    size_t i2 = 0;
    for(size_t i = 0; i < N; i++){
        size_t e1 = toks[9].find_first_of(',', i1);
        size_t e2 = toks[10].find_first_of(',', i2);
        if((i + 1) < N){
            toks[9][e1] = '\0';
            toks[10][e2] = '\0';
        }
        e.introns[i].lft = std::atoi(toks[9].c_str() + i1);
        e.introns[i].rgt = std::atoi(toks[10].c_str() + i2);
        i1 = e1 + 1;
        i2 = e2 + 1;
    }

    return e;
}
std::vector<gwsc::GeneEntry> gwsc::parse_genes(const std::string & prefix) {
    std::vector<GeneEntry> entries;
    ParserTokens toks;
    FileWrapper in(prefix + "_genes.txt.gz");
    in.tokenize_line(toks);
    while(in.tokenize_line(toks) >= 0){
        entries.push_back(parse_gene(toks));
    }

    return entries;
}

gwsc::TranscriptEntry parse_transcript(gwsc::ParserTokens & toks) {
    gwsc::TranscriptEntry e;
    e.tid = std::stoi(toks[1]);
    e.gid = std::stoi(toks[2]);
    e.txid = std::stoi(toks[3]);
    e.transcript_id = toks[4];
    e.transcript_name = toks[5];
    e.lft = std::stoi(toks[6]);
    e.rgt = std::stoi(toks[7]);
    e.coding_start = std::stoi(toks[12]);
    e.coding_end = std::stoi(toks[13]);
    e.strand = toks[8][0];
    unsigned int N = std::stoi(toks[9]);
    e.exons.resize(N);
    e.rexons.resize(N);
    size_t i1 = 0;
    size_t i2 = 0;
    for(size_t i = 0; i < N; i++){
        size_t e1 = toks[10].find_first_of(',', i1);
        size_t e2 = toks[11].find_first_of(',', i2);
        if((i + 1) < N){
            toks[10][e1] = '\0';
            toks[11][e2] = '\0';
        }
        e.exons[i].lft = std::atoi(toks[10].c_str() + i1);
        e.exons[i].rgt = std::atoi(toks[11].c_str() + i2);
        i1 = e1 + 1;
        i2 = e2 + 1;
    }

    unsigned int rlft = 0;
    if(e.strand == '-') std::reverse(e.exons.begin(), e.exons.end());
    e.isizes.resize(N - 1);
    for(size_t i = 0; i < N; i++){
        e.rexons[i].lft = rlft;
        e.rexons[i].rgt = rlft + e.exons[i].length() - 1;
        rlft = e.rexons[i].rgt + 1;
        if(i > 0){
            if(e.strand == '-'){
                e.isizes[i - 1] = e.exons[i - 1].lft - e.exons[i].rgt - 1;
            }else{
                e.isizes[i - 1] = e.exons[i].lft - e.exons[i - 1].rgt - 1;
            }
        }
    }
    return e;
}

std::vector<gwsc::TranscriptEntry> gwsc::parse_transcripts(const std::string & prefix, std::vector<GeneEntry> & genes) {
    std::vector<TranscriptEntry> entries;
    ParserTokens toks;
    FileWrapper in(prefix + "_transcripts.txt.gz");
    in.tokenize_line(toks);
    uint32_t lgid = std::numeric_limits<uint32_t>::max();
    while(in.tokenize_line(toks) >= 0){
        entries.push_back(parse_transcript(toks));
        if(entries.back().gid != lgid){
            if(lgid != std::numeric_limits<uint32_t>::max()){
                genes[lgid].tend = entries.back().txid;
            }
            lgid = entries.back().gid;
            genes[lgid].tstart = entries.back().txid;
        }
    }
    return entries;
}

gwsc::TaskLog gwsc::tout = gwsc::TaskLog();
