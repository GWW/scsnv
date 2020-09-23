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
#include "barcodes.hpp"
#include "sequence.hpp"
#include "read_buffer.hpp"
#include "gzstream.hpp"
#include <fstream>
#include <unordered_map>

using namespace gwsc;

const char* CBWhiteListShort::header[] = { "total", "correct", "ambig", NULL};

void CBWhiteListShort::load(const std::string & wlist) {
    FileWrapper in(wlist);
    ParserTokens toks;
    AlignSummary::bint idx = 0;
    while(in.tokenize_line(toks) >= 0){
        uint64_t code = 0;
        if(!seq2int<gwsc::ADNA4, uint64_t>(toks[0], code)) continue;
        auto res = hash_.insert(std::make_pair(code, 0));
        if(res.second){
            res.first->second = idx++;
            barcodes_.push_back(toks[0]);
            counts_.push_back(toks.size() == 1 ? 0 : std::stoi(toks[1]));
        }
    }
    tout << "Loaded " << idx << " known barcodes from " << wlist << "\n";
}

void CBWhiteListShort::merge(const std::string & wlist) {
    FileWrapper in(wlist);
    ParserTokens toks;
    AlignSummary::bint idx = 0;
    while(in.tokenize_line(toks) >= 0){
        uint64_t code = 0;
        if(!seq2int<gwsc::ADNA4, uint64_t>(toks[0], code)) continue;
        //auto res = hash_.insert(std::make_pair(code, 0));
        auto it = hash_.find(code);
        if(it != hash_.end()){
            counts_[it->second] += std::stoi(toks[1]);
            idx++;
        }
    }
    tout << "Merged " << idx << " known barcodes from " << wlist << "\n";
}

void CBWhiteListShort::reset() {
    counts_.clear();
    counts_.resize(barcodes_.size());
}

void CBWhiteListShort::copy(const CBWhiteListShort & src) {
    hash_ = src.hash_;
    barcodes_ = src.barcodes_;
    reset();
}

int CBWhiteListShort::count(const std::string & bc) {
    uint64_t code = 0;
    bool res = seq2int<gwsc::ADNA4, uint64_t>(bc, code);
    if(!res) {
        return 1;
    }
    auto it = hash_.find(code);
    if(it != hash_.end()){
        counts_[it->second]++;
        return 0;
    }
    return 2;
}


int CBWhiteListShort::count(const std::string & bc, CountSummary & summary) {
    uint64_t code = 0;
    summary.total++;
    bool res = seq2int<gwsc::ADNA4, uint64_t>(bc, code);
    if(!res) {
        summary.ambig++;
        return 1;
    }
    auto it = hash_.find(code);
    if(it != hash_.end()){
        counts_[it->second]++;
        summary.correct++;
        return 0;
    }
    return 2;
}

AlignSummary::bint CBWhiteListShort::bid(const std::string & bc) const {
    uint64_t code = 0;
    seq2int<gwsc::ADNA4, uint64_t>(bc, code);
    auto it = hash_.find(code);
    return it->second;
}

int CBWhiteListShort::correct(std::string & bc, AlignSummary::bint & index) const {
    uint64_t code = 0;
    bool res = seq2int<gwsc::ADNA4, uint64_t>(bc, code);
    if(!res) return 1;
    auto it = hash_.find(code);
    if(it != hash_.end()) {
        index = it->second;
        return 0;
    }

    //Need to try every combination of mismatch
    unsigned int maxc = 0;
    index = 0;
    for(size_t i = 0; i < bc.size(); i++){
        uint64_t mask = code & ~getmask<uint64_t>(i, 2);
        //Zero out the two bits for base i
        for(uint64_t b = 0; b < 4; b++){
            //shift over the base and set the two bits in the mask
            uint64_t m = mask | (b << (i * 2));
            it = hash_.find(m);
            if(it != hash_.end() && counts_[it->second] > maxc){
                maxc = counts_[it->second];
                index = it->second;
            }
        }
    }
    if(maxc > 0){
        // Correct the barcode string
        bc = barcode(index);
    }

    return maxc > 0 ? 2 : 1;
}

void CBWhiteListShort::write(const std::string & out) const {
    gzofstream os(out);
    size_t cnt = 0;
    for(size_t i = 0; i < barcodes_.size(); i++)
        if(counts_[i] > 0) {
            os << barcodes_[i] << "\t" << counts_[i] << "\n";
            cnt++;
        }
    tout << "Wrote " << cnt << " non-zero barcodes\n";
}

void gwsc::CBWhiteListShort::parse_summary(const std::string & fname, FastqPairs & reads, std::vector<CountSummary> & counts, size_t start) {
    FileWrapper in(fname);
    ParserTokens toks;
    std::unordered_map<std::string, uint32_t> cmap;
    in.tokenize_line(toks);
    uint32_t idx = 0;
    for(auto & f : reads) cmap[f.prefix] = idx++;
    while(in.tokenize_line(toks) >= 0){
        idx = cmap[toks[1]];
        counts[start + idx].total = stoul(toks[2]);
        counts[start + idx].correct = stoul(toks[3]);
        counts[start + idx].ambig = stoul(toks[4]);
    }
}

uint64_t gwsc::CBWhiteListShort::find_total_reads(const std::string & prefix, FastqPairs & reads) {
    FileWrapper in(prefix + "_totals.txt");
    ParserTokens toks;
    std::unordered_map<std::string, uint32_t> cmap;
    in.tokenize_line(toks);
    while(in.tokenize_line(toks) >= 0){
        uint32_t cnt = stoul(toks[2]);
        std::string tt = toks[0];
        tt += "_";
        tt += toks[1];
        cmap[tt] = cnt;
    }
    uint32_t total = 0;
    for(auto & f : reads){
        auto const & res = cmap[f.dir + "_" + f.prefix];
        f.total = res;
        total += res;
    }
    return total;
}
