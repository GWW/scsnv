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

#include <vector>
#include <utility>
#include <numeric>
#include <algorithm>
#include <iostream>
#include "gzstream.hpp"
#include "parallel-hashmap/parallel_hashmap/phmap.h"

namespace gwsc {

struct ReadBaseHist{
    ReadBaseHist(unsigned int max_barcodes){
        bmaxes.resize(max_barcodes);
    }

    void add_count(unsigned int barcode, unsigned int reads, unsigned int dups, unsigned int bases){
        if(bases >= bmaxes[barcode].first){
            if(bases == bmaxes[barcode].first) {
                bmaxes[barcode].second = std::min(bmaxes[barcode].second, reads);
            }else{
                bmaxes[barcode].second = reads;
                bmaxes[barcode].first = bases;
            }
        }
        reads = std::min((unsigned int)std::numeric_limits<uint16_t>::max(), reads);
        dups = std::min((unsigned int)std::numeric_limits<uint16_t>::max(), dups);
        uint32_t key = (static_cast<uint32_t>(reads) << 16) | static_cast<uint32_t>(dups);
        bases = std::min(bases, (unsigned int)1000);
        bases--;
        auto & b = counts[key];
        if(bases >= b.size()){
            b.resize(bases + 1);
        }
        b[bases]++;
    }

    void merge(const ReadBaseHist & h){
        for(size_t i = 0; i < bmaxes.size(); i++){
            if(h.bmaxes[i].first < bmaxes[i].first) continue;
            if(h.bmaxes[i].first == bmaxes[i].first) {
                bmaxes[i].second = std::min(bmaxes[i].second, h.bmaxes[i].second);
            }else{
                bmaxes[i] = h.bmaxes[i];
            }
        }

        for(auto & k : h.counts){
            auto & b = counts[k.first];
            if(k.second.size() > b.size()){
                b.resize(k.second.size());
            }
            for(size_t i = 0; i < std::min(b.size(), k.second.size()); i++) b[i] += k.second[i];
        }
    }

    void write(const std::string & out, const std::vector<std::string> & barcodes){
        {
            gzofstream os(out + "base_counts.txt.gz");
            std::vector<uint64_t> keys;
            for(auto & k : counts) keys.push_back(k.first);
            std::sort(keys.begin(), keys.end());
            os << "reads\tpcr_dups\tbases\tcount\n";
            for(auto k : keys){
                uint16_t reads = (k >> 16);
                uint16_t dups = (uint32_t)0x0000FFFF & k;
                auto & b = counts[k];
                for(size_t i = 0; i < b.size(); i++){
                    if(b[i] > 0)
                        os << reads << "\t" << dups << "\t" << (i + 1) << "\t" << b[i] << "\n";
                }
            }
        }
        {
            gzofstream os(out + "barcode_bases.txt.gz");
            os << "barcode\tmax_bases\tmin_reads\n";
            for(size_t i = 0; i < barcodes.size(); i++){
                os << barcodes[i] << "\t" << bmaxes[i].first << "\t" << bmaxes[i].second << "\n";
            }
        }

    }

    phmap::flat_hash_map<uint64_t, std::vector<uint32_t>> counts;
    std::vector<std::pair<uint32_t, uint32_t>> bmaxes;
};

}
