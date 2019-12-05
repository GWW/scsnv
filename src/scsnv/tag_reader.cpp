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
#include "tag_reader.hpp"
#include "aux.hpp"
#include <numeric>

using namespace gwsc;


TagReader::TagReader (){
    start_ = tout.seconds();
}

void TagFile::load(const std::string & prefix, std::vector<uint32_t> & arates, std::vector<std::string> & barcodes) {
    // Merge the ECC ids from this sample into the global index and build a map
    {
        zin_ = gzopen((prefix + "_tags.gz").c_str(), "rb");
        gzFile zin = gzopen((prefix + "_tags_idx.gz").c_str(), "rb");

        AlignSummary::bint bc = 0, bcp = 0;
        size_t i = 0, ip = 0;
        size_t c = 0;

        if(gzread(zin, reinterpret_cast<char*>(&bcp), sizeof(bcp)) == 0 || 
                gzread(zin, reinterpret_cast<char*>(&ip), sizeof(ip)) == 0){
            std::cerr << "Error reading tag index from " << prefix << "\n";
            bindex_.clear();
            curr_ = 1;
        }

        while(gzread(zin, reinterpret_cast<char*>(&bc), sizeof(bc)) != 0 && 
                gzread(zin, reinterpret_cast<char*>(&i), sizeof(i)) != 0)
        {
            //std::cout << "bcp = " << bcp << " bc = " << bc << " i = " << i << " ip = " << ip << "\n";
            if(bindex_.size() <= bcp){
                bindex_.resize(bcp + 1);
            }
            bindex_[bcp] = i - ip;
            ip = i;
            bcp = bc;
            c++;
        }
        gzclose(zin);
        tout << "Loaded " << c << " barcodes from " << prefix << "\n";
        expected_ = i;
    }


    {
        ParserTokens toks;
        FileWrapper in(prefix + "_barcode_rates.txt.gz");
        in.tokenize_line(toks);
        while(in.tokenize_line(toks) >= 0){
            AlignSummary::bint bid = std::stoul(toks[0]);
            if(barcodes.size() <= bid){
                barcodes.resize(bid + 1);
            }
            barcodes[bid] = toks[1];
            size_t MR = (bid + 1) * AlignGroup::ELEM_COUNT;
            if(arates.size() <= MR) arates.resize(MR);
            MR -= AlignGroup::ELEM_COUNT;
            for(size_t i = 0; i < AlignGroup::ELEM_COUNT; i++){
                arates[MR + i] += std::stoul(toks[i + 2]);
            }
        }
    }
    {
        ParserTokens toks;
        FileWrapper in(prefix + "_alignment_summary.txt");
        in.tokenize_line(toks);
        while(in.tokenize_line(toks) >= 0){
            totals.push_back(std::stoul(toks[1]));
        }
    }
}

void TagReader::combine_totals() {
    totals.resize(AlignGroup::ELEM_COUNT * fnames.size());
    size_t idx = 0;
    for(auto & f : files_){
        for(size_t i = 0; i < AlignGroup::ELEM_COUNT; i++){
            totals[idx * AlignGroup::ELEM_COUNT + i] = f.totals[i];
        }
        idx++;
    }
}

void TagFile::read(AlignSummary::bint bi, std::vector<AlignSummary> & aligns) {
    if(curr_ >= bindex_.size()) return;
    if(bi < curr_){
        std::cerr << "Problem with tag reading out of sync!";
        exit(1);
    }
    size_t cnt = bindex_[bi];
    for(size_t i = 0; i < cnt; i++){
        AlignSummary s;
        auto ret = gzread(zin_, reinterpret_cast<char*>(&s), sizeof(s));
        if(ret == 0){
            std::cerr << "Less tags than expected\n";
            exit(1);
        }
        aligns.push_back(s);
    }
    curr_ = bi;
}

void TagFile::read_range(AlignSummary::bint bstart, AlignSummary::bint bend, std::vector<AlignSummary> & aligns) {
    //std::cout << "bstart = " << bstart << " bend = " << bend << " curr = " << curr_ << "\n";
    if(curr_ >= bindex_.size()) return;
    if(bstart < curr_){
        std::cerr << "Problem with tag reading out of sync!";
        exit(1);
    }
    size_t cnt = std::accumulate(bindex_.begin() + bstart, bindex_.begin() + bend, 0UL);
    for(size_t i = 0; i < cnt; i++){
        AlignSummary s;
        auto ret = gzread(zin_, reinterpret_cast<char*>(&s), sizeof(s));
        if(ret == 0){
            std::cerr << "Less tags than expected\n";
            exit(1);
        }
        aligns.push_back(s);
        total_++;
    }
    curr_ = bend;
}

void TagReader::add(const std::string & file) {
    files_.emplace_back(file, arates, barcodes);
    size_t m = files_.back().bsize();
    if(m > N_){
        N_ = m;
        for(auto & f : files_) 
            f.set_bsize(N_);
    }
    fnames.push_back(file);
}

unsigned int TagReader::read_N(size_t N, std::vector<AlignSummary> & aligns) {
    aligns.clear();
    size_t start = curr_;
    size_t end = std::min(curr_ + N, N_);
    //std::cout << "Read N " << N << " start = " << start << " end = " << end << " N_ = " << N_ << " count = " << (end - start) << "\n";
    for(auto & f : files_){
        f.read_range(start, end, aligns);
    }
    curr_ = end;
    unsigned int count = (end - start);
    total_ += count;
    if((ltotal_ + 20000) <= total_){
        size_t sec = tout.seconds();
        double ps = 1.0 * total_ / (sec - start_);
        size_t eta = (N_ - total_) / ps;
        int hours = eta / (60 * 60);
        eta -= (hours * 60 * 60);
        int minutes = eta / 60;
        eta -= (minutes * 60);
        tout << "Barcodes procesed: " <<  total_ << " / " << N_ << " [" << static_cast<int>(ps) << " / sec], ETA = "
              << hours << "h " << minutes << "m " << eta << "s\n";
        ltotal_ = total_;
    }
    return (end - start);
}

