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

#include "base.hpp"
#include "../scmap/barcodes.hpp"
#include "../scmap/reader.hpp"
#include <exception>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

namespace gwsc{

class ProgBarcodes : public ProgBase {
    public:
        argagg::parser parser() const;
        std::string    usage() const;
        void           load();
        int            run();

    private:
        template <typename R> 
        int run_10X_();

        std::string              known_;
        std::string              out_;
        std::string              lib_type_;
        std::vector<std::string> dirs_;
        FastqPairs               fastqs_;
};

template <typename R> 
int ProgBarcodes::run_10X_() {
    typename R::LibraryBarcode bcs, bcd;
    if(known_.empty()) return EXIT_FAILURE;
    bcs.load(known_);
    std::string bc, bq;
    size_t idx = 0;
    std::vector<typename R::LibraryBarcode::CountSummary> totals;
    for(auto & d : dirs_){
        tout << "Processing read directory " << d << "\n";
        auto fastqs = find_fastq_files(d);
        totals.resize(totals.size() + fastqs.size());
        std::string bout = d + "_barcode_counts.txt.gz";
        std::string sout = d + "_barcode_totals.txt";
        struct stat sb1 = {}, sb2 = {};   
        fastqs_.insert(fastqs_.end(), fastqs.begin(), fastqs.end());
        if(stat(bout.c_str(), &sb1) == 0 && stat(sout.c_str(), &sb2) == 0){
            R::LibraryBarcode::parse_summary(sout, fastqs, totals, idx);
            bcs.merge(bout);
            idx += fastqs.size();
        }else{
            std::ofstream out(sout);
            out << "dir\tfile";
            bcs.write_header(out);
            out << "\n";
            bcd.copy(bcs);
            for(auto & f : fastqs){
                R in;
                in.open(f.first, f.second);
                auto & cc = totals[idx];
                tout << "Counting barcodes for " << f.first << "\n";
                tout.flush();
                while(in.read_barcode(bc, bq)){
                    bcd.count(bc, cc);
                    bcs.count(bc);
                }
                out << f.dir << "\t" << f.prefix;
                cc.to_tab(out);
                out << "\n";
                idx++;
            }
            bcd.write(bout);
        }
    }
    auto total = typename R::LibraryBarcode::CountSummary(); 
    for(auto const & t : totals){
        total += t;
    }
    tout << total.correct << " barcodes matched out of " << total.total
        << " [" << std::fixed << std::setprecision(2) << (100.0 * total.correct / total.total) << "%]"
        << " and " << total.ambig << " contained at least one N"
        << " [" << std::fixed << std::setprecision(2) << (100.0 * total.ambig / total.total) << "%]"
        << "\n";
    {
        std::ofstream out(out_ + "_totals.txt");
        out << "dir\tfile";
        bcs.write_header(out);
        out << "\n";
        size_t i = 0;
        for(auto & f : fastqs_){
            auto cc = totals[i++];
            out << f.dir << "\t" << f.prefix << "\t" << cc.total << "\t" << cc.correct << "\t" << cc.ambig << "\n";
        }
    }

    bcs.write(out_ + "_counts.txt.gz");
    return EXIT_SUCCESS;
}

}
