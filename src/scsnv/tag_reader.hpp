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

#include <unordered_map>
#include <vector>
#include <iostream>
#include <list>
#include <zlib.h>
#include "align_aux.hpp"

namespace gwsc {
/*
TODO: 
1) Merge all the barcodes so that each thread can read X barcodes in order
2) Merge the gene classes and create maps for each file
*/


class TagFile {
    TagFile( const TagFile& ) = delete;
    TagFile& operator=(const TagFile&) = delete;
    public:
        TagFile(const std::string & prefix, std::vector<uint32_t> & arates, 
                std::vector<std::string> & barcodes) {
            load(prefix, arates, barcodes);
        }

        void load(const std::string & prefix, std::vector<uint32_t> & arates, 
                std::vector<std::string> & barcodes);

        size_t barcodes() const {
            return bindex_.size();
        }

        void set_bsize(size_t bsize){
            if(bindex_.size() < bsize)
                bindex_.resize(bsize);
        }

        size_t bsize() const {
            return bindex_.size();
        }

        size_t total() const {
            return total_;
        }

        size_t expected() const {
            return expected_;
        }

        size_t curr() const {
            return curr_;
        }

        void read(AlignSummary::bint bi, std::vector<AlignSummary> & aligns);
        void read_range(AlignSummary::bint bstart, AlignSummary::bint bend, std::vector<AlignSummary> & aligns);

        std::vector<uint32_t> totals;
    private:
        gzFile                                 zin_ = nullptr;
        std::vector<size_t>                    bindex_;
        size_t                                 curr_ = 0;
        size_t                                 total_ = 0;
        size_t                                 expected_ = 0;

};

class TagReader {
    TagReader( const TagReader& ) = delete;
    TagReader& operator=(const TagReader&) = delete;
    public:
        TagReader ();

        void add(const std::string & file);

        size_t total_barcodes() const {
            return N_;
        }

        unsigned int read_N(size_t N, std::vector<AlignSummary> & aligns);

        void combine_totals();

        void debug() const {
            size_t i = 0;
            for(auto & f : files_){
                std::cout << "Tag file " << fnames[i] << " total read = " << f.total() 
                          << " exptected = " << f.expected() << " curr = " << f.curr() 
                          << " bsize = " << f.bsize() << "\n";
                i++;
            }
        }

        std::vector<uint32_t>    arates;
        std::vector<uint32_t>    totals;
        std::vector<std::string> barcodes;
        std::vector<std::string> fnames;
    private:
        std::list<TagFile>    files_;
        size_t                N_ = 0;
        size_t                total_ = 0;
        size_t                ltotal_ = 0;
        size_t                curr_ = 0;
        size_t                start_ = 0;
};

}
