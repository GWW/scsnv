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

#include "aux.hpp"
#include "align_aux.hpp"
#include "sparsepp/sparsepp/spp.h"
#include <vector>

namespace gwsc {


//For short barcodes ie. less than 16 nucleotides?
class CBWhiteListShort{
    public:

        static const char * header[];

        struct CountSummary{
            unsigned int total = 0;
            unsigned int correct = 0;
            unsigned int ambig = 0;

            CountSummary & operator+=(const CountSummary & rhs) {
                total += rhs.total;
                correct += rhs.correct;
                ambig += rhs.ambig;
                return *this;
            }

            void to_tab(std::ostream & out){
                out << "\t" << total << "\t" << correct << "\t" << ambig;
            }
        };

        static uint64_t find_total_reads(const std::string & prefix, FastqPairs & reads);
        static void parse_summary(const std::string & fname, FastqPairs & reads, std::vector<CountSummary> & counts, size_t start);

        void         load(const std::string & wlist);
        void         merge(const std::string & wlist);
        int          count(const std::string & bc);
        int          count(const std::string & bc, CountSummary & summary);
        int          correct(std::string & bc, AlignSummary::bint & code) const;
        void         write(const std::string & out) const;
        void         copy(const CBWhiteListShort & src);
        void         reset();

        size_t       size() const {
            return barcodes_.size();
        }

        void write_header(std::ostream & out){
            size_t i = 0;
            while(header[i] != NULL){
                out << "\t" << header[i];
                i++;
            }
        }

        const std::string & barcode(AlignSummary::bint bc) const {
            return barcodes_[bc];
        }

        AlignSummary::bint bid(const std::string & bc) const;

    private:
        bool get_index_(const std::string & bc, AlignSummary::bint & code) const;
        std::vector<std::string>                       barcodes_;
        std::vector<unsigned int>                      counts_;
        spp::sparse_hash_map<AlignSummary::bint, AlignSummary::bint>       hash_;
};

class CBDeNovo{
    public:
    private:
};

}
