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

#include "index.hpp"
#include "reader.hpp"
#include "barcodes.hpp"
#include "gzstream.hpp"
#include "pmap.hpp"
#include "sparsepp/sparsepp/spp.h"
#include "sbam_writer.hpp"
#include "dust.hpp"
#include "transcript_align.hpp"
#include "genome_align.hpp"
#include "reader.hpp"
#include <exception>
#include <fstream>
#include <list>
#include <thread>
#include <functional>
#include <map>

namespace gwsc{

class MapWorker {
    MapWorker( const MapWorker& ) = delete;
    MapWorker& operator=(const MapWorker&) = delete;

    public:
        using cb_read = std::function<unsigned int(size_t, Reads &, const AlignGroup::ResultCounts &s)>;
        using cb_correct = std::function<int(std::string &, AlignSummary::bint &)>;

        MapWorker(unsigned int reads_per_step, cb_read in, cb_correct bc, StrandMode smode, double max_dust) 
            : counts{}, in_(in), bc_(bc), max_dust_(max_dust), rps_(reads_per_step), smode_(smode) {
        }

        ~MapWorker(){
            for(auto a : buff){
                bam_destroy1(a);
            }
            buff.clear();
            rcount = 0;
        }

        void start(){
            thread_ = std::thread(std::ref(*this));
        }

        void join() {
            thread_.join();
        }

        void prepare_bam() {
            write_bam_ = true;
            rcount = 0;
            bout->prepare_thread_buffer(buff);
        }

        void operator()();

        spp::sparse_hash_map<AlignSummary::bint, AlignGroup::ResultCounts> bc_rates;
        std::vector<AlignSummary>                                aligns;
        AlignGroup::ResultCounts                                 counts;
        SortedBamWriter::read_buffer                             buff;
        AlignGroup                                               data;

        // Could probably make these private but whatever
        SortedBamWriter                                        * bout = nullptr;
        const TranscriptAlign                                  * tx_align = nullptr;
        const GenomeAlign                                      * genome_align = nullptr;
        const TXIndex                                          * tx_idx = nullptr;
        

        unsigned int                                             barcode_correct = 0;
        unsigned int                                             barcode_corrected = 0;
        unsigned int                                             rcount = 0;

    private:
        void map_(Read & read);
        std::thread                      thread_;
        Dust                             dust_;
        Reads                            reads_;
        cb_read                          in_;
        cb_correct                       bc_;
        std::string                      bprefix_;
        double                           max_dust_;
        unsigned int                     rps_;
        StrandMode                       smode_;
        bool                             write_bam_ = false;
};


}
