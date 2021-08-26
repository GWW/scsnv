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

#include "bam_genes.hpp"
#include "fasta.hpp"
#include <thread>
#include "collapse_aux.hpp"
#include "parallel-hashmap/parallel_hashmap/phmap.h"

namespace gwsc {

class CollapseWorker {
    CollapseWorker( const CollapseWorker& ) = delete;
    CollapseWorker& operator=(const CollapseWorker&) = delete;
    public:
        using gene_hash = phmap::flat_hash_map<uint32_t, size_t>;
        using cb_bhash = std::function<AlignSummary::bint(std::string &)>;

        using clengths = std::vector<std::pair<unsigned int, unsigned int>>;
        CollapseWorker(BamOutputBuffer & cbuffer, const Fastas & genome)
            : cbuffer_(cbuffer), genome_(genome) {
        }

        void set_callback(cb_bhash cb) {
            bhash_ = cb;
        }

        ~CollapseWorker(){
            for(auto i : islands_){
                delete i;
            }
            islands_.clear();
            for(auto i : contigs_){
                delete i;
            }
            contigs_.clear();

            for(auto b : collapsed){
                delete b;
            }
            collapsed.clear();
        }

        void set_buffer(BamBuffer * buffer){
            buffer_ = buffer;
        }

        void start(BamBuffer * buffer){
            buffer_ = buffer;
            thread_ = std::thread(std::ref(*this));
        }

        void join() {
            thread_.join();
        }

        void operator()();
        void process_range(BamBuffer::rpair, unsigned int fno = std::numeric_limits<unsigned int>::max());
        void process_range_ds(BamBuffer::rpair, double ds);

        using cout_vect = std::vector<BamDetail*>;
        cout_vect                                 collapsed;

        unsigned int                              ccount = 0;
        unsigned int                              total = 0;
        unsigned int                              ambig = 0;
        unsigned int                              rlost = 0;
        unsigned int                              rreads = 0;
        unsigned int                              rdups = 0;
        unsigned int                              rcollapsed = 0;
        unsigned int                              creads = 0;

    private:
        void process_gene_(size_t start, size_t end, BamBuffer::rit it);
        void process_barcode_();
        void collapse_umi_();
        void build_contigs_();
        bool build_islands_();
        void make_read_(BamDetail & b, uint32_t NM);
        //void infer_qpos_cigar(BamDetail & b, std::stringstream & ss);

        std::thread                               thread_;
        std::string                               tmp_;
        std::vector<size_t>                       sidx_;
        phmap::flat_hash_map<uint64_t, size_t>    uhash_;
        std::vector<BamDetail*>                   barcodes_;
        std::vector<BamDetail*>                   umis_;
        std::vector<ReadContig*>                  contigs_;
        std::vector<ReadIsland*>                  islands_;
        std::vector<CollapseSplice>               fsplices_;
        std::string                               cgaps_;
        std::string                               fbases_;
        std::string                               fquals_;
        std::string                               fqname_;
        std::vector<uint8_t>                      bcov_;
        std::vector<std::pair<int, int>>          mpos_;
        CigarString                               fcigar_;
        CigarString                               ccigar_;
        cb_bhash                                  bhash_;

        BamOutputBuffer                         & cbuffer_;
        const Fastas                            & genome_;
        BamBuffer                               * buffer_;
        double                                    fds_filter = 1.1;
        unsigned int                              ccount_ = 0;
        unsigned int                              icount_ = 0;
        unsigned int                              fno_filter_ = std::numeric_limits<unsigned int>::max();
        uint32_t                                  fdups_ = 0;
        uint32_t                                  freads_ = 0;

};


}
