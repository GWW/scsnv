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

#include "../muts/bam_genes.hpp"
#include "../util/fasta.hpp"
#include "../sparsepp/sparsepp/spp.h"
#include <thread>
#include "collapse_aux.hpp"

namespace gwsc {

class CollapseWorker {
    CollapseWorker( const CollapseWorker& ) = delete;
    CollapseWorker& operator=(const CollapseWorker&) = delete;
    public:
        using gene_hash = spp::sparse_hash_map<uint32_t, size_t>;
        using cb_bhash = std::function<AlignSummary::bint(std::string &)>;

        using clengths = std::vector<std::pair<unsigned int, unsigned int>>;
        CollapseWorker(const std::vector<UMIMap> & umap, const gene_hash & gh, cb_bhash cb, 
                unsigned int umi_len, BamOutputBuffer & cbuffer, const Fastas & genome, bool cellranger)
            : bhash_(cb), umap_(umap), ghash_(gh), cbuffer_(cbuffer), genome_(genome), umi_len_(umi_len), cellranger_(cellranger) {
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
                bam_destroy1(b);
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


        std::vector<bam1_t*>                      collapsed;
        clengths                                  ldata;


        unsigned int                              ccount = 0;
        unsigned int                              total = 0;
        unsigned int                              ambig = 0;
        unsigned int                              rlost = 0;
        unsigned int                              rreads = 0;
        unsigned int                              rdups = 0;
        unsigned int                              rcollapsed = 0;
        unsigned int                              creads = 0;
        unsigned int                              corrected = 0;


    private:
        void process_gene_(size_t start, size_t end, BamBuffer::rit it);
        void process_barcode_(bool resort);
        void collapse_umi_();
        void build_contigs_();
        bool build_islands_();
        void make_read_(bam1_t * b, uint32_t NM);

        std::thread                               thread_;
        std::string                               tmp_;
        std::vector<size_t>                       sidx_;
        spp::sparse_hash_map<uint64_t, size_t>    uhash_;
        std::vector<BamDetail*>                   barcodes_;
        std::vector<BamDetail*>                   umis_;
        std::vector<ReadContig*>                  contigs_;
        std::vector<ReadIsland*>                  islands_;
        std::string                               fbases_;
        std::string                               fquals_;
        std::string                               fqname_;
        std::vector<uint8_t>                      bcov_;
        std::vector<unsigned int>                 fcoverage_;
        CigarString                               fcigar_;
        std::vector<CollapseSplice>               fsplices_;
        cb_bhash                                  bhash_;

        const std::vector<UMIMap>               & umap_;
        const gene_hash                         & ghash_;
        BamOutputBuffer                         & cbuffer_;
        const Fastas                            & genome_;
        BamBuffer                               * buffer_;
        unsigned int                              umi_len_;
        unsigned int                              ccount_ = 0;
        unsigned int                              icount_ = 0;
        uint32_t                                  fdups_ = 0;
        uint32_t                                  freads_ = 0;
        bool                                      cellranger_ = false;

};

}
