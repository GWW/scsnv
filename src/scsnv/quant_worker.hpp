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

#include "index.hpp"
#include "tag_reader.hpp"
#include "../progs/map.hpp"
#include "../sparsepp/sparsepp/spp.h"
#include "dups.hpp"
#include <exception>
#include <list>
#include <thread>
#include <functional>
#include <map>

namespace gwsc {

// In the future maybe make this a template?
class QuantBase {
    public:
        using group_hash = spp::sparse_hash_map<uint32_t, uint32_t>;

        QuantBase() {

        }

        void load_index(const std::string & txindex){
            index_.load(txindex);
        }

        void prepare_tags(const std::vector<std::string> & prefixes);
        void build_gene_groups(const std::string & fname);
        void run(unsigned int umi_len, unsigned int num_threads, unsigned int min_molecules, bool bam, const std::string & cmd);
        void write_output(const std::string & prefix, unsigned int min_molecules);

        std::vector<uint32_t>           molecule_counts;
        std::vector<uint32_t>           group_counts;
        std::vector<GeneCount>          gene_counts;

    private:
        unsigned int read_(size_t N, std::vector<AlignSummary> & aligns);
        void write_umi_map_(const std::string & out_file);

        group_hash                                  group_map_;
        std::map<uint32_t, uint32_t>                gene_id_map_;
        std::vector<std::string>                    group_names_;
        std::vector<UMIMap>                         umi_correct_;
        TagReader                                   in_;
        TXIndex                                     index_;
        std::string                                 cmd_;
        std::vector<size_t>                         bpos_;
        std::vector<std::string>                    prefixes_;
        std::mutex                                  mtx_read_;
        size_t                                      total_ = 0;
};

class QuantWorker {
    QuantWorker( const QuantWorker& ) = delete;
    QuantWorker& operator=(const QuantWorker&) = delete;

    public:
        friend class QuantBase;
        using cb_read = std::function<unsigned int(size_t, std::vector<AlignSummary> &)>;
        using dedup = Dedup<DupNode, AlignSummary>;
        QuantWorker(const TXIndex & idx, const QuantBase::group_hash & gh, 
                size_t NG, cb_read cb, unsigned int umi_len, unsigned int total_barcodes, bool bam) 

            : cb_(cb), dups_(umi_len, total_barcodes, NG, gh, bam), idx_(idx) 
        {

        }


        void start(){
            thread_ = std::thread(std::ref(*this));
        }

        void join() {
            thread_.join();
        }

        void operator()();

        size_t                    total_molecules = 0;
        size_t                    total_umis = 0;
        size_t                    total_reads = 0;

    private:
        void process_(std::vector<AlignSummary>::iterator bstart, std::vector<AlignSummary>::iterator bend);
        std::thread                                 thread_;
        cb_read                                     cb_;
        std::vector<AlignSummary>                   aligns_;
        dedup                                       dups_;
        const TXIndex                             & idx_;
        size_t                                      NG_;
};

}
