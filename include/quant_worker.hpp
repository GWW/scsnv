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
#include "pmap.hpp"
#include "parallel_hashmap/parallel_hashmap/phmap.h"
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
        using group_hash = phmap::flat_hash_map<uint32_t, std::vector<uint32_t> >;

        QuantBase() {

        }

        void set_index(TXIndex & index){
            index_ = &index;
        }

        void build_gene_groups(const std::string & fname);

        void set_data(std::vector<AlignSummary> & aligns, size_t total_barcodes){
            aligns_ = & aligns;
            total_barcodes_ = total_barcodes;
            astart_ = 0;
            start_ = tout.seconds();
        }

        void run(unsigned int umi_len, unsigned int num_threads, bool bam, const std::string & cmd);
        void write_output(const std::string & prefix, std::vector<std::string> & bnames, 
                std::map<AlignSummary::bint, AlignGroup::ResultCounts> & brates);
        void write_umi_map(const std::string & out_file);

        std::vector<uint32_t>           molecule_counts;
        std::vector<uint32_t>           group_counts;
        std::vector<GeneCount>          gene_counts;
        std::vector<UMIMap>             umi_correct;
        std::vector<UMIBad>             umi_bad;

    private:
        unsigned int read_(size_t N, std::vector<AlignSummary> & aligns);

        group_hash                                  group_map_;
        std::map<uint32_t, uint32_t>                gene_id_map_;
        std::vector<std::string>                    group_names_;
        std::string                                 cmd_;
        std::vector<size_t>                         bpos_;
        std::vector<std::string>                    prefixes_;
        std::mutex                                  mtx_read_;
        TXIndex                                   * index_;
        std::vector<AlignSummary>                 * aligns_;
        size_t                                      total_ = 0;
        size_t                                      ltotal_ = 0;
        size_t                                      atotal_ = 0;
        size_t                                      total_barcodes_ = 0;
        size_t                                      astart_ = 0;
        size_t                                      start_ = 0;
};

class QuantWorker {
    QuantWorker( const QuantWorker& ) = delete;
    QuantWorker& operator=(const QuantWorker&) = delete;

    public:
        friend class QuantBase;
        using cb_read = std::function<unsigned int(size_t, std::vector<AlignSummary> &)>;
        using dedup = Dedup<DupNode, AlignSummary>;
        QuantWorker(const TXIndex & idx, const QuantBase::group_hash & gh, 
                size_t NG, cb_read cb, unsigned int umi_len, unsigned int total_barcodes) 

            : cb_(cb), dups_(umi_len, total_barcodes, NG, gh), idx_(idx) 
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
