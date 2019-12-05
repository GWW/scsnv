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

#include "collapse.hpp"
#include "../scsnv/reader.hpp"
#include "../scsnv/align_aux.hpp"
#include "../util/misc.hpp"
#include "../util/gzstream.hpp"
#include "../sparsepp/sparsepp/spp.h"
#include "../muts/collapse_worker.hpp"
#include <list>
#include <thread>
#include <algorithm>
#include <numeric>

using namespace gwsc;

void ProgCollapse::load() {
    index_ = args_["txidx"].as<std::string>();
    out_ = args_["out"].as<std::string>();
    ref_ = args_["reference"].as<std::string>();
    umi_map_ = args_["umi_map"].as<std::string>("");
    lib_type_ = args_["library"].as<std::string>("V2");
    cellranger_ = args_["cellranger"];
    bam_write_threads_ = args_["bam_write"].as<unsigned int>(1);
    threads_ = args_["threads"].as<unsigned int>(1);
    bc_counts_ = args_["barcodes"].as<std::string>();
    if(args_.pos.size() != 1){
        throw std::runtime_error("Missing the prefix option");
    }

    tout << "Loading the genome\n";
    {
        FastaReader fr(ref_);
        fr.read_all(genome_);
    }

    if(!cellranger_){
        tout << "Finding bam files\n";
        for(auto & f : args_.pos){
            std::string d = f;
            if(d.back() != '/')
                d += '/';
            d += "scsnv_tmp_*.bam";
            auto bams = glob(d);
            std::cout << "  " << f << ": found " << bams.size() << " bam files\n";
            bam_files_.insert(bam_files_.end(), bams.begin(), bams.end());
        }
    }else{
        bam_files_.push_back(args_.pos[0]);
    }
    tout << "Total bam files found = " << bam_files_.size() << "\n";
}


template <typename T>
int ProgCollapse::run_wrap_(){
    std::vector<UMIMap> umi_map;
    tout << "Loading the transcriptome index\n";

    BamGeneReader<T> br(cellranger_);
    typename T::LibraryBarcode bc;
    bc.load(bc_counts_);
    br.add_bams(bam_files_.begin(), bam_files_.end());
    br.index.load(index_);
    br.prepare();


    spp::sparse_hash_map<uint32_t, size_t> ghash;
    if(!umi_map_.empty()){
        FileWrapper fin(umi_map_);
        ParserTokens toks;
        uint32_t lgid = std::numeric_limits<uint32_t>::max();
        while(fin.tokenize_line(toks) >= 0){
            AlignSummary::bint bid = bc.bid(toks[0]);
            uint32_t gid = std::stoul(toks[1]);
            uint32_t umif = std::stoul(toks[2]);
            uint32_t umit = std::stoul(toks[3]);

            if(gid != lgid){
                ghash[gid] = umi_map.size();
                lgid = gid;
            }
            umi_map.push_back(UMIMap(bid, gid, umif, umit));
        }
        tout << "Read " << umi_map.size() << " umi correct mappings\n";
    }

    CollapsedBamWriter bout(out_, bam_write_threads_, br.header());
    BamOutputBuffer cbuffer;

    BamBuffer * rbuffer = new BamBuffer();
    BamBuffer * wbuffer = new BamBuffer();
    unsigned int t = 0;
    CollapseWorker::cb_bhash cbhash = std::bind(&T::LibraryBarcode::bid, &bc, std::placeholders::_1);

    std::list<CollapseWorker> threads;


    //CollapseWorker cw(umi_map, ghash, cbhash, T::UMI_LEN, cbuffer);
    //cw.set_buffer(rbuffer);
    auto umi_len = T::UMI_LEN;

    for(size_t i = 0; i < threads_; i++){
        threads.emplace_back(umi_map, ghash, cbhash, umi_len, cbuffer, genome_, cellranger_);
    }

    bool started = false;

    unsigned int total = 0, ambig = 0, rlost = 0, rreads = 0, rdups = 0, rcollapsed = 0, creads = 0, corrected = 0;
    unsigned int lreads = 0;
    //unsigned int t2 = 0, t3 = 0, t4 = 0;
    //int lnum = -1;
    while((t = br.read_genes(*rbuffer, 100)) > 0){
        for(auto & t : threads){
            t.start(rbuffer);
        }
        total = 0; ambig = 0; rlost = 0; rreads = 0; rdups = 0; rcollapsed = 0; creads = 0; corrected = 0;
        for(auto & t : threads){
            t.join();
            total += t.total;
            ambig += t.ambig;
            rlost += t.rlost;
            rreads += t.rreads;
            rdups += t.rdups;
            rcollapsed += t.rcollapsed;
            creads += t.creads;
            corrected += t.corrected;
        }

        if((rreads - lreads) > 500000){
            tout << "Processed " << rreads << " reads from " << total << " groups,"
                << " Dups = " << std::fixed << std::setprecision(2) << (100.0 * rdups / rreads) << "%,"
                << " Collapsed = " << std::fixed << std::setprecision(2) << (100.0 * rcollapsed / rreads) << "%,"
                << " Lost = " << std::fixed << std::setprecision(6) << (100.0 * rlost / rreads) << "%"
                << " due to " << ambig << " ambiguous groups,"
                << " Collapsed Reads = " << creads << " Corrected UMIs = " << corrected << "\n"
                ;

            lreads = rreads;
        }



        if(started) bout.join();
        std::swap(rbuffer, wbuffer);
        cbuffer.ret(bout.collapsed);
        for(auto & t : threads){
            bout.collapsed.insert(bout.collapsed.end(), t.collapsed.begin(), t.collapsed.begin() + t.ccount);
            t.collapsed.erase(t.collapsed.begin(), t.collapsed.begin() + t.ccount);
            t.ccount = 0;
        }
        bout.start(wbuffer);
        started = true;
    }

    //std::cout << "t = " << t2 << " rb count = " << t3 << " buff total " << t4
    //    << " br total = " << br.total() << " bm total = " << br.mtotal() << "\n";

    CollapseWorker::clengths cl_out;
    for(auto & t : threads) cl_out.insert(cl_out.end(), t.ldata.begin(), t.ldata.end());

    std::sort(cl_out.begin(), cl_out.end());
    std::string flengths = out_ + "_lengths.txt.gz";
    gzofstream zout(flengths);
    zout << "reads\tbases\n";
    for(auto const & d : cl_out){
        zout << d.first << "\t" << d.second << "\n";
    }

    tout << "Total Reads Processed = " << std::fixed << rreads << " dups = " << std::fixed << rdups << " corrected = " << corrected << "\n";

    if(started) bout.join();
    delete rbuffer;

    if(!cellranger_){
        tout << "Deleting the temporary bam files\n";
        for(auto & b : bam_files_){
            unlink(b.c_str());
        }
    }



    return EXIT_SUCCESS;
}

int ProgCollapse::run() {
    if(lib_type_ == "V2"){
        return run_wrap_<Reader10X_V2>();
    }else if(lib_type_ == "V3"){
        return run_wrap_<Reader10X_V3>();
    }
    return EXIT_SUCCESS;
}
