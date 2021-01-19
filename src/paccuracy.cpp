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
#include "paccuracy.hpp"
#include "iprocessor.hpp"
#include "bam_genes_aux.hpp"
#include <unordered_map>
#include <algorithm>
#include <zlib.h>
#include <fstream>
#include "tokenizer.hpp"
#include "gzstream.hpp"

using namespace gwsc;

argagg::parser ProgAccuracy::parser() const {
    argagg::parser argparser {{
        { "barcodes", {"-b", "--barcodes"},
          "Barcode filter file", 1},
        { "txidx", {"-i", "--index"},
          "scSNV Transcript Index", 1},
        { "reference", {"-r", "--ref"},
          "Reference Genome File", 1},
        { "output", {"-o", "--output"},
          "Output prefix", 1},
        { "snvs", {"-s", "--snvs"},
          "Tab separated list of strand specific SNVs data from scsnvpy merge command", 1},
        { "window", {"-w", "--window"},
          "Group SNVs that are within a X bp window (Default: 100)", 1},
        { "library", {"-l", "--library"},
          "libary type (V2)", 1},
        { "aligner", {"-a", "--aligner"},
          "Aligner used the bam file (C=Cell Ranger / STAR Solo, S=scSNV Collapsed, G=Genome D=scSNV Merged Bam File) for example -a G", 1},
        { "name", {"-n", "--name"},
          "Name of the tool should be a column in the SNV list from the merge command", 1},
        { "threads", {"-t", "--threads"},
          "Processing threads (Default 1)", 1},
        { "minqual", {"--min-qual"},
          "Minimum base quality to count position for genome bam file (Default 20)", 1},
        { "help", {"-h", "--help"},
          "shows this help message", 1},
      }};
    return argparser;
}

void ProgAccuracy::read_passed_(unsigned int blength){
    FileWrapper in(bcin_);
    std::string line;
    in.get_line(line);
    size_t index = 0;
    while(in.get_line(line) > -1){
        barcodes_.push_back(line.substr(0, blength));
        bchash_[barcodes_.back()] = index++;
    }
    tout << "Read " << barcodes_.size() << " passed barcodes\n";
}

void ProgAccuracy::load() {
    isnvs_ = args_["snvs"].as<std::string>();
    threads_ = args_["threads"].as<unsigned int>(1);
    iprefix_ = args_["txidx"].as<std::string>();
    bcin_ = args_["barcodes"].as<std::string>();
    iref_ = args_["reference"].as<std::string>();
    window_ = args_["window"].as<unsigned int>(100);
    lib_ = args_["library"].as<std::string>("V2");
    outp_ = args_["output"].as<std::string>();
    if(args_.pos.size() != 1){
        throw std::runtime_error("Must specify a bam file");
    }


    iname_ = args_["name"].as<std::string>();
    aligners_ = args_["aligner"].as<std::string>();
    bam_ = args_.as<std::string>(0);

    if(aligners_.size() != 1){
        throw std::runtime_error("Must be a single bam file and single aligner string");
    }


    tout << "Loading the genome\n";
    {
        FastaReader fr(iref_);
        fr.read_all(genome_);
    }
}

template <typename B>
inline int ProgAccuracy::run_(){
    idx_.load(iprefix_);
    tout << "Loading the index\n";
    read_passed_(B::BARCODE_LEN);

    unsigned int min_alternative_ = 10;
    unsigned int min_qual_ = args_["minqual"].as<unsigned int>(20);
    unsigned int min_barcodes_ = 0;
    double min_af_ = 0.01;
    unsigned int min_edge_ = 5;
    unsigned int splice_win_ = 10;
 

    FileWrapper in(isnvs_);
    std::string header;
    in.get_line(header);
    std::cout << "Header: " << header << "\n";
    Tokenizer::tokens toks;
    size_t found = std::numeric_limits<size_t>::max();
    {
        std::string tmp = header;
        Tokenizer tk(tmp, '\t');
        tk.get_all(toks);
        //for(auto & t : toks) std::cout << "Token: " << t << "\n";
        for(uint32_t i = 5; i < toks.size(); i++){
            std::string name = toks[i];
            if(name == iname_){
                found = i;
            }
        }
    }

    if(found == std::numeric_limits<size_t>::max()){
        throw std::runtime_error("The name argument must match a column in the SNV file");
    }

    DataManager data;

    std::string tmp, line;
    unsigned int sidx = 0;
    while(in.get_line(line) >= 0){
        AccSNV snv;
        snv.line = line;
        Tokenizer tk(line, '\t');;
        tk.get_all(toks);
        std::string chrom = toks[0];
        unsigned int pos = std::stoul(toks[1]);
        snv.chrom = toks[0];
        snv.pos = pos;
        snv.ref = toks[2][0];
        snv.alt = toks[3][0];
        snv.strand = toks[4][0];
        snv.found = toks[found][0] == '1';
        snv.index = sidx++;
        data.snvs.push_back(snv);
    }


    {
        data.groups.push_back(AccSNVGroup(data.snvs.begin(), data.snvs.begin(), data.snvs.front().pos));
        for(auto it = std::next(data.snvs.begin()); it != data.snvs.end(); it++){
            auto & s = *it;
            if(s.chrom != data.groups.back().start->chrom || (it->pos - data.groups.back().end_pos) > window_){
                data.groups.back().end = it;
                data.groups.push_back(AccSNVGroup(it, it, it->pos));
            }else{
                data.groups.back().end_pos = it->pos;
            }
        }
        data.groups.back().end = data.snvs.end();
    }

    /*
    for(size_t i = 1; i < data.groups.size(); i++){
        if(data.groups[i - 1].start->chrom == data.groups[i].start->chrom){
            int d = data.groups[i].start->pos  - data.groups[i - 1].end_pos;
            if(d < (int)window_){
                std::cout << "Group problem d = " << d << " vs " << window_ << "\n"
                    << "  " << data.groups[i - 1].start->pos << " - " << data.groups[i - 1].end_pos << " snvs: ";
                for(auto it = data.groups[i - 1].start; it != data.groups[i - 1].end; it++){
                    std::cout << it->pos << " ";
                }
                std::cout << "\n"
                    << "  " << data.groups[i].start->pos << " - " << data.groups[i].end_pos << " snvs: ";
                for(auto it = data.groups[i].start; it != data.groups[i].end; it++){
                    std::cout << it->pos << " ";
                }
                std::cout << "\n";
            }
        }
    }
    */

    std::cout << "Total SNVs: " << data.snvs.size() << " in " << data.groups.size() << " overlap groups\n";
    data.processors.resize(threads_);
    data.name = iname_;
    for(size_t i = 0; i < threads_; i++){
        ProcessorBase * p = nullptr;
        if(aligners_[0] == 'S'){
            p = new BamScSNVProcessor(idx_, B::LibraryStrand);
        }else if(aligners_[0] == 'D'){
            p = new BamScSNVProcessor(idx_, B::LibraryStrand);
        }else if(aligners_[0] == 'C'){
            p = new BamCellRangerProcessor(idx_, B::LibraryStrand);
        }else if(aligners_[0] == 'G'){
            BamGenomeProcessor * x = new BamGenomeProcessor(idx_, B::LibraryStrand);
            x->set_mapq(20);
            p = x;
        }else{
            std::cout << "Invalid aligner key " << aligners_ << "\n";
            return EXIT_FAILURE;
        }


        IndexProcessor * pp = new IndexProcessor(iname_, p, idx_, genome_, aligners_[0] != 'D', bam_, bchash_, barcodes_.size(), aligners_[0] != 'G');
        pp->make_targets(data.snvs);
        pp->worker().set_params(min_alternative_, min_qual_, min_barcodes_, min_edge_, splice_win_, min_af_, &pp->targets());
        data.processors[i] = pp;
    }

    std::cout << "Running\n";
    data.run();
    data.write(outp_, header, barcodes_);
    return 0;
}

int ProgAccuracy::run(){
        if(lib_ == "V2"){
            return run_<Reader10X_V2>();
        }else if(lib_ == "V3"){
            return run_<Reader10X_V3>();
        }else if(lib_ == "V2_5P"){
            return run_<Reader10X_V2_5P>();
        }else if(lib_ == "V3_5P"){
            return run_<Reader10X_V3_5P>();
        }
    return 0;
}

