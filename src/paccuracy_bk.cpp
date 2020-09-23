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
          "Output file", 1},
        { "snvs", {"-s", "--snvs"},
          "Tab separated list of strand specific SNVs data from scsnvpy merge command", 1},
        { "library", {"-l", "--library"},
          "libary type (V2)", 1},
        { "aligner", {"-a", "--aligner"},
          "Aligner used the bam file (C=Cell Ranger / STAR Solo, S=scSNV, G=Genome) for example -a G", 1},
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
    iprefix_ = args_["txidx"].as<std::string>();
    bcin_ = args_["barcodes"].as<std::string>();
    iref_ = args_["reference"].as<std::string>();
    lib_ = args_["library"].as<std::string>("V2");
    outp_ = args_["output"].as<std::string>();
    if(args_.pos.size() < 1){
        throw std::runtime_error("Missing the dsc bam files");
    }


    inames_ = args_["barcodes"].as<std::string>();
    aligners_ = args_["aligners"].as<std::string>();
    for(size_t i = 0; i < args_.count(); i++){
        bams_.push_back(args_.as<std::string>(i));
    }

    if(bams_.size() != aligners_.size()){
        throw std::runtime_error("The number of dsc bam files, aligners and names must all be the same size");
    }


    tout << "Loading the genome\n";
    {
        FastaReader fr(iref_);
        fr.read_all(genome_);
    }
}


unsigned int base2int(char base){
    switch(base){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
    };
    return 4;
}

struct GetGroups{
    void operator()(int offset, int k) {
        if (k == 0) {
            combos.push_back(combination);
            return;
        }
        for (int i = offset; i <= (int)idx.size() - k; ++i) {
            combination.push_back(idx[i]);
            (*this)(i+1, k-1);
            combination.pop_back();
        }
    }

    void reset(){
        combination.clear();
        combos.clear();
    }

    std::vector<unsigned int> combination;
    std::vector<std::vector<unsigned int>> combos;
    std::vector<unsigned int> idx;
};

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
    {
        std::string tmp = header;
        Tokenizer tk(tmp, '\t');
        tk.get_all(toks);
        //for(auto & t : toks) std::cout << "Token: " << t << "\n";
        for(uint32_t i = 5; i < toks.size(); i++){
            std::string name = toks[i];
            names_.push_back(name);
        }
    }

    std::cout << "Names:\n";
    for(auto & n : names_) std::cout << "  " << n << "\n";
    std::cout << "Bams:\n";
    for(auto & n : bams_) std::cout << "  " << n << "\n";
    if(names_.size() != bams_.size()){
        throw std::runtime_error("Number of names in snv file is different from number of bam files");
    }
    std::vector<AccSNV> snvs;

    std::string tmp, line;
    std::unordered_map<uint32_t, uint32_t> count;
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
        //std::cout << snvs.size() << " " << snv.line << " tokens = " << toks.size() << " Test: ";
        for(uint32_t i = 5, j = 0; i < toks.size(); i++, j++){
            if(toks[i][0] == '1'){
                //std::cout << " " << names_[j];
                snv.founds |= 1 << j;
            }
        }
        //std::cout << " key = " << snv.founds << "\n";
        count[snv.founds]++;
        snvs.push_back(snv);
    }

    std::cout << "Total SNVs: " << snvs.size() << "\n";
    GetGroups gg;
    gg.idx.resize(names_.size());
    std::iota (std::begin(gg.idx), std::end(gg.idx), 0);

    for(size_t i = 0; i < names_.size(); i++){
        uint32_t k = (1 << i);
        std::cout << "  {" << names_[i] << "} count = " << count[k] << "\n";
    }

    for(size_t k = 2; k < names_.size(); k++){
        gg.reset();
        gg(0, k);
        for(auto & g : gg.combos){
            uint32_t x  = 0;
            for(auto i : g){
                x |= (1 << i);
            }
            std::cout << "  {";
            bool first = true;
            for(size_t i = 0; i < names_.size(); i++){
                if((x >> i) & 1){
                    if(!first) std::cout << ", ";
                    std::cout << names_[i];
                    first = false;
                }
            }
            std::cout << "} count = " << count[x] << "\n";
        }
    }

    {
        uint32_t k = 0;
        for(size_t i = 0; i < names_.size(); i++){
            k |= (1 << i);
        }
        std::cout << "  {";
        bool first = true;
        for(size_t i = 0; i < names_.size(); i++){
            if((k >> i) & 1){
                if(!first) std::cout << ", ";
                std::cout << names_[i];
                first = false;
            }
        }
        std::cout << "} count = " << count[k] << "\n";
    }

    std::vector<IndexProcessor *> processors;

    for(size_t i = 0; i < bams_.size(); i++){
        ProcessorBase * p = nullptr;
        if(aligners_[i] == 'S'){
            p = new BamScSNVProcessor(idx_, B::LibraryStrand);
            std::cout << names_[i] << " scsnv aligner = " << aligners_[i] << "\n";
        }else if(aligners_[i] == 'C'){
            p = new BamCellRangerProcessor(idx_, B::LibraryStrand);
            std::cout << names_[i] << " cellranger aligner = " << aligners_[i] << "\n";
        }else if(aligners_[i] == 'G'){
            std::cout << names_[i] << " genome aligner = " << aligners_[i] << "\n";
            BamGenomeProcessor * x = new BamGenomeProcessor(idx_, B::LibraryStrand);
            x->set_mapq(30);
            p = x;
        }else{
            std::cout << "Invalid aligner key " << aligners_[i] << "\n";
            return EXIT_FAILURE;
        }
        IndexProcessor * pp = new IndexProcessor(names_[i], p, idx_, genome_, true, bams_[i], bchash_, aligners_[i] != 'G');
        pp->make_targets(snvs);
        pp->worker().set_params(min_alternative_, min_qual_, min_barcodes_, min_edge_, splice_win_, min_af_, &pp->targets());
        processors.push_back(pp);
    }

    gzofstream zout(outp_);

    std::vector<unsigned int> oks(processors.size());
    std::vector<unsigned int> tots(processors.size());
    std::vector<unsigned int> others(processors.size());

    zout << header;
    for(auto p : processors){
        p->header(zout);
        p->start();
    }
    zout << "\n";

    std::cout << "Processors built starting to process SNVs\n";

    size_t cnt = 0, last = 0;
    for(auto & s : snvs){
        auto refi = base2int(s.ref); auto alti = base2int(s.alt);

        /*
        // Read each one separately to avoid disk issues
        for(size_t i = 0; i < processors.size(); i++){
            auto & p = *processors[i];
            p.read_position(s.chrom, s.pos);
        }
        */

        //Let each thread run its own pileup
        for(size_t i = 0; i < processors.size(); i++){
            auto & p = *processors[i];
            p.run_position(s.chrom, s.pos);
        }


        zout << s.line;
        for(size_t i = 0; i < processors.size(); i++){
            auto & p = *processors[i];
            // Wait until each thread is done
            p.wait();
            p.write(zout, refi, alti);
            bool check = (s.founds >> i) & 1;
            oks[i] += p.valid() && check;
            tots[i] += check;
            others[i] += (p.valid() && !check);
        }
        zout << "\n";
        cnt++;
        if((cnt - last) == 200){
            last = cnt;
            tout << "Finished " << cnt << " out of " << snvs.size() << " SNV sites\n";
            for(size_t i = 0; i < processors.size(); i++){
                std::cout << "    " << processors[i]->name() << " positions found and valid: " << oks[i] << " out of " << tots[i] << " non-found positions that were valid " << others[i] << "\n";
            }
        }
    }
    
    zout.close();
    tout << "Finished " << snvs.size() << " out of " << snvs.size() << " SNV sites\n";
    tout.flush();
    for(size_t i = 0; i < processors.size(); i++){
        std::cout << "  " << processors[i]->name() << " positions found and valid: " << oks[i] << " out of " << tots[i] << " non-found positions that were valid " << others[i] << "\n";
        processors[i]->done();
        delete processors[i];
    }

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

