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

#include "ppileup.hpp"
#include "reader.hpp"
#include "align_aux.hpp"
#include "misc.hpp"
#include "read_buffer.hpp"
#include "gzstream.hpp"
#include "pileup_worker.hpp"
#include "h5misc.hpp"
#include <list>
#include <fstream>
#include <thread>
#include <algorithm>
#include <numeric>

using namespace gwsc;

void ProgPileup::load() {
    index_ = args_["txidx"].as<std::string>();
    out_ = args_["out"].as<std::string>();
    ref_ = args_["reference"].as<std::string>();
    lib_type_ = args_["library"].as<std::string>("V2");
    threads_ = args_["threads"].as<unsigned int>(1);
    rthreads_ = args_["rthreads"].as<unsigned int>(1);

    //min_coverage_ = args_["mincov"].as<unsigned int>(min_coverage_);
    min_alternative_ = args_["minalt"].as<unsigned int>(min_alternative_);
    min_qual_ = args_["minqual"].as<unsigned int>(min_qual_);
    min_barcodes_ = args_["minbarcodes"].as<unsigned int>(min_barcodes_);
    min_af_ = args_["minaf"].as<double>(min_af_);
    min_edge_ = args_["minedge"].as<unsigned int>(min_edge_);
    splice_win_ = args_["splicewin"].as<unsigned int>(splice_win_);
    passed_ = args_["passed"].as<std::string>();
    read_genes_ = args_["genes"].as<unsigned int>(read_genes_);
    if(args_["snvlist"])
        snvlist_ = args_["snvlist"].as<std::string>();
    cellranger_ = args_["cellranger"];
    dups_ = args_["dups"];
    if(cellranger_ && !dups_){
        tout << "WARNING: Processing Cell Ranger without the duplicate flag enabled will produce a massive false negative rate for SNV calling!!!!";
    }
    if(args_.pos.size() != 1){
        throw std::runtime_error("Missing the bam file argument");
    }

    tout << "Loading the genome\n";
    {
        FastaReader fr(ref_);
        fr.read_all(genome_);
    }
    bam_file_ = args_.pos[0];

    if(!snvlist_.empty()) parse_targets_();
}

void ProgPileup::read_passed_(){
    FileWrapper in(passed_);
    std::string line;
    in.get_line(line);
    size_t index = 0;
    while(in.get_line(line) > -1){
        barcodes_.push_back(line);
        bchash_[line] = index++;
    }
    tout << "Read " << barcodes_.size() << " passed barcodes\n";
}

template <typename T, typename P>
int ProgPileup::run_wrap_(){
    read_passed_();
    tout << "Loading the transcriptome index\n";

    typename BamGeneReaderFiltered<T, BamReader, P>::bcfilter_func fp = std::bind(&ProgPileup::filter_func, *this, std::placeholders::_1, std::placeholders::_2);
    BamGeneReaderFiltered<T, BamReader, P> br(fp); //br(cellranger_);

    br.index.load(index_);
    br.index.build_splice_site_index();
    std::cout << "Splice site index built\n";

    br.set_bam(bam_file_);
    if(rthreads_ > 1){
        br.set_threads(rthreads_);
    }
    br.prepare();

    BamBuffer * rbuffer = new BamBuffer();
    BamBuffer * pbuffer = new BamBuffer();
    std::vector<PileupWorker*> threads;
    if(threads_ < 2) threads_ = 2;
    for(size_t i = 0; i < threads_ - 1; i++){
        threads.push_back(new PileupWorker(barcodes_.size(), br.index, genome_, dups_));
        threads.back()->set_params(min_alternative_, min_qual_, min_barcodes_, min_edge_, splice_win_, min_af_, 
                (targets_.empty() ? nullptr : &targets_));
    }

    unsigned int pbases = 0, bases = 0, reads = 0;
    unsigned int plus_bases = 0, minus_bases = 0;
    unsigned int lreads = 0;
    std::vector<PositionCount*> buffer;
    //std::vector<PositionCoverage*> cbuffer;
    gzofstream os(out_ + ".txt.gz");
    std::vector<BBout> bouts;
    bouts.resize(4);

    //cos << "chrom\tpos\tplus_coverage\tminus_coverage\n";
    os << "chrom\tpos\tcoverage\tbarcodes\tfailed_count\tref\tmax_non_ref\tplus_base\tminus_base\tplus_donor_dist\tplus_acceptor_dist\tminus_donor_dist\tminus_acceptor_dist";
    for(size_t i = 0; i < 4; i++){
        char bs = "ACGT"[i];
        os << "\t" 
            << bs << "_total\t" 
            << bs << "_total_barcodes";
    }
    for(auto s : {"_plus", "_minus"}){
        for(size_t i = 0; i < 4; i++){
            std::string bs(1, "ACGT"[i]);
            bs += s;
            os << "\t" 
                << bs << "_counts\t"
                << bs << "_end_dist";
        }
    }
    os << "\n";

    std::cout << "Min af = " << min_af_ << "\n";
    unsigned int tot = br.read_genes(*rbuffer, read_genes_);
    std::vector<unsigned int> barcode_molecules(barcodes_.size());
    auto rall = rbuffer->get_all_unsafe();
    for(auto it = rall.first; it != rall.second; it++) barcode_molecules[(*it)->barcode]++;

    auto start_time = tout.seconds();

    while(tot > 0){
        std::swap(rbuffer, pbuffer);
        for(auto it = threads.begin(); it != threads.end(); it++){
            auto & t = *(*it);
            t.set_buffer(pbuffer);
            t.start();
        }

        // While the threads are processing the buffer read in the next set to speed things up
        tot = br.read_genes(*rbuffer, read_genes_);
        rall = rbuffer->get_all_unsafe();
        for(auto it = rall.first; it != rall.second; it++) barcode_molecules[(*it)->barcode]++;

        pbases = 0; bases = 0; reads = 0;
        plus_bases = 0, minus_bases = 0;
        for(auto it = threads.begin(); it != threads.end(); it++){
            (*it)->join();
        }

        buffer.clear();
        //cbuffer.clear();
        for(auto tp : threads){
            auto & t = *(tp);
            for(size_t i = 0; i < t.pcount; i++){
                buffer.push_back(&t.positions[i]);
            }

            coverage_.insert(coverage_.end(), t.coverage.begin(), t.coverage.end());
            bases += t.bases;
            plus_bases += t.plus_bases;
            minus_bases += t.minus_bases;
            reads += t.reads;
            pbases += t.pbases;
        }

        //std::cout << "Sorting buffer size = " << buffer.size() << "\n";
        std::sort(buffer.begin(), buffer.end(), [](const PositionCount * p1, PositionCount * p2) { return (*p1) < (*p2); });
        for(uint32_t i = 0; i < buffer.size(); i++){
            auto const & p = *buffer[i];
            size_t snp_id = tids_.size();
            tids_.push_back(p.tid);
            pos_.push_back(p.pos);
            refs_.push_back(p.ref);
            for(auto & b : p.bcounts){
                for(size_t k = 0; k < 4; k++){
                    if((b.pbases[k] + b.mbases[k]) > 0){
                        auto & bb = bouts[k];
                        bb.snps.push_back(snp_id);
                        bb.barcodes.push_back(b.barcode);
                        bb.plus.push_back(b.pbases[k]);
                        bb.minus.push_back(b.mbases[k]);
                    }
                }
            }
        }

        for(auto ptr : buffer){
            auto const & p = *ptr;
            os << genome_[p.tid].name << "\t" << p.pos << "\t" 
               << p.coverage << "\t" << p.barcodes << "\t" << p.ambig << "\t" 
               << p.ref << "\t" << p.max_nr << "\t" << p.pbase << "\t" << p.mbase
               << "\t" << p.sdists[0] << "\t" << p.sdists[1] << "\t" << p.sdists[3] << "\t" << p.sdists[2];

            for(size_t i = 0; i < 4; i++){
                auto const & b = p.bases[i];
                os << "\t" << (b.m_count + b.p_count) << "\t" << b.t_barcodes;
            }
            for(size_t i = 0; i < 4; i++){
                auto const & b = p.bases[i];
                os 
                    << "\t" << b.p_count
                    << "\t" << b.p_edge_dist;
            }
            for(size_t i = 0; i < 4; i++){
                auto const & b = p.bases[i];
                os 
                    << "\t" << b.m_count
                    << "\t" << b.m_edge_dist;
            }
            os << "\n";
        }

                /*
                std::cout
                    << "Piled up ref = " << p.ref << " " << " reads = " << (out_.end() - out_.begin())
                    << "  Position: " << p.tid << ": " << p.pos << " coverage = " << p.coverage << " ambig = " << p.ambig << "\n";
                for(size_t i = 0; i < 4; i++){
                    auto & b = p.bases[i];
                    if(b.count > 0){
                        std::cout << std::setprecision(2) << std::fixed
                            << "    base = " << "ACGT"[i] << " count = " << b.count 
                            << " gaps = " << (1.0 * b.gaps / b.count) << " NM = " << (1.0 * b.NM / b.count) 
                            << " barcodes = " << b.barcodes << " lens = " << (1.0 * b.lens / b.count) << "\n";
                    }
                }
                */

        if((reads - lreads) > 500000 && !buffer.empty()){
            size_t sec = tout.seconds();
            double ps = 1.0 * reads / (sec - start_time);
            tout << "Processed " << reads << " reads [" << ps << " / second], bases with min barcodes = " << bases << " plus = " << plus_bases << " minus = " << minus_bases << ", total passed bases = " << pbases << " current ref = " << genome_[buffer.back()->tid].name << ": " << buffer.back()->pos << "\n";
        }
    }

    tout << "Finished. Processed  " << reads << " reads, bases with min barcodes =  " << bases << " plus = " << plus_bases << " minus = " << minus_bases << ", total passed bases = " << pbases << "\n";

    delete rbuffer;
    delete pbuffer;

    std::vector<unsigned int> barcode_coverage(barcodes_.size());
    std::vector<unsigned int> barcode_bases(barcodes_.size());
    for(size_t i = 0; i < threads.size(); i++){
        std::transform (barcode_coverage.begin(), barcode_coverage.end(), threads[i]->bcoverage.begin(), barcode_coverage.begin(), std::plus<unsigned int>());
        std::transform (barcode_bases.begin(), barcode_bases.end(), threads[i]->bbases.begin(), barcode_bases.begin(), std::plus<unsigned int>());
        delete threads[i];
    }

    write_h5f_(out_ + "_barcode_matrices.h5", bouts);
    {
        gzofstream ofz(out_ + "_barcodes.txt.gz");
        ofz << "barcode\tmolecules\tbases_covered\tbases\n";
        for(size_t i = 0; i < barcodes_.size(); i++){
            ofz << barcodes_[i] << "\t" << barcode_molecules[i] << "\t" << barcode_coverage[i] << "\t" << barcode_bases[i] << "\n"; 

        }
    }
    bouts.clear();
    threads.clear();


    return EXIT_SUCCESS;
}

uint32_t count_bigger(const std::vector<uint32_t> & elems) {
    return std::count_if(elems.begin(), elems.end(), [](int c){return c > 0;});
}

void ProgPileup::write_h5f_(const std::string & out, std::vector<BBout> & bouts){
    using namespace H5;
    H5File file(out, H5F_ACC_TRUNC);
    //H5::Group group(file.createGroup("/barcode_rates"));
    std::vector<const char *> ctmp;
    for(auto & b : barcodes_) ctmp.push_back(b.c_str());
    write_h5_string("barcodes", ctmp, file);
    write_h5_numeric("refs", refs_, file, PredType::NATIVE_UINT8);
    write_h5_numeric("tids", tids_, file, PredType::NATIVE_UINT32);
    write_h5_numeric("pos", pos_, file, PredType::NATIVE_UINT32);
    for(size_t bi = 0; bi < 4; bi++){
        std::string gs = "/base_";
        auto & bb = bouts[bi];
        gs+="ACGT"[bi];
        H5::Group group(file.createGroup(gs));

        write_h5_numeric("snps", bb.snps, group, PredType::NATIVE_UINT32);
        write_h5_numeric("barcode_ids", bb.barcodes, group, PredType::NATIVE_UINT32);
        write_h5_numeric("plus", bb.plus, group, PredType::NATIVE_UINT16);
        write_h5_numeric("minus", bb.minus, group, PredType::NATIVE_UINT16);
    }

    {

        H5::Group group(file.createGroup("coverage"));
        std::vector<int32_t> ti;
        std::vector<uint32_t> tu;
        std::sort(coverage_.begin(), coverage_.end());
        ctmp.clear();

        for(size_t i = 0; i < genome_.size(); i++) ctmp.push_back(genome_[i].name.c_str());
        write_h5_string("chroms", ctmp, group);
        for(auto c : coverage_) ti.push_back(c.tid);
        write_h5_numeric("tid", ti, group, PredType::NATIVE_INT32);
        for(auto c : coverage_) tu.push_back(c.pos);
        write_h5_numeric("pos", tu, group, PredType::NATIVE_UINT32);

        tu.clear();
        for(auto c : coverage_) tu.push_back(c.pcoverage);
        std::cout << "plus > 0: " << count_bigger(tu) << " / " << tu.size() << "\n";
        write_h5_numeric("plus", tu, group, PredType::NATIVE_UINT32);


        tu.clear();
        for(auto c : coverage_) tu.push_back(c.mcoverage);
        std::cout << "minus > 0: " << count_bigger(tu) << " / " << tu.size() << "\n";
        write_h5_numeric("minus", tu, group, PredType::NATIVE_UINT32);

        tu.clear();
        for(auto c : coverage_) tu.push_back(c.tbarcodes);
        std::cout << "total_barcodes > 0: " << count_bigger(tu) << " / " << tu.size() << "\n";
        write_h5_numeric("total_barcodes", tu, group, PredType::NATIVE_UINT32);

        tu.clear();
        for(auto c : coverage_) tu.push_back(c.pbarcodes);
        std::cout << "plus_barcodes > 0: " << count_bigger(tu) << " / " << tu.size() << "\n";
        write_h5_numeric("plus_barcodes", tu, group, PredType::NATIVE_UINT32);

        tu.clear();
        for(auto c : coverage_) tu.push_back(c.mbarcodes);
        std::cout << "minus_barcodes > 0: " << count_bigger(tu) << " / " << tu.size() << "\n";
        write_h5_numeric("minus_barcodes", tu, group, PredType::NATIVE_UINT32);
    }

    file.close();
}

void ProgPileup::parse_targets_(){
    targets_.resize(genome_.size());
    std::map<std::string, int32_t> tids;
    for(size_t i = 0; i < genome_.size(); i++)
        tids[genome_[i].name] = i;

    FileWrapper in(snvlist_);
    ParserTokens toks;
    in.tokenize_line(toks);
    while(in.tokenize_line(toks) > -1){
        if(toks.empty() || toks[0].empty()) continue;
        auto it = tids.find(toks[0]);
        if(it == tids.end()){
            std::cout << "Missing reference " << toks[0] << "\n";
            continue;
        }
        uint32_t pos = std::stoul(toks[1]);
        targets_[it->second].push_back(TargetFinder::Target(pos, toks[3][0]));
    }
    size_t tot = 0;
    for(auto & t : targets_) {
        std::sort(t.begin(), t.end());
        tot += t.size();
    }
    tout << "Loaded " << tot << " pileup targets\n";
}

int ProgPileup::run() {
    if(cellranger_){
        if(lib_type_ == "V2"){
            return run_wrap_<Reader10X_V2, BamCellRangerProcessor>();
        }else if(lib_type_ == "V3"){
            return run_wrap_<Reader10X_V3, BamCellRangerProcessor>();
        }else if(lib_type_ == "V2_5P"){
            return run_wrap_<Reader10X_V2_5P, BamCellRangerProcessor>();
        }else if(lib_type_ == "V3_5P"){
            return run_wrap_<Reader10X_V3_5P, BamCellRangerProcessor>();
        }
    }else{
        if(lib_type_ == "V2"){
            return run_wrap_<Reader10X_V2, BamScSNVProcessor>();
        }else if(lib_type_ == "V3"){
            return run_wrap_<Reader10X_V3, BamScSNVProcessor>();
        }else if(lib_type_ == "V2_5P"){
            return run_wrap_<Reader10X_V2_5P, BamScSNVProcessor>();
        }else if(lib_type_ == "V3_5P"){
            return run_wrap_<Reader10X_V3_5P, BamScSNVProcessor>();
        }
    }
    return EXIT_SUCCESS;
}
