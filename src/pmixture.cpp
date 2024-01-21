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

#include "pmixture.hpp"
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

void ProgMixture::load() {
    index_ = args_["txidx"].as<std::string>();
    out_ = args_["out"].as<std::string>();
    ref_ = args_["reference"].as<std::string>();
    lib_type_ = args_["library"].as<std::string>("V2");
    threads_ = args_["threads"].as<unsigned int>(1);

    //min_coverage_ = args_["mincov"].as<unsigned int>(min_coverage_);
    min_alternative_ = args_["minalt"].as<unsigned int>(min_alternative_);
    min_qual_ = args_["minqual"].as<unsigned int>(min_qual_);
    min_barcodes_ = args_["minbarcodes"].as<unsigned int>(min_barcodes_);
    min_af_ = args_["minaf"].as<double>(min_af_);
    min_edge_ = args_["minedge"].as<unsigned int>(min_edge_);
    splice_win_ = args_["splicewin"].as<unsigned int>(splice_win_);
    mixture_ = args_["mixture"].as<std::string>();
    if(args_["seed"]){
        seed_ = args_["seed"].as<unsigned int>();
    }
    ambient_ = args_["ambient"].as<double>(ambient_);
    //cellranger_ = args_["cellranger"];

    dups_ = true;
    if(args_.pos.size() < 2){
        for(auto & a : args_.pos){
            std::cout << a << ", ";
        }
        std::cout << "\n";
        throw std::runtime_error("Requires at least two bam file arguments");
    }

    tout << "Loading the genome\n";
    {
        FastaReader fr(ref_);
        fr.read_all(genome_);
    }

    for(auto & p : args_.pos){
        bam_files_.push_back(p);
    }
    fbarcodes_.resize(bam_files_.size());
    bam_counts_.resize(bam_files_.size());
    molecules_.resize(bam_files_.size());
}

void tokenize_barcodes(std::string bcs, std::vector<std::pair<unsigned int, std::string>> & out){
    auto start = 0U;
    out.clear();
    while(start < bcs.size()){
        auto end = bcs.find(',', start);
        if(end == std::string::npos) end = bcs.length();
        std::string s = bcs.substr(start, end - start);
        auto i = s.find(':');
        unsigned int fno = std::stoi(s.substr(0, i));
        s = s.substr(i + 1);
        out.push_back({fno, s});
        start = end + 1;
    }
}

void ProgMixture::read_mixture_(){
    FileWrapper in(mixture_);
    ParserTokens toks;
    in.tokenize_line(toks);
    size_t index = 0;
    size_t singlets = 0, doublets = 0;
    std::vector<std::pair<unsigned int, std::string>> bcs;
    while(in.tokenize_line(toks) > -1){
        tokenize_barcodes(toks[1], bcs);
        barcodes_.push_back(toks[0]);
        //bmolecules_.push_back(std::stoi(toks[2]));
        for(auto & b : bcs){
            fbarcodes_[b.first][b.second] = index;
        }
        if(bcs.size() > 1){ 
            doublets++;
        }else{
            singlets++;
        }
        index++;
    }
    tout << "Read " << singlets << " singlets and " << doublets << " doublets\n";
    cell_molecules_.resize(barcodes_.size());
    for(size_t i = 0; i < bam_files_.size(); i++){
        std::cout << "Using " << fbarcodes_[i].size() << " barcodes from " << bam_files_[i] << "\n";
    }
}

template <typename T, typename P>
int ProgMixture::run_wrap_(){
    read_mixture_();
    tout << "Loading the transcriptome index\n";

    typename BamGeneReaderFiltered<T, BamMerger, P>::bcfilter_func fp = std::bind(&ProgMixture::filter_func, *this, std::placeholders::_1, std::placeholders::_2);
    BamGeneReaderFiltered<T, BamMerger, P> br(fp);
    /*
    if(collapse_){
        if(use_seed_){
            std::cout << "Setting random seed to " << seed_ << "\n";
            br.set_seed(seed_);
        }else{
            std::cout << "Using random seed of " << br.get_seed() << "\n";
        }
    }
    */

    //BamOutputBuffer cbuffer;
    //CollapseWorker cw(cbuffer, genome_);

    br.index.load(index_);
    br.index.build_splice_site_index();

    br.add_bams(bam_files_.begin(), bam_files_.end());
    br.prepare();

    BamBuffer * rbuffer = new BamBuffer();
    BamBuffer * pbuffer = new BamBuffer();
    std::vector<PileupWorker*> threads;
    if(threads_ < 2) threads_ = 2;
    for(size_t i = 0; i < threads_ - 1; i++){
        threads.push_back(new PileupWorker(barcodes_.size(), br.index, genome_, dups_));
        threads.back()->set_params(min_alternative_, min_qual_, min_barcodes_, min_edge_, splice_win_, min_af_);
    }

    unsigned int pbases = 0, bases = 0, reads = 0, areads = 0, atotal = 0;
    unsigned int plus_bases = 0, minus_bases = 0;
    unsigned int lreads = 0;
    std::vector<PositionCount*> buffer;
    //std::vector<PositionCoverage*> cbuffer;
    std::vector<BarcodeCount*> bcbuffer;
    gzofstream os(out_ + ".txt.gz");
    std::vector<BBout> bouts;
    bouts.resize(4);

    //cos << "chrom\tpos\tplus_coverage\tminus_coverage\n";
    //os << "chrom\tpos\tcoverage\tbarcodes\tfailed_count\tref\tmax_non_ref\tplus_base\tminus_base";
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
    unsigned int tot = br.read_genes(*rbuffer, 250);
    tot = 0;
    std::mt19937 gen(seed_);
    std::uniform_int_distribution<unsigned int> dist(0, barcodes_.size() - 1); //barcodes_.begin(), barcodes_.end());
    std::uniform_real_distribution<double> pdist(0, 1);
    {
        auto ret = rbuffer->get_all_unsafe();
        //std::cout << "rbuff count = " << (ret.second - ret.first) << "\n";
        for(auto it = ret.first; it != ret.second; it++){
            auto gt = (*it)->filenum;
            if(ambient_ > 0.0 && pdist(gen) < ambient_){
                unsigned int nb = dist(gen);
                while(nb == (*it)->barcode){
                    nb = dist(gen);
                }
                cell_molecules_[(*it)->barcode].lost_molecules++;
                (*it)->barcode = nb;
                cell_molecules_[(*it)->barcode].ambient_molecules++;
                cell_molecules_[(*it)->barcode].gt_ambients[gt]++;
                areads++;
            }else{
                cell_molecules_[(*it)->barcode].kept_molecules++;
                cell_molecules_[(*it)->barcode].gt_totals[gt]++;
            }
            atotal++;
            //std::cout << " " << debug_read(*(*it)) << "\n";
            bam_counts_[(*it)->filenum]++;
            tot++;
        }
    }

    while(tot > 0){
        std::swap(rbuffer, pbuffer);
        for(auto it = threads.begin(); it != threads.end(); it++){
            auto & t = *(*it);
            t.set_buffer(pbuffer);
            t.start();
        }

        // While the threads are processing the buffer read in the next set to speed things up
        //std::cout << "pre read pbuffer = " << pbuffer->count() << " rbuffer = " << rbuffer->count() << "\n";
        tot = br.read_genes(*rbuffer, 250);


        /*
        if(collapse_){
            process_collapsed(rbuffer, cw, downsamples_, cstarts, molecules_);
        }
        */


        // Can probably collapse the reads here

        //std::cout << "post read pbuffer = " << pbuffer->count() << " rbuffer = " << rbuffer->count() << "\n";
        tot = 0;
        {
            auto ret = rbuffer->get_all_unsafe();
            for(auto it = ret.first; it != ret.second; it++){
                auto gt = (*it)->filenum;
                if(ambient_ > 0.0 && pdist(gen) < ambient_){
                    unsigned int nb = dist(gen);
                    while(nb == (*it)->barcode){
                        nb = dist(gen);
                    }
                    cell_molecules_[(*it)->barcode].lost_molecules++;
                    (*it)->barcode = nb;
                    cell_molecules_[(*it)->barcode].ambient_molecules++;
                    cell_molecules_[(*it)->barcode].gt_ambients[gt]++;
                    areads++;
                }else{
                    cell_molecules_[(*it)->barcode].kept_molecules++;
                    cell_molecules_[(*it)->barcode].gt_totals[gt]++;
                }
                atotal++;
                //std::cout << " " << debug_read(*(*it)) << "\n";
                bam_counts_[(*it)->filenum]++;
                tot++;
            }
        }

        pbases = 0; bases = 0; reads = 0;
        plus_bases = 0, minus_bases = 0;
        for(auto it = threads.begin(); it != threads.end(); it++){
            (*it)->join();
        }

        buffer.clear();
        //cbuffer.clear();
        bcbuffer.clear();
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
               << "\t" << p.sdists[0] << "\t" << p.sdists[1] << "\t" << p.sdists[2] << "\t" << p.sdists[3];

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

        if((reads - lreads) > 500000 && !buffer.empty()){
            tout << "Processed " << reads << " reads";
            if(areads > 0){
                std::cout << " swapped bacodes for " << areads << "  / " << atotal << " [" << std::fixed << std::setprecision(2) << (100.0 * areads / atotal) << "]";
            }
            std::cout << " current ref = " << genome_[buffer.back()->tid].name << ": " << buffer.back()->pos << "\n";
        }
    }

    {
        tout << "Finished. Processed  " << reads << " reads, bases with min barcodes =  " << bases << " plus = " << plus_bases << " minus = " << minus_bases << ", total passed bases = " << pbases << "\n";
        if(areads > 0){
            tout << "Swapped bacodes for " << areads << " [" << std::fixed << std::setprecision(2) << (100.0 * areads / reads) << "]";
        }
        std::cout << "\nFile Summary:\n";
        std::vector<unsigned int> ftot = br.file_counts();

        std::ofstream sout(out_ + "_mix_counts.txt");
        /*
        if(collapse_){
            sout << "bam_file\tmolecules\treads_used\ttotal\n";
            for(size_t i = 0; i < bam_files_.size(); i++){
                size_t rtot = br.ds_kept[i] + br.ds_skipped[i];
                std::cout << "  " << i << ": " << bam_files_[i] << " collapsed " << bam_counts_[i] << " molecules from " << br.ds_kept[i] << " reads out of " 
                    << rtot << " total reads matching barcods"
                    << std::fixed << std::setprecision(2) << " [" << std::fixed << (100.0 * br.ds_kept[i] / rtot) << "%]"
                    << " saturation = " << std::fixed << (100.0 - 100.0 * bam_counts_[i] / br.ds_kept[i]) << "%" 
                    << "\n";
                sout << bam_files_[i] << "\t" << bam_counts_[i] << "\t" << br.ds_kept[i] << "\t" << rtot << "\n";
            }
        }else{
        */
            sout << "bam_file\tused\ttotal\n";
            for(size_t i = 0; i < bam_files_.size(); i++){
                std::cout << "  " << i << ": " << bam_files_[i] << " used " << bam_counts_[i] << " reads out of " 
                    << std::fixed << std::setprecision(2) << ftot[i] << " [" << std::fixed << (100.0 * bam_counts_[i] / ftot[i]) << "%]\n";
                sout << bam_files_[i] << "\t" << bam_counts_[i] << "\t" << ftot[i] << "\n";
        //    }
        }
    }

    //BarcodeRate total_rates;
    //barcode_rates_.resize(barcodes_.size());

    std::vector<unsigned int> barcode_coverage(barcodes_.size());
    std::vector<unsigned int> barcode_bases(barcodes_.size());
    for(size_t i = 0; i < threads.size(); i++){
        std::transform (barcode_coverage.begin(), barcode_coverage.end(), threads[i]->bcoverage.begin(), barcode_coverage.begin(), std::plus<unsigned int>());
        std::transform (barcode_bases.begin(), barcode_bases.end(), threads[i]->bbases.begin(), barcode_bases.begin(), std::plus<unsigned int>());
        delete threads[i];
    }

    {
        gzofstream ofz(out_ + "_barcodes.txt.gz");
        ofz << "barcode\tkept_molecules\tambient_molecules\tlost_molecules\tbarcode_bases\tgt1_molecules\tgt2_molecules\tgt1_ambient\tgt2_ambient\n";
        for(size_t i = 0; i < barcodes_.size(); i++){
            ofz << barcodes_[i] << "\t" << cell_molecules_[i].kept_molecules << "\t" << cell_molecules_[i].ambient_molecules 
                << "\t" << cell_molecules_[i].lost_molecules << "\t" << barcode_bases[i] 
                << "\t" << cell_molecules_[i].gt_totals[0] << "\t" << cell_molecules_[i].gt_totals[1] 
                << "\t" << cell_molecules_[i].gt_ambients[0] << "\t" << cell_molecules_[i].gt_ambients[1]
                << "\n"; 
        }
    }

    delete rbuffer;
    delete pbuffer;


    /*
    {
        std::cout << "Base Rates (with at least " << min_barcodes_ << " overlapping barcodes)\n";
        std::ofstream out(out_ + "_base_rates.txt");
        out << "base\tcode";
        for(size_t i = 0; i < 4; i++) out << "\t" << "ACGT"[i] << "_plus\t" << "ACGT"[i] << "_minus";
        out << "\n";

        std::array<std::string, 6> bnames{{"coding", "non_coding", "utr3", "utr5", "intron", "intergenic"}};
        for(size_t ri = 0; ri < 6; ri++){
            //std::cout << "    " << bnames[ri] << " plus = " << total_rates.prates[ri] << " minus = " << total_rates.mrates[ri] << "\n";
            std::cout << bnames[ri];
            out << bnames[ri] << "\t" << "CE35NI"[ri];
            for(size_t i = 0; i < 4; i++){
                std::cout << "\t" << "ACGT"[i] << " = " << total_rates.prates[6 * i + ri] << "/" << total_rates.mrates[6 * i + ri];
                out << "\t" << total_rates.prates[6 * i + ri] << "\t" << total_rates.mrates[6 * i + ri];
            }
            out << "\n";
            std::cout << "\n";
        } 
    }
    */

    write_h5f_(out_ + "_barcode_matrices.h5", bouts);
    bouts.clear();
    threads.clear();

    return EXIT_SUCCESS;
}

void ProgMixture::write_h5f_(const std::string & out, std::vector<BBout> & bouts){
    using namespace H5;
    H5File file(out, H5F_ACC_TRUNC);
    //H5::Group group(file.createGroup("/barcode_rates"));
    std::vector<const char *> ctmp;
    for(auto & b : barcodes_) ctmp.push_back(b.c_str());
    write_h5_string("barcodes", ctmp, file);
    ctmp.clear();
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


    /*
    {
        H5::Group group(file.createGroup("barcode_rates"));
        std::vector<uint32_t> tr;
        std::array<std::string, 6> bnames{{"coding", "non_coding", "utr3", "utr5", "intron", "intergenic"}};
        for(size_t bi = 0; bi < 6; bi++){
            for(size_t i = 0; i < 4; i++){
                tr.clear();
                for(auto const & r : barcode_rates_){
                    tr.push_back(r.prates[bi + i * 6]);
                }
                std::string k("plus_");
                k += "ACGT"[i];
                write_h5_numeric(k + "_" + bnames[bi], tr, group, PredType::NATIVE_UINT32);
            }
        }
        for(size_t bi = 0; bi < 6; bi++){
            for(size_t i = 0; i < 4; i++){
                tr.clear();
                for(auto const & r : barcode_rates_){
                    tr.push_back(r.prates[bi + i * 6]);
                }
                std::string k("minus_");
                k += "ACGT"[i];
                write_h5_numeric(k + "_" + bnames[bi], tr, group, PredType::NATIVE_UINT32);
            }
        }
    }
    */

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
        write_h5_numeric("pos", tu, group, PredType::NATIVE_INT32);
        tu.clear();
        for(auto c : coverage_) tu.push_back(c.pcoverage);
        write_h5_numeric("plus", tu, group, PredType::NATIVE_INT32);
        tu.clear();
        for(auto c : coverage_) tu.push_back(c.mcoverage);
        write_h5_numeric("minus", tu, group, PredType::NATIVE_INT32);
        //genome_[p.tid].name
    }

    file.close();
}

int ProgMixture::run() {
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
