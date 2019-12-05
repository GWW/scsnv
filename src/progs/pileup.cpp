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

#include "pileup.hpp"
#include "../scmap/reader.hpp"
#include "../scmap/align_aux.hpp"
#include "../util/misc.hpp"
#include "../util/gzstream.hpp"
#include "../muts/pileup_worker.hpp"
#include "../util/h5misc.hpp"
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

    //min_coverage_ = args_["mincov"].as<unsigned int>(min_coverage_);
    min_alternative_ = args_["minalt"].as<unsigned int>(min_alternative_);
    min_qual_ = args_["minqual"].as<unsigned int>(min_qual_);
    min_barcodes_ = args_["minbarcodes"].as<unsigned int>(min_barcodes_);
    min_af_ = args_["minaf"].as<double>(min_af_);

    bc_counts_ = args_["barcodes"].as<std::string>();
    cellranger_ = args_["cellranger"];
    dups_ = args_["dups"];
    if(args_.pos.size() != 1){
        throw std::runtime_error("Missing the bam file argument");
    }

    tout << "Loading the genome\n";
    {
        FastaReader fr(ref_);
        fr.read_all(genome_);
    }
    //tout << "Skipping genome load\n";


    bam_files_.push_back(args_.pos[0]);
}

template <typename T>
int ProgPileup::run_wrap_(){
    std::vector<UMIMap> umi_map;
    tout << "Loading the transcriptome index\n";

    BamGeneReader<T> br(false); //br(cellranger_);
    typename T::LibraryBarcode bc;
    bc.load(bc_counts_);

    br.add_bams(bam_files_.begin(), bam_files_.end());
    br.index.load(index_);
    br.prepare();

    BamBuffer * rbuffer = new BamBuffer();
    BamBuffer * pbuffer = new BamBuffer();
    PileupWorker::cb_bhash cbhash = std::bind(&T::LibraryBarcode::bid, &bc, std::placeholders::_1);
    std::vector<PileupWorker*> threads;
    if(threads_ < 2) threads_ = 2;
    for(size_t i = 0; i < threads_ - 1; i++){
        threads.push_back(new PileupWorker(br.index, genome_, cbhash, dups_, cellranger_));
        threads.back()->set_params(min_alternative_, min_qual_, min_barcodes_, min_af_);
    }

    unsigned int pbases = 0, bases = 0, reads = 0;
    unsigned int plus_bases = 0, minus_bases = 0;
    unsigned int lreads = 0;
    std::vector<PositionCount*> buffer;
    //std::vector<PositionCoverage*> cbuffer;
    std::vector<BarcodeCount*> bcbuffer;
    gzofstream os(out_ + ".txt.gz");
    std::vector<BBout> bouts;
    bouts.resize(4);

    //cos << "chrom\tpos\tplus_coverage\tminus_coverage\n";
    os << "chrom\tpos\tcoverage\tbarcodes\tfailed_count\tref\tplus_base\tminus_base";
    for(size_t i = 0; i < 4; i++){
        char bs = "ACGT"[i];
        os << "\t" 
            << bs << "_total\t" 
            << bs << "_total_barcodes";
    }

    for(size_t i = 0; i < 4; i++){
        std::string bs(1, "ACGT"[i]);
        bs += "_plus";
        os << "\t" 
            << bs << "_counts" ;
            //<< bs << "_gaps\t" 
            //<< bs << "_mismatches\t" 
            //<< bs << "_tot_lens";
    }
    for(size_t i = 0; i < 4; i++){
        std::string bs(1, "ACGT"[i]);
        bs += "_minus";
        os << "\t" 
            << bs << "_counts" ;
            //<< bs << "_gaps\t" 
            //<< bs << "_mismatches\t" 
            //<< bs << "_tot_lens";
    }
    os << "\n";
    unsigned int tot = br.read_genes(*rbuffer, 250);

    while(tot > 0){
        std::swap(rbuffer, pbuffer);
        for(auto it = threads.begin(); it != threads.end(); it++){
            auto & t = *(*it);
            t.set_buffer(pbuffer);
            t.start();
        }

        // While the threads are processing the buffer read in the next set to speed things up
        tot = br.read_genes(*rbuffer, 250);

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
        //std::sort(cbuffer.begin(), cbuffer.end(), [](const PositionCoverage * p1, PositionCoverage * p2) { return (*p1) < (*p2); });
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

        /*
        for(auto ptr : cbuffer){
            auto const & p = *ptr;
            cos << genome_[p.tid].name << "\t" << p.pos << "\t" 
               << p.pcoverage << "\t" << p.mcoverage << "\n";
        }
        */
        for(auto ptr : buffer){
            auto const & p = *ptr;
            os << genome_[p.tid].name << "\t" << p.pos << "\t" 
               << p.coverage << "\t" << p.barcodes << "\t" << p.ambig << "\t" << p.ref << "\t" << p.pbase << "\t" << p.mbase;
            for(size_t i = 0; i < 4; i++){
                auto const & b = p.bases[i];
                os << "\t" << (b.m_count + b.p_count) << "\t" << b.t_barcodes;
            }
            for(size_t i = 0; i < 4; i++){
                auto const & b = p.bases[i];
                os 
                    << "\t" << b.p_count;
                    //<< "\t" << b.p_gaps
                    //<< "\t" << b.p_NM
                    //<< "\t" << b.p_lens;
            }
            for(size_t i = 0; i < 4; i++){
                auto const & b = p.bases[i];
                os 
                    << "\t" << b.m_count;
                    //<< "\t" << b.m_gaps
                    //<< "\t" << b.m_NM
                    //<< "\t" << b.m_lens;
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
            tout << "Processed " << reads << " reads, bases with min barcodes = " << bases << " plus = " << plus_bases << " minus = " << minus_bases << ", total passed bases = " << pbases << " current ref = " << genome_[buffer.back()->tid].name << ": " << buffer.back()->pos << "\n";

        }
    }

    tout << "Finished. Processed  " << reads << " reads, bases with min barcodes =  " << bases << " plus = " << plus_bases << " minus = " << minus_bases << ", total passed bases = " << pbases << " current ref = " << genome_[buffer.back()->tid].name << ": " << buffer.back()->pos << "\n";

    delete rbuffer;
    delete pbuffer;
    BarcodeRate total_rates;

    barcode_rates_.resize(bc.size());

    for(size_t i = 0; i < threads.size(); i++){
        total_rates += threads[i]->total_rates;
        for(size_t j = 0; j < threads[i]->barcode_rates.size(); j++){
            barcode_rates_[j] += threads[i]->barcode_rates[j];
        }
        delete threads[i];
    }


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

    write_h5f_(out_ + "_barcode_matrices.h5", bc, bouts);
    bouts.clear();
    threads.clear();

    return EXIT_SUCCESS;
}

template <typename TL>
void ProgPileup::write_h5f_(const std::string & out, TL & bc, std::vector<BBout> & bouts){
    using namespace H5;
    H5File file(out, H5F_ACC_TRUNC);
    //H5::Group group(file.createGroup("/barcode_rates"));
    std::vector<const char *> ctmp;
    for(size_t i = 0; i < bc.size(); i++) ctmp.push_back(bc.barcode(i).c_str());
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

int ProgPileup::run() {
    if(lib_type_ == "V2"){
        return run_wrap_<Reader10X_V2>();
    }else if(lib_type_ == "V3"){
        return run_wrap_<Reader10X_V3>();
    }
    return EXIT_SUCCESS;
}
