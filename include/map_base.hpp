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
#include "reader.hpp"
#include "barcodes.hpp"
#include "map_worker.hpp"
#include "gzstream.hpp"
#include "pmap.hpp"
#include "sparsepp/sparsepp/spp.h"
#include "sbam_writer.hpp"
#include "transcript_align.hpp"
#include "genome_align.hpp"
#include <exception>
#include <fstream>
#include <list>
#include <functional>
#include <map>
#include <unistd.h>

namespace gwsc{

template <typename T>
class MapBase {
    public:
        MapBase() : counts{}, index_(), bout_(), txa_(index_, T::LibraryStrand), gna_(index_, T::LibraryStrand) {
            smode_ = T::LibraryStrand;
        }

        void load_index(const std::string & txindex, unsigned int min_overhang, const std::string & gindex){
            index_.load(txindex);
            index_.build_splice_index(min_overhang);
            txa_.load(txindex, min_overhang);
            AlignScore ascore;
            gna_.ascore = ascore;
            txa_.ascore = ascore;
            gna_.load(gindex, min_overhang);
            txa_.set_gidx(&gna_);
            //index_.load(txindex, exonic_intergenic_ratio, true);
            //if(!gindex.empty()) index_.load_genome(gindex);
        }

        void unload(){
            gna_.unload();
            txa_.unload();
        }

        size_t total_reads() const {
            return btotal_;
        }

        void load_barcode_counts(const std::string & barcodes, FastqPairs & fastqs){
            bc_.load(barcodes + "_counts.txt.gz");
            btotal_ = T::LibraryBarcode::find_total_reads(barcodes, fastqs);
        }

        void get_barcodes(std::vector<std::string> & bnames) const {
            for(size_t i = 0; i < bc_.size(); i++){
                bnames.push_back(bc_.barcode(i));
            }
        }

        void prepare_bam(const std::string & prog_name, unsigned int bam_per_thread, unsigned int bam_per_file, 
                const std::string & prefix, unsigned int threads);

        void run(unsigned int num_threads, FastqPairs & fastqs, double dust, bool internal, size_t downsample, size_t seed);
        void write_output(const std::string & prefix);
        void write_tags(const std::string & prefix, const std::vector<UMIMap> & correct, const std::vector<uint32_t> & bidx);

        TXIndex & index() {
            return index_;
        }

        htsThreadPool & pool() {
            return bout_.pool();
        }

        AlignGroup::ResultCounts   counts;
        std::vector<AlignSummary>  aligns;
        std::map<AlignSummary::bint, AlignGroup::ResultCounts> brates;
        unsigned int               barcode_correct = 0;
        unsigned int               barcode_corrected = 0;

    private:
        unsigned int read_(size_t N, Reads & reads, const AlignGroup::ResultCounts & counts);

        using lib_bc = typename T::LibraryBarcode;

        MultiReader<T>             in_;
        lib_bc                     bc_;
        TXIndex                    index_;
        SortedBamWriter            bout_;
        TranscriptAlign            txa_;
        GenomeAlign                gna_;
        std::mutex                 mtx_read_;
        std::string                bam_tmp_;
        double                     ds_ = 0.0;
        size_t                     start_;
        size_t                     total_ = 0;
        size_t                     ltotal_ = 0;
        size_t                     rtotal_ = 0;
        size_t                     btotal_ = 0;
        int                        write_threads_ = 0;
        StrandMode                 smode_;
        bool                       bam_ = false;
        bool                       internal_ = false;
};

template <typename T>
inline void MapBase<T>::run(unsigned int num_threads, FastqPairs & fastqs, double dust, bool internal, size_t downsample, size_t seed) {
    static const unsigned int READS_PER_STEP = 2000;
    internal_ = internal;
    start_ = tout.seconds();
    std::list<MapWorker> threads;
    in_.set_files(fastqs);


    if(downsample > 0){
        double perc = 1.0 * downsample / btotal_;
        ds_ = perc;
        in_.set_downsample(perc, seed);
        tout << "Downsampling to " << downsample << " reads [" << std::fixed << std::setprecision(2) << (100.0 * perc) << "% of total reads]\n";
        tout << "Estimated memory needed for alignment tags: " << std::setprecision(2) << std::fixed
            << (1.0 * sizeof(AlignSummary) * downsample  / (1024 * 1024 * 1024)) << " GB\n";
    }else{
        tout << "Estimated memory needed for alignment tags: " << std::setprecision(2) << std::fixed
            << (1.0 * sizeof(AlignSummary) * btotal_ / (1024 * 1024 * 1024)) << " GB\n";
    }

    using namespace std::placeholders;
    MapWorker::cb_correct correct_cb = std::bind(&lib_bc::correct, &bc_, _1, _2);
    MapWorker::cb_read read_cb = std::bind(&MapBase<T>::read_, this, _1, _2, _3);

    for(size_t i = 0; i < num_threads; i++){
        threads.emplace_back(READS_PER_STEP, read_cb, correct_cb, smode_, dust);
        threads.back().tx_align = &txa_;
        threads.back().genome_align = &gna_;
        threads.back().tx_idx = &index_;
        if(bam_) {
            threads.back().bout = &bout_;
            threads.back().prepare_bam();
        }
    }
    auto it = threads.begin();
    while(++it != threads.end()){
        it->start();
    }
    threads.front()();
    it = threads.begin();
    while(++it != threads.end()){
        it->join();
    }

    it = threads.begin();
    size_t test = 0;
    while(it != threads.end()){
        if(bam_) {
            bout_.merge_buffer(it->buff, it->rcount);
        }
        barcode_corrected += it->barcode_corrected;
        barcode_correct += it->barcode_correct;
        for(auto & g : it->aligns){
            aligns.push_back(g);
        }
        for(auto & r : it->bc_rates){
            for(size_t i = 0; i < r.second.size(); i++) {
                brates[r.first][i] += r.second[i];
                test += r.second[i];
            }
        }
        it = threads.erase(it);
    }
    if(bam_) bout_.force_write();
    tout << "Merged alignments size = " << aligns.size() << " Merged barcode rate size = " << brates.size() << " sum = " << test << "\n";
    tout << "Wrote " << bout_.total_reads() << " sorted alignments across " << bout_.total_files() << " bam files\n";
    tout << "Sorting the alignment tags\n";
    std::sort(aligns.begin(), aligns.end());
}

template <typename T>
inline void MapBase<T>::write_output(const std::string & prefix) {
{
        std::ofstream out(prefix + "alignment_summary.txt");
        out << "type\tcount\tpercent\n";
        std::cout << "Alignment summary total reads = " << btotal_ << "\n";
        for(size_t i = 0; i < AlignGroup::ELEM_COUNT; i++){
            std::cout << std::setfill(' ') << std::setw(35) << std::left <<
                AlignGroup::alignres2str(static_cast<AlignGroup::Result>(i), true)  << std::right
                << "\t" << std::setfill(' ') << std::setw(9) << counts[i] << "\t" 
                << std::fixed << std::setprecision(2) << std::setw(6) << 100.0 * counts[i] / total_ << "%\n";
            out << AlignGroup::alignres2str(static_cast<AlignGroup::Result>(i), false) 
                << "\t" << counts[i] << "\t" 
                << std::fixed << std::setprecision(2) << 100.0 * counts[i] / total_ << "%\n";
        }
        std::cout << "\n";
        std::cout << std::setfill(' ') << std::setw(35) << std::left <<
            "Barcodes Correct" << std::right
            << "\t" << std::setfill(' ') << std::setw(9) << barcode_correct << "\t" 
            << std::fixed << std::setprecision(2) << std::setw(6) << 100.0 * barcode_correct / total_ << "%\n";
        std::cout << std::setfill(' ') << std::setw(35) << std::left <<
            "Barcodes Corrected" << std::right
            << "\t" << std::setfill(' ') << std::setw(9) << barcode_corrected << "\t" 
            << std::fixed << std::setprecision(2) << std::setw(6) << 100.0 * barcode_corrected / total_ << "%\n";
        out << "barcodes_correct" 
            << "\t" << barcode_correct << "\t" 
            << std::fixed << std::setprecision(2) << 100.0 * barcode_correct / total_ << "%\n";
        out << "barcodes_corrected" 
            << "\t" << barcode_corrected << "\t" 
            << std::fixed << std::setprecision(2) << 100.0 * barcode_corrected / total_ << "%\n";
    }

    /*

    {
        gzofstream os(prefix + "_barcode_rates.txt.gz");
        os << "barcode_id\tbarcode\tbarcode_corrected";
        for(size_t i = 1; i < AlignGroup::ELEM_COUNT; i++)
            os << "\t" << AlignGroup::alignres2str(static_cast<AlignGroup::Result>(i), false);

        os << "\n";
        for(auto const & p : brates){
            os << p.first << "\t" << bc_.barcode(p.first);
            for(size_t i = 0; i < AlignGroup::ELEM_COUNT; i++)
                os << "\t" << p.second[i];
            os << "\n";
        }

    }

    */
}
template <typename T>
inline void MapBase<T>::write_tags(const std::string & prefix, const std::vector<UMIMap> & correct, const std::vector<uint32_t> & bidx) {
    gzFile zout = gzopen((prefix + "_tags.gz").c_str(), "wb");
    gzFile zout2 = gzopen((prefix + "_tags_idx.gz").c_str(), "wb");

    tout << "Correcting the tag data\n";
    //std::ofstream os(prefix + "_test.txt");
    std::vector<bool> ccheck(correct.size());
    auto rit = correct.begin() + bidx[0];
    auto rend = correct.begin() + bidx[1];
    auto bstart = aligns.begin();
    bool corrected = false;
    for(auto it = aligns.begin(); it != aligns.end(); it++){
        if(bstart->barcode != it->barcode){
            if(corrected){
                std::sort(bstart, it,
                    [&](const AlignSummary & m1, AlignSummary & m2) { 
                    return 
                        std::tie(m1.gene_id, m1.umi) < 
                        std::tie(m2.gene_id, m2.umi);
                    }
                );
            }

            rit = correct.begin() + bidx[it->barcode];
            rend = correct.begin() + bidx[it->barcode + 1];
            bstart = it;
            corrected = false;
        }

        while(rit != rend && std::tie(rit->gene_id, rit->umi_from) < std::tie(it->gene_id, it->umi)){
            rit++;
        }
        if(rit != rend && rit->gene_id == it->gene_id && rit->umi_from == it->umi){
            size_t ridx = std::distance(correct.begin(), rit);
            ccheck[ridx] = true;
            corrected = true;
            it->umi = rit->umi_to;
        }
    }
    size_t missing = 0;
    for(auto c : ccheck){
        if(!c){
            missing++;
        }
    }

    if(missing > 0)
        std::cout << "Missed " << missing << " out of " << correct.size() << " UMI corrections\n";
    tout << "Writing the tag data\n";
    AlignSummary::bint lbarcode = std::numeric_limits<AlignSummary::bint>::max();

    size_t i = 0;
    for(auto & a : aligns){
        if(a.barcode != lbarcode){
            lbarcode = a.barcode;
            gzwrite(zout2, reinterpret_cast<char*>(&lbarcode), sizeof(lbarcode));
            gzwrite(zout2, reinterpret_cast<char*>(&i), sizeof(i));
        }

        AlignTagOut tag(a.barcode, a.gene_id, a.umi, a.intronic);
        gzwrite(zout, reinterpret_cast<char*>(&tag), sizeof(AlignTagOut));
        i++;
    }
    lbarcode = std::numeric_limits<AlignSummary::bint>::max();
    gzwrite(zout2, reinterpret_cast<char*>(&lbarcode), sizeof(lbarcode));
    gzwrite(zout2, reinterpret_cast<char*>(&i), sizeof(i));
    gzclose(zout);
    gzclose(zout2);
}

template <typename T>
unsigned int MapBase<T>::read_(size_t N, Reads & reads, const AlignGroup::ResultCounts & rcounts) {
    std::lock_guard<std::mutex> lock(mtx_read_);
    // Update alignment counts
    for(size_t i = 0; i < rcounts.size(); i++) {
        counts[i] += rcounts[i];
        total_ += rcounts[i];
    }
    if((ltotal_ + 2500000) <= total_){
        size_t sec = tout.seconds();
        double ps = 1.0 * total_ / (sec - start_);
        size_t eta = (btotal_ - total_) / ps;
        int hours = eta / (60 * 60);
        int minutes = (eta - (hours * 60 * 60)) / 60;

        double punque = 100.0 * counts[AlignGroup::CDNA] / total_;
        double piunque = 100.0 * counts[AlignGroup::INTRONIC] / total_;
        double pinter = 100.0 * counts[AlignGroup::INTERGENIC] / total_;
        double pambig = 100.0 * counts[AlignGroup::AMBIGUOUS] / total_;
        double punmapped = 100.0 * (counts[AlignGroup::UNMAPPED] + counts[AlignGroup::UMI_FAIL] + counts[AlignGroup::BARCODE_FAIL]) / total_;
        double pantisense = 100.0 * counts[AlignGroup::ANTISENSE] / total_;
        double pmulti = 100.0 * counts[AlignGroup::MULTIMAPPED] / total_;
        double pfail = 100.0 * counts[AlignGroup::TAG_FAIL] / total_;
        tout << "Processed " << total_ << " / " << btotal_ << " [" << static_cast<int>(ps) << " / sec], ETA = " 
            << hours << "h " << minutes << "m "
            << " CDNA: " << std::setprecision(2) << punque
            << " PRE: " << std::setprecision(2) << piunque
            << " INT: " << std::setprecision(2) << pinter
            << " AMB: " << std::setprecision(2) << pambig
            << " ANT: " << std::setprecision(2) << pantisense
            << " MUL: " << std::setprecision(2) << pmulti
            << " UNM: " << std::setprecision(2) << punmapped
            << " QA: " << std::setprecision(2) << pfail;
        if(ds_ > 0.0){
            double ps = 100.0 * in_.skipped() / (in_.skipped() + total_);
            std::cout << " DS SKIPPED: " << std::setprecision(2) << ps << " GOAL = " << (100.0 * ds_);

        }
        std::cout << "\n";
        ltotal_ = total_;
    }
    //if(total_ >= 100000) return 0;
    auto read = in_.read_N(N, reads);
    rtotal_ += read;
    return read;
}

template <typename T>
void MapBase<T>::prepare_bam(const std::string & prog_name, unsigned int bam_per_thread, unsigned int bam_per_file, 
        const std::string & prefix, unsigned int threads)
{
    tout << "Preparing bam header for " << index_.refs().size() << " references\n";
    bam_tmp_ = prefix;
    bam_ = true;
    auto & bh = bout_.bh;
    for(auto & r : index_.refs()){
        bh.add_name(r.name, r.len, r.tid);
    }
    bh.build_header();
    std::stringstream lines;
    lines << "@HD\tSO:coordinate\n@PG\tID:scsnv\tCL:" << prog_name << "\tVN:1.0\n";
    bh.set_text(lines.str());
    bout_.set_thread_buffer(bam_per_file, bam_per_thread);
    bout_.set_prefix(prefix);
    if(threads > 1) bout_.make_pool(threads);
    write_threads_ = threads;
}

}
