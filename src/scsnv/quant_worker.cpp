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

#include "quant_worker.hpp"
#include "../util/gzstream.hpp"
#include "../util/h5misc.hpp"
#include <algorithm>
#include <cstdio>
#include <H5Cpp.h>
#include <numeric>

using namespace gwsc;


void QuantBase::prepare_tags(const std::vector<std::string> & prefixes) {
    prefixes_ = prefixes;
    for(auto & p : prefixes){
        in_.add(p);
    }
}

void QuantBase::build_gene_groups(const std::string & fname) {
    if(fname.empty()) return;
    ParserTokens toks;
    FileWrapper in(fname);
    in.tokenize_line(toks);
    std::unordered_map<std::string, uint32_t> gid_map;
    std::unordered_map<std::string, uint32_t> group_map;
    for(size_t i = 0; i < index_.max_gid(); i++){
        gid_map[index_.gene(i).gene_id] = i;
    }
    size_t skipped = 0;
    while(in.tokenize_line(toks) >= 0){
        auto git = gid_map.find(toks[0]);
        if(git == gid_map.end()){
            //std::cout << "Skipping gene_id " << toks[0] << " from group " << toks[1] << "\n";
            skipped++;
            continue;
        }
        auto grit = group_map.insert(std::make_pair(toks[1], group_map.size()));
        group_map_[git->second] = grit.first->second;
    }

    group_names_.resize(group_map.size());
    std::cout << "Skipping " << skipped << " genes from the group index\n";
    for(auto & g : group_map){
        std::cout << g.second << " - " << g.first << "\n";
        group_names_[g.second] = g.first;
    }
}

// Sort barcodes by their index
//From https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes
template <typename T>
std::vector<size_t> sort_molecule_indexes(const std::vector<T> &v, size_t NG) {
    // initialize original index locations
    std::vector<size_t> idx(v.size() / NG);
    std::iota(idx.begin(), idx.end(), 0);

    sort(idx.begin(), idx.end(),
            [&v, &NG](size_t i1, size_t i2) {return v[i1 * NG] > v[i2 * NG];});

    return idx;
}

void QuantBase::run(unsigned int umi_len, unsigned int num_threads, unsigned int min_molecules, bool bam, const std::string & cmd){
    using namespace std::placeholders;
    cmd_ = cmd;

    tout << "Quantifying UMIs\n";

    QuantWorker::cb_read read_cb = std::bind(&QuantBase::read_, this, _1, _2);
    std::list<QuantWorker> threads;
    for(size_t i = 0; i < num_threads; i++){
        threads.emplace_back(index_, group_map_, group_names_.size(), read_cb, umi_len, in_.total_barcodes(), bam);
    }

    size_t NG = group_names_.size();

    for(auto it = std::next(threads.begin()); it != threads.end(); it++){
        it->start();
    }

    threads.front()();

    for(auto it = std::next(threads.begin()); it != threads.end(); it++){
        it->join();
    }

    in_.debug();

    tout << "Summarizing\n";
    molecule_counts.resize(in_.total_barcodes() * QuantWorker::dedup::ELEM_COUNT);
    group_counts.resize(in_.total_barcodes() * NG);

    std::vector<size_t> bcounts(in_.total_barcodes());
    size_t N = 0;
    size_t R = 0, R2 = 0, R3 = 0, R4 = 0;
    for(auto it = threads.begin(); it != threads.end(); it++){
        R += it->total_reads;
        R2 += it->dups_.total_reads;
        R3 += it->dups_.total_reads2;
        R4 += it->dups_.total_reads3;
        for(size_t i = 0; i < molecule_counts.size(); i++) molecule_counts[i] += it->dups_.molecule_counts[i];
        for(size_t i = 0; i < group_counts.size(); i++) group_counts[i] += it->dups_.group_counts[i];
        for(auto & c : it->dups_.gene_counts) {
            bcounts[c.barcode]++;
            N++;
        }
        if(bam) umi_correct_.insert(umi_correct_.end(), it->dups_.umi_correct.begin(), it->dups_.umi_correct.end());
    }
    std::cout << "R = " << R << " R2 = " << R2 << " R3 = " << R3 << " R4 = " << R4 << "\n";

    //tout << "Index sorting the molecule counts N = " << N << " MC = " << molecule_counts.size() << "\n";
    //Indexes of where each barcode should go to sort them by their total molecules
    bpos_ = sort_molecule_indexes(molecule_counts, QuantWorker::dedup::ELEM_COUNT);

    // This tells me the order of the barcode items ie. the one with the highest counts is bpos[0]
    std::vector<size_t> bindex(in_.total_barcodes() + 1);
    std::vector<size_t> brev(in_.total_barcodes());
    for(size_t i = 1; i < bindex.size(); i++){
        bindex[i] = bcounts[bpos_[i - 1]] + bindex[i - 1];
        brev[bpos_[i - 1]] = i - 1;
    }
    //tout << "Index N = " << bindex.back() << "\n";
    // Just write each tag to a gzipped file
    gene_counts.resize(bindex.back());
    //tout << "Merging the tags\n";
    auto it = threads.begin();
    total_ = 0;
    while(it != threads.end()){
        auto bit = it->dups_.gene_counts.begin();
        while(bit != it->dups_.gene_counts.end()){
            //auto lbarcode = bit->barcode;
            auto idx = brev[bit->barcode];
            bool cc = molecule_counts[bit->barcode * QuantWorker::dedup::ELEM_COUNT] >= min_molecules;
            total_ += cc;
            for(size_t i = bindex[idx]; i < bindex[idx + 1]; i++){
                //if(lbarcode != bit->barcode) std::cout << "Error with logic lgid = " << lbarcode << " gid = " << bit->barcode << "\n";
                if(cc) {
                    gene_id_map_.insert(std::make_pair(bit->gene_id, 0));
                    //if(bit->molecules == 1) barcode_low_gene_counts[bit->barcode]++;
                }
                gene_counts[i] = *bit++;
            }
        }
        it = threads.erase(it);
    }

    size_t idx = 0;
    for(auto & m : gene_id_map_) m.second = idx++;

    tout << "Total genes detected " << gene_id_map_.size() << " from " << total_ << " barcodes with at least " << min_molecules << " molecules\n";
}

unsigned int QuantBase::read_(size_t N, std::vector<AlignSummary> & aligns) {
    std::lock_guard<std::mutex> lock(mtx_read_);
    return in_.read_N(N, aligns);
}

void QuantWorker::operator()() {
    unsigned int N = cb_(100, aligns_);
    while(N > 0){
        std::sort(aligns_.begin(), aligns_.end());
        auto start = aligns_.begin();
        for(auto it = std::next(start); it != aligns_.end(); it++){
            // If the next barcode is different or we are on the last barcode
            if(it->barcode != start->barcode){
                dups_.process(start, it);
                total_reads += (it - start);
                start = it;
            }
        }
        if(start != aligns_.end()) {
            dups_.process(start, aligns_.end());
            total_reads += (aligns_.end() - start);
        }

        N = cb_(100, aligns_);
    }
}

void QuantBase::write_umi_map_(const std::string & out_file){
    if(umi_correct_.empty()) return;

    tout << "Total umi mappings " << umi_correct_.size() << "\n";
    // Sort the mappings by their tid, position, and barcode to make deconvoluting the bam file easier
    std::sort(umi_correct_.begin(), umi_correct_.end(), 
        [&](const UMIMap & m1, UMIMap & m2) { 
        return 
            std::tie(m1.gene_id, m1.barcode, m1.umi_from) < 
            std::tie(m2.gene_id, m2.barcode, m2.umi_from);
        }
    );
    tout << "Writing the umi map data\n";

    gzofstream zout(out_file + "_umi_map.txt.gz");
    for(auto & u : umi_correct_){
        zout << in_.barcodes[u.barcode] << "\t" << u.gene_id << "\t" << u.umi_from << "\t" << u.umi_to << "\n";
    }
    zout.close();
}

void QuantBase::write_output(const std::string & out_file, unsigned int min_molecules){
    write_umi_map_(out_file);
    std::vector<uint32_t> barcode_map;
    barcode_map.resize(in_.total_barcodes());
    std::vector<const char *> gene_names;
    std::vector<const char *> gene_ids;
    std::vector<const char *> barcodes;
    std::vector<uint32_t>     barcode_ids;
    size_t idx = 0;
    for(auto i : bpos_){
        if(molecule_counts[i * QuantWorker::dedup::ELEM_COUNT] < min_molecules) continue;
        barcodes.push_back(in_.barcodes[i].c_str());
        barcode_ids.push_back(i);
        barcode_map[i] = idx++;
    }

    for(auto const & m : gene_id_map_) {
        gene_names.push_back(index_.gene(m.first).gene_name.c_str());
        gene_ids.push_back(index_.gene(m.first).gene_id.c_str());
    }

    using namespace H5;
    H5File file(out_file + ".h5", H5F_ACC_TRUNC);
    {
        // Write out all the molecule counts, group counts, and alignment rates for each barcode
        tout << "Writing barcode rates\n";
        H5::Group group(file.createGroup("/barcode_rates"));
        std::array<double, QuantWorker::dedup::ELEM_COUNT> totals{};
        std::vector<uint32_t> t1;
        std::string k;
        std::vector<std::string> order;
        for(size_t i = 0; i < QuantWorker::dedup::ELEM_COUNT; i++){
            totals[i] = 0;
            t1.clear();
            for(auto b : barcode_ids) { 
                t1.push_back(molecule_counts[b * QuantWorker::dedup::ELEM_COUNT + i]);
                totals[i] += t1.back();
            }
            k = QuantWorker::dedup::ctype2str(static_cast<QuantWorker::dedup::CountType>(i));
            order.push_back(k);
            write_h5_numeric(k, t1, group, PredType::NATIVE_UINT32);
        }
        std::cout 
            << "  cDNA UMI Duplicate Rate =     " << (100.0 - 100.0 * totals[QuantWorker::dedup::MOLECULES] / totals[QuantWorker::dedup::READS]) << "\n"
            << "  cDNA PCR Duplicate Rate =     " << (100.0 * totals[QuantWorker::dedup::PCR_DUPS] / totals[QuantWorker::dedup::READS]) << "%\n"
            << "  Intronic UMI Duplicate Rate = " << (100.0 - 100.0 * totals[QuantWorker::dedup::INTRONIC] / totals[QuantWorker::dedup::IREADS]) << "%\n"
            << "  Intronic PCR Duplicate Rate = " << (100.0 * totals[QuantWorker::dedup::IPCR_DUPS] / totals[QuantWorker::dedup::IREADS]) << "%\n"
            << "  Total UMI Duplicate Rate =    " << (100.0 - 100.0 * (totals[QuantWorker::dedup::MOLECULES] + totals[QuantWorker::dedup::INTRONIC]) / (totals[QuantWorker::dedup::READS] + totals[QuantWorker::dedup::IREADS])) << "%\n"
            << "  Total PCR Duplicate Rate =    " << (100.0 * (totals[QuantWorker::dedup::PCR_DUPS] + totals[QuantWorker::dedup::IPCR_DUPS]) / (totals[QuantWorker::dedup::READS] + totals[QuantWorker::dedup::IREADS])) << "%\n";

        std::cout << "  Total Reads Used = " << std::fixed << (totals[QuantWorker::dedup::READS] + totals[QuantWorker::dedup::IREADS]) << "\n";
        std::cout << "  Total Molecules  = " << std::fixed << (totals[QuantWorker::dedup::MOLECULES] + totals[QuantWorker::dedup::INTRONIC]) << "\n";
        std::cout << "  Total PCR Dups   = " << std::fixed << (totals[QuantWorker::dedup::IPCR_DUPS] + totals[QuantWorker::dedup::PCR_DUPS]) << "\n";

        auto NG = group_names_.size();
        for(size_t i = 0; i < NG; i++){
            t1.clear();
            for(auto b : barcode_ids) t1.push_back(group_counts[b * NG + i]);
            k = "group_" + group_names_[i];
            write_h5_numeric(k, t1, group, PredType::NATIVE_UINT32);
            order.push_back(k);
        }


        t1.clear();
        for(auto i : barcode_ids){
            size_t t = 0;
            for(size_t j = 1; j < AlignGroup::ELEM_COUNT; j++) t += in_.arates[i * AlignGroup::ELEM_COUNT + j];
            t1.push_back(t);
        }
        order.push_back("align_total_reads");
        write_h5_numeric("align_total_reads", t1, group, PredType::NATIVE_UINT32);

        for(size_t i = 0; i < AlignGroup::ELEM_COUNT; i++){
            t1.clear();
            for(auto b : barcode_ids) t1.push_back(in_.arates[b * AlignGroup::ELEM_COUNT + i]);
            if(static_cast<AlignGroup::Result>(i) == AlignGroup::BARCODE_FAIL){
                order.push_back("align_barcode_corrected");
            }else{
                std::string o = "align_";
                order.push_back(o + AlignGroup::alignres2str(static_cast<AlignGroup::Result>(i), false));
            }
            write_h5_numeric(order.back(), t1, group, PredType::NATIVE_UINT32);
        }

        std::vector<const char *> corder;
        for(auto & o : order) corder.push_back(o.c_str());
        write_h5_string("field_order", corder, group);
        corder.clear();
    }


    {
        // Write the file specific rates
        tout << "Writing file rates\n";
        in_.combine_totals();
        H5::Group group(file.createGroup("/file_rates"));
        std::vector<uint32_t> t1;
        std::vector<const char *> fnames;
        for(auto & f : in_.fnames) fnames.push_back(f.c_str());
        write_h5_string("file_prefixes", fnames, group);
        std::vector<std::string> order;
        for(size_t j = 0; j < AlignGroup::ELEM_COUNT; j++){
            t1.clear();
            for(size_t i = 0; i < in_.fnames.size(); i++){
                t1.push_back(in_.arates[i * AlignGroup::ELEM_COUNT + j]);
            }
            order.push_back(AlignGroup::alignres2str(static_cast<AlignGroup::Result>(j), false));
            write_h5_numeric(order.back(), t1, group, PredType::NATIVE_UINT32);
        }
        std::vector<const char *> corder;
        for(auto & o : order) corder.push_back(o.c_str());
        write_h5_string("field_order", corder, group);
    }

    write_h5_numeric("barcode_ids", barcode_ids, file, PredType::NATIVE_UINT32);
    write_h5_string("gene_names", gene_names, file);
    write_h5_string("gene_ids", gene_ids, file);
    write_h5_string("barcodes", barcodes, file);

    //Build the molecule count matrix
    std::vector<uint32_t> rows, cols, data;
    {
        for(auto & g : gene_counts){
            if(molecule_counts[g.barcode * QuantWorker::dedup::ELEM_COUNT] < min_molecules || g.molecules == 0) continue;
            rows.push_back(gene_id_map_[g.gene_id]);
            cols.push_back(barcode_map[g.barcode]);
            data.push_back(g.molecules);
        }
        H5::Group group(file.createGroup("/exonic"));

        tout << "Writing cDNA counts\n";
        write_h5_numeric("data", data, group, PredType::NATIVE_UINT32);
        write_h5_numeric("rows", rows, group, PredType::NATIVE_UINT32);
        write_h5_numeric("cols", cols, group, PredType::NATIVE_UINT32);
    }
    {
        rows.clear(); cols.clear(); data.clear();
        for(auto & g : gene_counts){
            if(molecule_counts[g.barcode * QuantWorker::dedup::ELEM_COUNT] < min_molecules || g.intronic == 0) continue;
            rows.push_back(gene_id_map_[g.gene_id]);
            cols.push_back(barcode_map[g.barcode]);
            data.push_back(g.intronic);
        }
        if(!rows.empty()){
            tout << "Writing intronic counts\n";
            H5::Group group(file.createGroup("/intronic"));
            write_h5_numeric("data", data, group, PredType::NATIVE_UINT32);
            write_h5_numeric("rows", rows, group, PredType::NATIVE_UINT32);
            write_h5_numeric("cols", cols, group, PredType::NATIVE_UINT32);
        }
    }

    file.close();
    tout << "Done\n";
}

