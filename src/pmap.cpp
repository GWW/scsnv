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
#include "pmap.hpp"
#include <memory>
#include <unordered_map>
#include "map_base.hpp"
#include "quant_worker.hpp"
#include "sbam_merge.hpp"
#include "bam_genes.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include "misc.hpp"
#include <unistd.h>
#include "htslib/htslib/hts.h"
#include "htslib/htslib/thread_pool.h"
#include "htslib/htslib/sam.h"

using namespace gwsc;

argagg::parser ProgMap::parser() const {
    argagg::parser argparser {{
        { "help", {"-h", "--help"},
          "shows this help message", 0},
        { "txidx", {"-i", "--index"},
          "Transcript Index", 1},
        { "genome", {"-g", "--genome"},
          "Genome BWA mem index", 1},
        { "barcodes", {"-b", "--barcodes"},
          "Barcode count prefix", 1},
        { "threads", {"-t", "--threads"},
          "Number of threads", 1},
        { "qthreads", {"-q", "--quant-threads"},
          "Number of threads to use for quantification", 1},
        { "output", {"-o", "--output"},
          "Output prefix", 1},
        { "library", {"-l", "--library"},
          "libary type (V2)", 1},
        //{ "internal", {"--internal"},
        //  "Look for internal polyA sites with a length of at least 25 bp with 1 mismatch. Useful when the read lengths are longer than the normal 98 bp", 1},
        { "dust", {"-d", "--dust"},
          "Dust complexity cutoff (-1 to disable, default 4)", 1},
        { "overhang", {"--overhang"},
          "Trim reads with X bp or less of overlap with an exon-intron junction or reads with a terminal splice site of X bp or less", 1},
        { "no_bam", {"--no-bam"},
          "Disable writing the sorted bam files of the countable (ie. uniquely mapped reads)", 0},
        //{ "wtags", {"--tags"},
        //  "Write alignment and UMI correction tags (for debugging purposes)", 0},
        { "bam_tmp", {"--bam-tmp"},
          "Temporary directory to store sorted bam files (Default: {out_prefix}_btmp", 1},
        { "cgroups", {"-c", "--count-groups"},
          "Gene Groups for cell quantification", 1},
        { "bam_per_thread", {"--bam-thread"},
          "Number of output reads to buffer for each thread (Default: 50000)", 1},
        { "bam_per_file", {"--bam-file"},
          "Number of output reads per file (Default: 5000000)", 1},
        { "bam_write", {"--bam-write"},
          "Number of writer threads to use when emitting sorted bam files (Default 1)", 1},
      }};
    return argparser;
}

std::string ProgMap::usage() const {
    return "scsnv map -i <transcript index prefix> -g <genome bwa index> -b <barcode prefix> -o <out prefix> <fastq folder 1> ... <fastq folder N>";
}

void ProgMap::load() {
    tx_idx_ = args_["txidx"].as<std::string>();
    bc_counts_ = args_["barcodes"].as<std::string>();
    out_prefix_ = args_["output"].as<std::string>();
    lib_type_ = args_["library"].as<std::string>("V2");
    tmp_bam_ = args_["bam_tmp"].as<std::string>("");
    gene_groups_ = args_["cgroups"].as<std::string>("");
    threads_ = args_["threads"].as<unsigned int>(1);
    qthreads_ = args_["qthreads"].as<unsigned int>(1);
    min_overhang_ = args_["overhang"].as<unsigned int>(5);
    dust_ = args_["dust"].as<double>(4.0);
    bam_write_threads_ = args_["bam_write"].as<unsigned int>(1);
    bam_ = !args_["no_bam"];
    //internal_ = args_["internal"];
    //write_tags_ = args_["wtags"];
    if(bam_){
        if(tmp_bam_.empty()){
            tmp_bam_ = out_prefix_ + "_btmp";
        }
        struct stat st = {};
        if (stat(tmp_bam_.c_str(), &st) == -1) {
            mkdir(tmp_bam_.c_str(), 0700);
        }

        auto bams = glob(tmp_bam_ + "/scsnv_tmp_*.bam");
        if(!bams.empty()){
            tout << "Removing " << bams.size() << " temporary bam files from " << tmp_bam_ << "\n";
            for(auto & f : bams){
                unlink(f.c_str());
            }
        }
        bam_per_thread_ = args_["bam_per_thread"].as<unsigned int>(50000);
        if(bam_per_thread_ < 500){
            std::cout << "Minimum bam-thread is 500\n";
            bam_per_thread_ = 500;
        }

        bam_per_file_ = args_["bam_per_file"].as<unsigned int>(5000000);
        //if(bam_per_file_ < 100000){
        //    std::cout << "Less than 100000 reads per file is not recommended setting to 100000\n";
        //    bam_per_file_ = 100000;
        //}

        if(bam_per_file_ % bam_per_thread_ != 0){
            std::cout << "--bam-file must be a multiple of --bam-thread";
            exit(1);
        }
    }

    genome_idx_ = args_["genome"].as<std::string>();
    //if(args_["genome"]){
    //}

    if(args_.pos.size() == 0){
        throw std::runtime_error("Missing fastq folder argument(s)");
    }
    for(auto & f : args_.pos){
        std::string d = f;
        while(d.rbegin() != d.rend() && *d.rbegin() == '\\') d.pop_back();
        dirs_.push_back(d);
    }
    std::sort(dirs_.begin(), dirs_.end(), ReadStringCmp());
    fastqs_ = find_fastq_files(dirs_);
}

struct BamData{
    BamData(){
        b = bam_init1();
    }

    ~BamData(){
        if(b != nullptr) bam_destroy1(b);
        b = nullptr;
    }
    bam1_t  * b = nullptr;

    bool operator<(const BamData & d) const {
        return std::tie(gid, barcode, umi, pos) < 
            std::tie(d.gid, d.barcode, d.umi, d.pos);
    }

    bool operator==(const BamData & d) const {
        return std::tie(gid, barcode, umi, pos) ==
            std::tie(d.gid, d.barcode, d.umi, d.pos);
    }

    uint64_t           pos = 0;
    int                tid = 0;
    AlignSummary::bint barcode = 0;
    uint32_t           umi = 0;
    uint32_t           gid = 0;
    char               ra = 0;
};

size_t process_dups(std::vector<BamData *>::const_iterator start, std::vector<BamData *>::const_iterator end, 
        std::vector<BamData*> & tmp){
    tmp.clear();
    size_t dups = 0;
    for(auto it = start; it != end; it++){
        if((*it)->ra == 'E' || (*it)->ra == 'N'){
            tmp.push_back(*it);
        }
    }

    if(tmp.size() < 2) return dups;

    std::sort(tmp.begin(), tmp.end(), 
            [](const BamData * p1, const BamData * p2) {
                return *p1 < *p2; 
            }
    );

    auto check = tmp.begin();
    for(auto it = std::next(tmp.begin()); it != tmp.end(); it++){
        if((**it) == (**check)){
            (*check)->b->core.flag |= BAM_FDUP;
            dups++;
        }else{
            check = it;
        }
    }

    return dups;
}

void ProgMap::write_progress_(){
    if((written_ - lwritten_) < 5000000) 
        return;

    lwritten_ = written_;

    size_t sec = tout.seconds();
    double ps = 1.0 * written_ / (sec - start_);
    size_t eta = (total_ - written_) / ps;
    int hours = eta / (60 * 60);
    int minutes = (eta - (hours * 60 * 60)) / 60;
    tout << "Merged " << written_ << " / " << total_ << " [" << static_cast<int>(ps) << " / sec], ETA = " 
        << hours << "h " << minutes << "m\n";

}

template <typename T>
int ProgMap::run_wrap_(){
    std::vector<std::string> bstrings;
    MapBase<T> base;
    base.load_barcode_counts(bc_counts_, fastqs_);
    base.load_index(tx_idx_, min_overhang_, genome_idx_);
    if(bam_){
        base.prepare_bam(full_cmd_, bam_per_thread_, bam_per_file_, tmp_bam_, bam_write_threads_);
    }
    base.run(threads_, fastqs_, dust_, internal_);
    base.write_output(out_prefix_);
    base.unload();
    total_ = base.total_reads();
    std::cout << "\n";
    tout << "Beginning quantification step\n";
    gwsc::QuantBase qbase;
    qbase.set_index(base.index());
    base.get_barcodes(bstrings);
    qbase.set_data(base.aligns, bstrings.size());
    qbase.build_gene_groups(gene_groups_);
    qbase.run(T::UMI_LEN, qthreads_, bam_, full_cmd_);
    qbase.write_output(out_prefix_ + "summary.h5", bstrings, base.brates);
    tout << "Quantification done\n\n";

    if(!bam_) {
        tout << "Done\n";
        return EXIT_SUCCESS;
    }

    auto & umi_correct = qbase.umi_correct;
    auto & umi_bad = qbase.umi_bad;
    std::sort(umi_correct.begin(), umi_correct.end());
    std::sort(umi_bad.begin(), umi_bad.end());
    std::vector<uint32_t> bindex(bstrings.size() + 1);
    std::vector<uint32_t> bad_index(bstrings.size() + 1);
    for(auto & u : umi_correct) bindex[u.barcode + 1]++;
    for(auto & u : umi_bad) bad_index[u.barcode + 1]++;
    for(size_t i = 1; i < bindex.size(); i++) bindex[i] += bindex[i - 1];
    for(size_t i = 1; i < bad_index.size(); i++) bad_index[i] += bad_index[i - 1];

    tout << "Merging bam files and correcting UMIs\n";
    tout << "Total UMI corrections " << umi_correct.size() << "\n";
    tout << "Total bad UMI combinations " << umi_bad.size() << "\n";

    /*
    if(write_tags_){
        base.write_tags(out_prefix_ + "align", umi_correct, bindex);
        qbase.write_umi_map(out_prefix_ + "tag_corrections.gz");
    }
    */


    start_ = tout.seconds();

    spp::sparse_hash_map<std::string, size_t> bhash;
    std::vector<bool> ccheck(umi_correct.size());
    for(size_t i = 0; i < bstrings.size(); i++) bhash[bstrings[i]] = i;
    std::string  d = tmp_bam_ + "/scsnv_tmp_*.bam";
    auto bam_files = glob(d);
    tout << "Found " << bam_files.size() << " bam files to merge\n";
    BamMerger bm;
    bm.add_bams(bam_files.begin(), bam_files.end());

    tout << "Merging, writing and correcting alignments\n";
    std::string outf = out_prefix_ + "merged.bam";
    samFile * bam_out = sam_open(outf.c_str(), "wb");

    if(bam_write_threads_ > 1){
        hts_set_threads(bam_out, bam_write_threads_);
    }

    std::cout << bm.header() << " " << bm.header()->n_targets << "\n";
    if(sam_hdr_write(bam_out, bm.header()) < 0) {
        std::cerr << "Error writing header\n";
        exit(1);
    }

    BamData * next = new BamData;
    std::vector<BamData*> reads, treads;
    for(size_t i = 0; i < 20; i++) reads.push_back(new BamData);
    size_t bidx = 0;
    int ltid = -1;
    int lpos = 0;
    std::string key;
    bool has_unmapped = false;
    size_t total_dups = 0;
    size_t total_counted = 0;
    size_t discarded = 0;
    while(bm.next(next->b) != nullptr){
        BamData & bd = *next;
        if(bd.b->core.flag & BAM_FUNMAP){
            has_unmapped = true;
            break;
        }else{
            bd.tid = bd.b->core.tid;
            bd.gid = std::numeric_limits<uint32_t>::max();
            bd.ra = bam_aux2A(bam_aux_get(bd.b, "RA"));
            bd.barcode = bhash[bam_aux2Z(bam_aux_get(bd.b, "CB"))];
            bd.pos = (static_cast<uint64_t>(bd.b->core.pos) << 32) | (bam_endpos(bd.b) - 1);
            key = bam_aux2Z(bam_aux_get(bd.b, "UB"));
            seq2int<gwsc::ADNA4, uint32_t>(key, bd.umi);

            if(bd.ra == 'E' || bd.ra == 'N'){
                total_counted++;
                bd.gid = base.index().gene_from_id(bam_aux2Z(bam_aux_get(bd.b, "XG"))).gid;
                UMIMap m;
                m.gene_id = bd.gid;
                m.barcode = bd.barcode;
                m.umi_from = bd.umi;
                auto bs = umi_correct.begin() + bindex[bd.barcode];
                auto be = umi_correct.begin() + bindex[bd.barcode + 1];
                if(bs != be){
                    auto rit = std::lower_bound(bs, be, m);
                    if(rit != be && (*rit) == m){
                        size_t ridx = std::distance(umi_correct.begin(), rit);
                        ccheck[ridx] = true;
                        uint32_t mask = (1 << ADNA4::size_) - 1;
                        uint8_t * cumi = bam_aux_get(bd.b, "UB") + 1;
                        uint32_t val = rit->umi_to;
                        bd.umi = val;
                        for(size_t i = 0; i < T::UMI_LEN; i++){
                            cumi[T::UMI_LEN - i - 1] = ADNA4::alphabet_str_[val & mask];
                            val >>= ADNA4::size_;
                        }
                        bam_aux_append(bd.b, "UR", 'Z', key.length() + 1, reinterpret_cast<uint8_t*>(const_cast<char*>(key.c_str())));
                    }
                }
                auto bad_s = umi_bad.begin() + bad_index[bd.barcode];
                auto bad_e = umi_bad.begin() + bad_index[bd.barcode + 1];
                UMIBad bm(bd.barcode, bd.gid, bd.umi);
                auto brit = std::lower_bound(bad_s, bad_e, bm);
                if(brit != bad_e && (*brit) == bm){
                    bd.ra = '?';
                    uint8_t * rptr = bam_aux_get(bd.b, "RA") + 1;
                    (*rptr) = '?';
                    discarded++;
                }
            }

            if(bd.tid != ltid || bd.b->core.pos != lpos){
                if(bidx > 0){
                    total_dups += process_dups(reads.begin(), reads.begin() + bidx, treads);
                    for(size_t i = 0; i < bidx; i++){
                        if(sam_write1(bam_out, bm.header(), reads[i]->b) < 0){
                            std::cerr << "Error writing sam\n"; exit(1);
                        }
                        written_++;
                        write_progress_();
                    }
                }
                ltid = bd.tid;
                lpos = bd.b->core.pos;
                bidx = 0;
            }
            if(bidx >= reads.size()){
                for(size_t i = 0; i < 20; i++) reads.push_back(new BamData);
            }
            std::swap(reads[bidx], next);
            bidx++;

            /*
            if(sam_write1(bam_out, bm.header(), bd.b) < 0){
                std::cerr << "Error writing sam\n";
                exit(1);
            }
            */
        }
    }


    if(bidx > 0){
        total_dups += process_dups(reads.begin(), reads.begin() + bidx, treads);
        for(size_t i = 0; i < bidx; i++){
            if(sam_write1(bam_out, bm.header(), reads[i]->b) < 0){
                std::cerr << "Error writing sam\n"; exit(1);
            }
            written_++;
            write_progress_();
        }
    }
    tout << "Done writing mapped reads\n";
    if(has_unmapped){
        tout << "Writing unmapped reads\n";
        written_++;
        if(sam_write1(bam_out, bm.header(), next->b) < 0){
            std::cerr << "Error writing sam\n"; exit(1);
        }
        while(bm.next(next->b) != nullptr){
            if(sam_write1(bam_out, bm.header(), next->b) < 0){
                std::cerr << "Error writing sam\n"; exit(1);
            }
            written_++;
            write_progress_();
        }
        tout << "Done writing unmapped reads\n";
    }

    size_t missing = 0;
    for(auto c : ccheck){
        if(!c){
            missing++;
        }
    }
    if(missing > 0){
        std::cout << "Missed " << missing << " out of " << umi_correct.size() << " UMI corrections\n";
    }
    //std::cout << "Duplicate check " << (100.0 * total_dups / total_counted) << "\n";
    tout << "Done writing " << written_ << " reads and " << discarded << " marked as discarded\n";

    //std::cout << "Read destroyed\n";
    //std::cout << "Bam Closed\n";
    sam_close(bam_out);

    tout << "Deleting the temporary bam files\n";
    for(auto & b : bam_files){
        unlink(b.c_str());
    }

    tout << "Done\n";

    return EXIT_SUCCESS;
}

int ProgMap::run() {
    if(lib_type_ == "V2"){
        return run_wrap_<Reader10X_V2>();
    }else if(lib_type_ == "V3"){
        return run_wrap_<Reader10X_V3>();
    }else if(lib_type_ == "V2_5P"){
        return run_wrap_<Reader10X_V2_5P>();
    }else if(lib_type_ == "V3_5P"){
        return run_wrap_<Reader10X_V3_5P>();
    }
    return EXIT_SUCCESS;
}

