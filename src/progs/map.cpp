#include "map.hpp"
#include "../scmap/map_base.hpp"
#include <sys/types.h>
#include <sys/stat.h>
#include "../util/misc.hpp"
#include <unistd.h>

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
        { "output", {"-o", "--output"},
          "Output prefix", 1},
        { "library", {"-l", "--library"},
          "libary type (V2)", 1},
        { "dust", {"-d", "--dust"},
          "Dust complexity cutoff (-1 to disable, default 4)", 1},
        { "overhang", {"--overhang"},
          "Trim reads with X bp or less of overlap with an exon-intron junction or reads with a terminal splice site of X bp or less", 1},
        { "no_bam", {"--no-bam"},
          "Disable writing the sorted bam files of the countable (ie. uniquely mapped reads)", 0},
        { "bam_tmp", {"--bam-tmp"},
          "Temporary directory to store sorted bam files (Default: {out_prefix}_tmp", 1},
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
    return "scmap map -i <transcript index prefix> -g <genome bwa index> -b <barcode prefix> -o <out prefix> <fastq folder 1> ... <fastq folder N>";
}

void ProgMap::load() {
    tx_idx_ = args_["txidx"].as<std::string>();
    bc_counts_ = args_["barcodes"].as<std::string>();
    out_prefix_ = args_["output"].as<std::string>();
    lib_type_ = args_["library"].as<std::string>("V2");
    tmp_bam_ = args_["bam_tmp"].as<std::string>("");
    threads_ = args_["threads"].as<unsigned int>(1);
    min_overhang_ = args_["overhang"].as<unsigned int>(5);
    dust_ = args_["dust"].as<double>(4.0);
    bam_write_threads_ = args_["bam_write"].as<unsigned int>(1);
    bam_ = !args_["no_bam"];
    if(bam_){
        if(tmp_bam_.empty()){
            tmp_bam_ = out_prefix_ + "_tmp";
        }
        struct stat st = {};
        if (stat(tmp_bam_.c_str(), &st) == -1) {
            mkdir(tmp_bam_.c_str(), 0700);
        }

        auto bams = glob(tmp_bam_ + "/scmap_tmp_*.bam");
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

template <typename T>
int ProgMap::run_wrap_(){
    MapBase<T> base;
    base.load_barcode_counts(bc_counts_, fastqs_);
    base.load_index(tx_idx_, min_overhang_, genome_idx_);
    if(bam_){
        base.prepare_bam(full_cmd_, bam_per_thread_, bam_per_file_, tmp_bam_, bam_write_threads_);
    }
    base.run(threads_, fastqs_, dust_);
    base.write_output(out_prefix_);
    return EXIT_SUCCESS;
}

int ProgMap::run() {
    if(lib_type_ == "V2"){
        return run_wrap_<Reader10X_V2>();
    }else if(lib_type_ == "V3"){
        return run_wrap_<Reader10X_V3>();
    }
    return EXIT_SUCCESS;
}

