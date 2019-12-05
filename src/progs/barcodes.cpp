#include "barcodes.hpp"
#include "../util/gzstream.hpp"
using namespace gwsc;

argagg::parser ProgBarcodes::parser() const {
    argagg::parser argparser {{
        { "known", {"-k", "--known-barcodes"},
          "Known Barcode File", 1},
        { "out", {"-o", "--out"},
          "Barcode count output file", 1},
        { "library", {"-l", "--library"},
          "libary type (V2, V3)", 1},
        { "help", {"-h", "--help"},
          "shows this help message", 0},
      }};
    return argparser;
}

std::string ProgBarcodes::usage() const {
    return "scmap count -k known_barcodes.txt -o barcodes.gz <fastq folder 1> <fastq folder 2> ...";
}

void ProgBarcodes::load() {
    known_ = args_["known"].as<std::string>("");
    out_ = args_["out"].as<std::string>();
    lib_type_ = args_["library"].as<std::string>("V2");
    if(args_.pos.size() == 0){
        throw std::runtime_error("Missing fastq folder argument(s)");
    }
    for(auto & f : args_.pos){
        std::string d = f;
        while(d.rbegin() != d.rend() && *d.rbegin() == '\\') d.pop_back();
        dirs_.push_back(d);
    }
    std::sort(dirs_.begin(), dirs_.end(), ReadStringCmp());
}

int ProgBarcodes::run() {
    if(lib_type_ == "V2"){
        return run_10X_<Reader10X_V2>();
    }else if(lib_type_ == "V3"){
        return run_10X_<Reader10X_V3>();
    }
    return EXIT_FAILURE;
}
