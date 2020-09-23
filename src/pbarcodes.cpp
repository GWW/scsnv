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
#include "pbarcodes.hpp"
#include "gzstream.hpp"
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
    return "scsnv count -k known_barcodes.txt -o barcodes.gz <fastq folder 1> <fastq folder 2> ...";
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
    }else if(lib_type_ == "V2_5P"){
        return run_10X_<Reader10X_V2_5P>();
    }else if(lib_type_ == "V3_5P"){
        return run_10X_<Reader10X_V3_5P>();
    }
    return EXIT_FAILURE;
}
