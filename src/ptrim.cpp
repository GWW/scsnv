
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
#include "ptrim.hpp"
#include "gzstream.hpp"
#include <zlib.h>
using namespace gwsc;

argagg::parser ProgTrim::parser() const {
    argagg::parser argparser {{
        { "tag", {"-t", "--tag-len"},
          "Known Barcode File", 1},
        { "out", {"-o", "--out"},
          "Directory to write trimmed fastqs", 1},
        { "library", {"-l", "--library"},
          "libary type (V2, V3)", 1},
        { "downsample", {"-d", "--downsample"},
          "downsample to X reads (Default 200000000)", 1},
        { "help", {"-h", "--help"},
          "shows this help message", 0},
      }};
    return argparser;
}

std::string ProgTrim::usage() const {
    return "scsnv trim -o <out_folder>/ <fastq folder 1> <fastq folder 2> ...";
}

void ProgTrim::load() {
    tag_ = args_["tag"].as<unsigned int>(98);
    out_ = args_["out"].as<std::string>();
    downsample_ = args_["downsample"].as<unsigned int>(200000000);
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

void trim_fastq(std::string & in, std::string & out, unsigned int len, size_t & total, size_t downsample){
    FileWrapper rin(in);
    std::string line;
    gzFile zout = gzopen(out.c_str(), "wb");
    while(true){
        if(!rin.get_line(line)){
            break;
        }
        line += '\n';
        gzwrite(zout, line.c_str(), line.size());
        if(!rin.get_line(line)){
            break;
        }
        line[len] = '\n';
        gzwrite(zout, line.c_str(), len + 1);
        if(!rin.get_line(line)){
            break;
        }
        line += '\n';
        gzwrite(zout, line.c_str(), line.size());
        if(!rin.get_line(line)){
            break;
        }
        line += '\n';
        line[len] = '\n';
        gzwrite(zout, line.c_str(), len + 1);
        total++;
        if(total % 10000000 == 0) tout << "  Processed " << total << " reads out of " << downsample << "\n";
        if(total >= downsample) break;
    }

    gzclose(zout);
}

template <typename R> 
int ProgTrim::run_10X_() {
    Read read;
    size_t r1_len = R::BARCODE_LEN + R::UMI_LEN;
    size_t r1_total = 0;
    size_t r2_total = 0;
    for(auto & d : dirs_){
        tout << "Processing read directory " << d << "\n";
        auto fastqs = find_fastq_files(d);
        for(auto & f : fastqs){
            std::string bn1 = out_ + f.first.substr(f.first.find_last_of('/'));
            std::string bn2 = out_ + f.second.substr(f.second.find_last_of('/'));
            tout << "Trimming " << f.first << " to " << r1_len << " bp to " << bn1 << "\n";
            trim_fastq(f.first, bn1, r1_len, r1_total, downsample_);
            tout << "Trimming " << f.second << " to " << tag_ << " bp to " << bn2 << "\n";
            trim_fastq(f.second, bn2, tag_, r2_total, downsample_);
            if(r1_total >= downsample_) break;
        }
    }
    return EXIT_SUCCESS;
}

int ProgTrim::run() {
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
