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

#include "reader.hpp"
#include <iostream>
#include <mutex>

using namespace gwsc;

gwsc::FastqPairs gwsc::find_fastq_files(std::vector<std::string> & dirs) {
    FastqPairs files;
    for(auto & d : dirs){
        auto r1 = glob(d + "/*_R1_[0-9][0-9][0-9].fastq.gz");
        auto r2 = glob(d + "/*_R2_[0-9][0-9][0-9].fastq.gz");
        std::sort(r1.begin(), r1.end(), ReadStringCmp());
        std::sort(r2.begin(), r2.end(), ReadStringCmp());
        for(size_t i = 0; i < r1.size(); i++){
            files.push_back({r1[i], r2[i], d});
        }
    }
    std::sort(files.begin(), files.end());
    return files;
}

gwsc::FastqPairs gwsc::find_fastq_files(const std::string & dir) {
    FastqPairs files;
    auto r1 = glob(dir + "/*_R1_[0-9][0-9][0-9].fastq.gz");
    auto r2 = glob(dir + "/*_R2_[0-9][0-9][0-9].fastq.gz");
    auto i1 = glob(dir + "/*_I1_[0-9][0-9][0-9].fastq.gz");
    auto i2 = glob(dir + "/*_I2_[0-9][0-9][0-9].fastq.gz");
    std::sort(r1.begin(), r1.end(), ReadStringCmp());
    std::sort(r2.begin(), r2.end(), ReadStringCmp());
    for(size_t i = 0; i < r1.size(); i++){
        files.push_back({r1[i], r2[i], dir});
        if(!i1.empty()) files.back().i1 = i1[i];
        if(!i2.empty()) files.back().i2 = i2[i];
    }
    std::sort(files.begin(), files.end());
    return files;
}


unsigned int ReaderBase::read_N(size_t N, Reads & reads) {
    if(reads.size() < N) reads.resize(N);    
    size_t i = 0;
    while(i < N && read(reads[i])){
        i++;
    }
    return i;
}


bool Reader10X_V2::read(Read & r) {
    bool res1 = f1_.read(r1_);
    bool res2 = f2_.read(r.name, r.comment, r.tag, r.q_tag);
    if(res1 != res2){
        std::cout << "Error! Inconsistent number of reads between the read1 and read 2 files\n";
        exit(0);
    }else if(!res1){
        return false;
    }
    r.barcode.assign(r1_.seq.begin(), r1_.seq.begin() + 16);
    r.umi.assign(r1_.seq.begin() + 16, r1_.seq.begin() + 26);
    r.obarcode.clear();
    r.q_barcode.assign(r1_.qual.begin(), r1_.qual.begin() + 16);
    r.q_umi.assign(r1_.qual.begin() + 16, r1_.qual.begin() + 26);

    return true; 
}

bool Reader10X_V2::read_barcode(std::string & bc, std::string & qbc) {
    bool res1 = f1_.read(r1_);
    if(!res1) return false;
    bc.assign(r1_.seq.begin(), r1_.seq.begin() + 16);
    qbc.assign(r1_.qual.begin(), r1_.qual.begin() + 16);
    return true; 
}


bool Reader10X_V3::read(Read & r) {
    bool res1 = f1_.read(r1_);
    bool res2 = f2_.read(r.name, r.comment, r.tag, r.q_tag);
    if(res1 != res2){
        std::cout << "Error! Inconsistent number of reads between the read1 and read 2 files\n";
        exit(0);
    }else if(!res1){
        return false;
    }
    r.barcode.assign(r1_.seq.begin(), r1_.seq.begin() + 16);
    r.obarcode.clear();
    r.umi.assign(r1_.seq.begin() + 16, r1_.seq.begin() + 28);
    r.q_barcode.assign(r1_.qual.begin(), r1_.qual.begin() + 16);
    r.q_umi.assign(r1_.qual.begin() + 16, r1_.qual.begin() + 28);

    return true; 
}

bool Reader10X_V3::read_barcode(std::string & bc, std::string & qbc) {
    bool res1 = f1_.read(r1_);
    if(!res1) return false;
    bc.assign(r1_.seq.begin(), r1_.seq.begin() + 16);
    qbc.assign(r1_.qual.begin(), r1_.qual.begin() + 16);
    return true; 
}
