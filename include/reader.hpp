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

#include <string>
#include <mutex>
#include <thread>
#include <random>
#include "fastq.hpp"
#include "misc.hpp"
#include "aux.hpp"
#include "barcodes.hpp"

namespace gwsc {

struct Read {
    std::string name;
    std::string comment;
    std::string barcode;
    std::string obarcode;
    std::string q_barcode;
    std::string umi;
    std::string q_umi;
    std::string tag;
    std::string q_tag;
    unsigned int tend = 0;
};

using Reads = std::vector<Read>;

struct BarcodeRead{
    std::string name;
    std::string comment;
    std::string barcode;
    std::string q_barcode;
};


FastqPairs find_fastq_files(const std::string & dir);
FastqPairs find_fastq_files(std::vector<std::string> & dirs);

template <typename T>
class MultiReader{
    MultiReader( const MultiReader& ) = delete;
    MultiReader& operator=(const MultiReader&) = delete;
    public:
        using ReaderType = T;
        MultiReader(bool verbose = true) : verbose_(verbose) {

        }

        void reset() {
            if(files_.empty()) return;
            curr_ = 0;
            if(verbose_) tout << "Processing " << files_[curr_].first << " and " 
                              << files_[curr_].second << " total reads = " << files_[curr_].total  << "\n";
            in_.open(files_[curr_].first, files_[curr_].second);
            read_.clear();
            read_.resize(files_.size());
        }

        bool read(Read & r);
        unsigned int read_N(size_t N, Reads & reads);
        unsigned int read_N_safe(size_t N, Reads & reads);
        bool read_barcode(std::string & bc, std::string & qbc);

        const FastqPair & file(size_t i) const {
            return files_[i];
        }

        size_t current() const {
            return curr_;
        }

        size_t count(size_t i) const {
            return read_[i];
        }

        size_t files() const {
            return files_.size();
        }

        void set_files(const FastqPairs & files){
            files_ = files;
            curr_ = std::numeric_limits<size_t>::max();
            reset();
        }

        void set_downsample(double perc, size_t seed){
            downsample_ = perc;
            mt_ = std::mt19937(seed);
            dist_ = std::uniform_real_distribution<double>(0, 1);
        }

        size_t skipped() const {
            return skipped_;
        }

    private:
        FastqPairs                 files_;
        T                          in_;
        std::mt19937               mt_;
        std::uniform_real_distribution<double> dist_;
        size_t                     curr_ = std::numeric_limits<size_t>::max();
        size_t                     skipped_ = 0;
        double                     downsample_ = 0.0;
        std::vector<size_t>        read_;
        std::mutex                 mutex_;
        bool                       verbose_;

};

class ReaderBase{
    ReaderBase( const ReaderBase& ) = delete;
    ReaderBase& operator=(const ReaderBase&) = delete;

    public:

        ReaderBase() {

        }

        ReaderBase(const std::string & f1, const std::string & f2) {
            open(f1, f2);
        }

        virtual ~ReaderBase(){

        }

        void open(const std::string & f1, const std::string & f2){
            f1_.open(f1);
            f2_.open(f2);
        }

        virtual bool read(Read & r) = 0;
        virtual bool read_barcode(std::string & bc, std::string & qbc) = 0;
        unsigned int read_N(size_t N, Reads & reads);

    protected:
        FastqReader f1_;
        FastqReader f2_;
};

class Reader10X_V2 : public ReaderBase {
    public:
        const static StrandMode LibraryStrand = TAG_FWD;
        using LibraryBarcode = CBWhiteListShort;
        const static unsigned int UMI_LEN = 10;
        const static unsigned int BARCODE_LEN = 16;

        Reader10X_V2() {

        }

        virtual ~Reader10X_V2(){

        }

        Reader10X_V2(const std::string & f1, const std::string & f2) : ReaderBase(f1, f2){

        }

        virtual bool read(Read & r);
        virtual bool read_barcode(std::string & bc, std::string & qbc);

    protected:
        Fastq       r1_;
};

class Reader10X_V2_5P : public Reader10X_V2 {
    public:
        const static StrandMode LibraryStrand = TAG_REV;

        virtual ~Reader10X_V2_5P() {

        }
};

class Reader10X_V3 : public ReaderBase {
    public:
        const static StrandMode LibraryStrand = TAG_FWD;
        using LibraryBarcode = CBWhiteListShort;
        const static unsigned int UMI_LEN = 12;
        const static unsigned int BARCODE_LEN = 16;

        Reader10X_V3() {

        }

        ~Reader10X_V3(){

        }

        Reader10X_V3(const std::string & f1, const std::string & f2) : ReaderBase(f1, f2){

        }

        virtual bool read(Read & r);
        virtual bool read_barcode(std::string & bc, std::string & qbc);

    protected:
        Fastq       r1_;
};

class Reader10X_V3_5P : public Reader10X_V3 {
    public:
        const static StrandMode LibraryStrand = TAG_REV;

        virtual ~Reader10X_V3_5P() {

        }
};


template <typename T>
inline unsigned int MultiReader<T>::read_N(size_t N, Reads & reads) {
    if(reads.size() < N) reads.resize(N);    
    size_t i = 0;
    while(i < N && read(reads[i])){
        if(downsample_ > 0.0 && dist_(mt_) > downsample_){
            skipped_++;
            continue;
        }
        i++;
    }
    reads.resize(i);
    return i;
}

template <typename T>
inline bool MultiReader<T>::read(Read & read) {
    if(curr_ >= files_.size()){
        return false;
    }
    while(!in_.read(read)){
        curr_++;
        if(curr_ < files_.size()){
            in_.open(files_[curr_].first, files_[curr_].second);
            if(verbose_) tout << "Processing " << files_[curr_].first << " and " << files_[curr_].second << " total reads = " << files_[curr_].total << "\n";
        }else{
            return false;
        }
    }
    read_[curr_]++;
    return true;
}

template <typename T>
inline bool MultiReader<T>::read_barcode(std::string & bc, std::string & qbc) {
    if(curr_ >= files_.size()){
        return false;
    }
    while(!in_.read_barcode(bc, qbc)){
        curr_++;
        if(curr_ < files_.size()){
            in_.open(files_[curr_].first, files_[curr_].second);
        }else{
            return false;
        }
    }
    read_[curr_]++;
    return true;
}

template <typename T>
inline unsigned int MultiReader<T>::read_N_safe(size_t N, Reads & reads) {
    std::lock_guard<std::mutex> lock(mutex_);
    return read_N(N, reads);
}

}
