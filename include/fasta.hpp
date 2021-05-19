/*
Copyright (c) 2018-2021 Gavin W. Wilson

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


#pragma once

#include <string>
#include <cctype>
#include <algorithm>
#include <iostream>
#include "read_buffer.hpp"
#include "misc.hpp"
#include "sequence.hpp"

namespace gwsc {

struct Fasta{
    std::string name;
    std::string comment;
    Sequence seq;
    int      tid = -1;
    void clear() {
        name.clear();
        seq.clear();
        comment.clear();
        tid = -1;
    }
};

typedef std::vector<Fasta> Fastas;

class FastaReader {
    public:

        FastaReader(const std::string & fin) {
            open(fin);
        }

        FastaReader() {

        }

        bool open(const std::string & fin){
            std::string gzp = ".gz";
            if(fin.size() > 3 && std::equal(gzp.rbegin(), gzp.rend(), fin.rbegin())){
                buffer = new GzipBuffer;
            }else{
                buffer = new FileBuffer;
            }
            bool flag = buffer->bopen(fin);
            if(!flag){
                std::cerr << "Could not open " << fin << " for reading\n";
                exit(1);
            }
            return flag;
        }

        void close(){
            if(buffer != nullptr) delete buffer;
            buffer = nullptr;
        }

        ~FastaReader(){
            if(buffer != std::nullptr_t()) delete buffer;
            buffer = nullptr;
        }

        bool read(std::string & name, std::string & comment, Sequence & seq){
            name.clear(); comment.clear(); seq.clear();
            if(buffer->last_char == 0){
                int c;
                while((c = buffer->get_char()) != -1 && c != '>'){
                    if(c == -1) return false;
                    buffer->last_char = c;
                }
            }
            if(buffer->get_until(0, name) < 0){
                return false;
            }

            char prev = *std::next(buffer->begin, -1);
            if(prev != '\n' && prev != '\r'){
                buffer->get_until(-1, comment);
            }

            char c = buffer->peek_char();
            while(c != -1 && c != '>'){
                buffer->get_until(-1, seq);
                c = buffer->peek_char();
            }
            return true;
        }

        bool read(Fasta & fa){
            return read(fa.name, fa.comment, fa.seq);
        }

        void read_all(Fastas & fas, bool sort = false){
            fas.push_back(Fasta());
            while(read(fas.back())){
                fas.push_back(Fasta());
            }
            fas.pop_back();
            if(sort){
                auto rs = gwsc::ReadStringCmp();
                std::sort(fas.begin(), fas.end(), [&rs](Fasta & r1, Fasta & r2) { return rs(r1.name, r2.name); });
            }
        }

    private:
        BufferBase            * buffer = nullptr;
};

}
