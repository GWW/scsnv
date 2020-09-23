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


#pragma once

#include <string>
#include <cctype>
#include <iostream>
#include <algorithm>
#include "read_buffer.hpp"
#include "sequence.hpp"

namespace gwsc {

struct Fastq{
    std::string name;
    std::string comment;
    Sequence seq;
    std::string qual;
    void clear() {
        name.clear();
        seq.clear();
        comment.clear();
        qual.clear();
    }
};


class FastqReader {
    public:
        FastqReader(const std::string & fin) {
            open(fin);
        }

        FastqReader() {

        }

        bool open(const std::string & fin){
            close();
            static const std::string gzp = ".gz";
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

        ~FastqReader(){
            if(buffer != nullptr) delete buffer;
            buffer = nullptr;
        }

        void close(){
            if(buffer != nullptr) delete buffer;
            buffer = nullptr;
        }

        bool read(Fastq & fa){
            return read(fa.name, fa.comment, fa.seq.str(), fa.qual);
        }

        bool read(std::string & name, std::string & comment, std::string & seq, std::string & qual){
            name.clear();
            comment.clear();
            seq.clear();
            qual.clear();
            if(buffer->last_char == 0){
                int c;
                while((c = buffer->get_char()) != -1 && c != '@'){
                    if(c == -1) return false;
                    buffer->last_char = c;
                }
            }
            if(buffer->get_until(0, name) < 0){
                return false;
            }

            if(*std::next(buffer->begin, -1) != '\n'){
                buffer->get_until('\n', comment);
            }

            buffer->get_until('\n', seq);
            char c = 0;
            while((c = buffer->get_char()) != -1 && c != '\n');
            buffer->get_until('\n', qual);
            buffer->last_char = 0;
            return true;
        }

        BufferBase          * buffer = nullptr;
};

}
