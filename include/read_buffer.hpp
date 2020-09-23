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


/* Based on the file parser written by Heng Li, but converted to C++ */

#pragma once

#include <string>
#include <zlib.h>
#include <cstdio>
#include <vector>
#include <functional>
#include <iostream>
#include <algorithm>


namespace gwsc {

typedef std::vector<std::string> ParserTokens;

class BufferBase {
    BufferBase( const BufferBase& ); // non construction-copyable
    BufferBase& operator=( const BufferBase& ); // non copyable

    protected:

    std::string             buffer;
    bool eof  = false;

    public:
    std::string::iterator   begin;
    size_t BUFF_SIZE = 4096;
    int    last_char = 0;

    virtual bool bopen(const std::string & fin) = 0;
    virtual size_t bread() = 0;
    virtual void bclose() = 0;

    BufferBase(){

    }

    virtual ~BufferBase() {

    }

    void init(){
        eof = false;
        last_char = 0;
        buffer.resize(BUFF_SIZE);
        begin = buffer.end();
    }

    bool read_more(){
        std::size_t end = bread();
        buffer.resize(end);
        begin = buffer.begin();
        if(end < BUFF_SIZE){
            eof = true;
        }
        if(end == 0) return false;
        return true;
    }

    int get_char(){
        if(eof && begin == buffer.end()){
            return -1;
        }
        if(begin == buffer.end() && !read_more()){
            return -1;
        }
        return *begin++;
    }

    int peek_char(){
        if(eof && begin == buffer.end()){
            return -1;
        }
        if(begin == buffer.end() && !read_more()){
            return -1;
        }
        return *begin;
    }

    /* 
     c = 0 indicates find until a isspace
    */

    template <typename T>
    int get_until(const char c, T & out){
        if(begin == buffer.end() && eof){
            return -1;
        }
        while(true){
            if(begin == buffer.end()){
                if(eof || !read_more()){
                    break;
                }
            }
            std::string::iterator it;
            if(c != 0){
                it = std::find_if(begin, buffer.end(), [c](char v){ return c == v;});
            }else{
                it = std::find_if(begin, buffer.end(), [](char v){ return isspace(v);});
            }
            out.append(begin, it);
            if(it != buffer.end()) {
                begin = ++it;
                break;
            }else{
                begin = buffer.end();
            }
        }
        return buffer.size();
    }

    template <typename T>
    int get_until(std::function<bool(char)> & fn, T & out){
        if(begin == buffer.end() && eof){
            return - 1;
        }
        while(true){
            if(begin == buffer.end()){
                if(eof || !read_more()){
                    break;
                }
            }
            std::string::iterator it;
            it = std::find_if(begin, buffer.end(), fn);
            out.append(begin, it);
            if(it != buffer.end()) {
                begin = ++it;
                break;
            }else{
                begin = buffer.end();
            }
        }
        return buffer.size();
    }

    int tokenize_line(ParserTokens & toks, char delim='\t'){
        // Not very memory efficient, need to preserve tokens maybe?
        // Use an index and add additional strings as needed ?
        if(begin == buffer.end() && eof){
            return - 1;
        }
        if(toks.empty()){
            toks.push_back(std::string());
        }else{
            toks[0].clear();
        }

        size_t i = 0;
        while(true){
            if(begin == buffer.end()){
                if(eof || !read_more()){
                    break;
                }
            }
            std::string::iterator it;
            it = std::find_if(begin, buffer.end(), [delim](char v){ return v == delim || v == '\n';});
            toks[i].append(begin, it);
            if(it != buffer.end()){
                if(*it == '\n') {
                    begin = ++it;
                    break;
                }
                begin = ++it;
                i++;
                if(i >= toks.size()){
                    toks.push_back(std::string());
                }else{
                    toks[i].clear();
                }
            }else{
                begin = buffer.end();
            }
        }
        return buffer.size();
    }
};

class GzipBuffer : public BufferBase {
    public:

    gzFile      fp;

    bool bopen(const std::string & fin){
        fp = gzopen(fin.c_str(), "r");
        buffer.resize(BUFF_SIZE);
        init();
        if(!fp) {
            eof = true;
            return false;
        }
        return true;
    }

    virtual ~GzipBuffer(){
        bclose();
    }

    void bclose(){
        if(fp != nullptr) gzclose(fp);
        fp = nullptr;
    }

    size_t bread() {
        return gzread(fp, &buffer[0], BUFF_SIZE);
    }

};

class FileBuffer : public BufferBase {
    public:

    FileBuffer() : BufferBase() {

    }

    virtual ~FileBuffer(){
        bclose();
    }

    bool bopen(const std::string & fin){
        fp = fopen(fin.c_str(), "r");
        init();
        if(!fp) {
            eof = true;
            return false;
        }
        return true;
    }
    
    void bclose(){
        if(fp != nullptr) fclose(fp);
        fp = nullptr;
    }

    size_t bread() {
        return fread(&buffer[0], 1, BUFF_SIZE, fp);
    }

    protected:
    FILE      * fp;
};

/*
TODO: Add better error handling
*/
class FileWrapper {
    public:
    FileWrapper(const std::string & fin) : buffer(nullptr){
        open(fin);
    }

    FileWrapper() : buffer(nullptr) {

    }

    ~FileWrapper() {
        if(buffer != nullptr) delete buffer;
        buffer = nullptr;
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
            std::cout << "Could not open " << fin << " for reading\n";
            exit(1);
            //return false;
        }
        return true;
    }

    int get_char(){
        return buffer->get_char();
    }

    int peek_char(){
        return buffer->peek_char();
    }

    template <typename T>
    int get_until(char c, T & out){
        return buffer->get_until(c, out);
    }

    template <typename T>
    int get_until(std::function<bool(char)> & fn, T & out){
        return buffer->get_until(fn, out);
    }

    template <typename T>
    int get_line(T & out){
        out.clear();
        return buffer->get_until('\n', out);
    }

    int tokenize_line(ParserTokens & toks, char delim='\t'){
        return buffer->tokenize_line(toks, delim);
    }

    private:

    BufferBase * buffer;
};

}
