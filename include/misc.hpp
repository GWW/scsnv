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

#include <vector>
#include <sstream>
#include <cstring>
#include <glob.h>
#include <algorithm>
#include <streambuf>
#include <string>
#include <ostream> 
#include <iostream>

namespace gwsc {

struct ReadStringCmp {
    /**
     * From samtools by Heng Li
     */
    bool operator()(const std::string & s1, const std::string & s2) const {
        const char *a = s1.c_str(), *b = s2.c_str();
        const char *pa = a;
        const char *pb = b;
        while (*pa && *pb) {
                if (isdigit(*pa) && isdigit(*pb)) {
                        long ai, bi;
                        ai = std::strtol(pa, const_cast<char **>(&pa), 10);
                        bi = std::strtol(pb, const_cast<char **>(&pb), 10);
                        if(ai < bi)      return true;
                        else if(ai > bi) return false;
                } else {
                        if (*pa != *pb) break;
                        ++pa; ++pb;
                }
        }
        if (*pa == *pb)
            return (pa-a) < (pb-b);
        return *pa<*pb;
    }
};

inline std::vector<std::string> glob(const std::string& pat){
    glob_t glob_result;
    glob(pat.c_str(),GLOB_TILDE,NULL,&glob_result);
    std::vector<std::string> ret;
    for(unsigned int i=0;i<glob_result.gl_pathc;++i){
        ret.push_back(std::string(glob_result.gl_pathv[i]));
    }
    globfree(&glob_result);
    std::sort(ret.begin(), ret.end());
    return ret;
}

class existing_string_buf : public std::streambuf
{
public:
    // Somehow store a pointer to to_append.
    explicit existing_string_buf(std::string &to_append) : 
        m_to_append(&to_append){}

    virtual int_type overflow (int_type c) {
        if (c != EOF) {
            m_to_append->push_back(c);
        }
        return c;
    }

    virtual std::streamsize xsputn (const char* s, std::streamsize n) {
        m_to_append->insert(m_to_append->end(), s, s + n);                                                                                 
        return n;
    }

private:
    std::string *m_to_append;
};


}
