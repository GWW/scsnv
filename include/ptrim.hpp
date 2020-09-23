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

#include "pbase.hpp"
#include "reader.hpp"
#include <exception>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

namespace gwsc{

class ProgTrim : public ProgBase {
    public:
        argagg::parser parser() const;
        std::string    usage() const;
        void           load();
        int            run();

    private:
        template <typename R> 
        int run_10X_();

        std::string              out_;
        std::string              lib_type_;
        std::vector<std::string> dirs_;
        FastqPairs               fastqs_;
        unsigned int             tag_ = 98;
        unsigned int             downsample_ = 200000000;
};

}

