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
#include <array>
#include <string>
#include "sequence.hpp"
namespace gwsc {

class Dust{
    public:
        double calculate(const std::string & r, unsigned int end);

    private:
        std::array<unsigned int, 64> hash;

};

inline double Dust::calculate(const std::string & r, unsigned int end){
    hash.fill(0);
    for(size_t i = 0; i < end - 2; i++){
        if(r[i] == 'N' || r[i + 1] == 'N' || r[i + 2] == 'N'){
            continue;
        }
        uint32_t k = (ADNA4::ltable_[static_cast<unsigned int>(r[i])] << 4) 
                        | (ADNA4::ltable_[static_cast<unsigned int>(r[i + 1])] << 2)
                        | (ADNA4::ltable_[static_cast<unsigned int>(r[i + 2])]);
        hash[k]++;
    }

    double score = 0.0;
    unsigned int t = 0;
    for(size_t i = 0; i < hash.size(); i++){
        t += hash[i];
        score += (hash[i] - 1) * (hash[i]) / 2.0;
    }
    score /= t;
    return score;
}

}
