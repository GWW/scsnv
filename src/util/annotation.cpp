/*
Copyright (c) 2018-2019 Gavin W. Wilson

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

#include "annotation.hpp"
#include <algorithm>
#include <iostream>
#include <numeric>

using namespace std;
using namespace gwsc;

void gwsc::sort_model(GeneModel & m){
    for(auto & c : m.chroms){
        m.chrom_order.push_back(c.first);
        for(auto & g : c.second){
            g.lft = g.children.front().children.front().lft;
            g.rgt = g.children.front().children.front().rgt;
            for(auto & t : g.children){
                std::sort(t.children.begin(), t.children.end());
                int tlft = 0;
                for(auto & e : t.children){
                    e.tlft = tlft;
                    e.trgt = tlft + e.size() - 1;
                    tlft = e.trgt + 1;
                }
                t.lft = t.children.front().lft;
                t.rgt = t.children.back().rgt;
                g.lft = min(t.lft, g.lft);
                g.rgt = max(t.rgt, g.rgt);
            }
            std::sort(g.children.begin(), g.children.end());
        }
        std::sort(c.second.begin(), c.second.end());
    }
    std::sort(m.chrom_order.begin(), m.chrom_order.end(), ReadStringCmp());
}

