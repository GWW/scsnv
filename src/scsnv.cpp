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

#include <string>
#include <iostream>
#include <getopt.h>
#include "scsnv/build.hpp"
#include "scsnv/index.hpp"
#include "progs/index.hpp"
#include "progs/map.hpp"
#include "progs/barcodes.hpp"
#include "progs/quant.hpp"
#include "progs/collapse.hpp"
#include "progs/pileup.hpp"

using namespace std;
using namespace gwsc;

int main(int argc, char * argv[]){
    if(argc < 2) return 1;
    string cmd = argv[1];
    argc--;
    argv++;
    if(cmd == "index"){
        ProgIndex prog;
        prog.parse(argc, argv);
    }else if(cmd == "map"){
        ProgMap prog;
        prog.parse(argc, argv);
    }else if(cmd == "count"){
        ProgBarcodes prog;
        prog.parse(argc, argv);
    }else if(cmd == "quant"){
        ProgQuant prog;
        prog.parse(argc, argv);
    }else if(cmd == "collapse"){
        ProgCollapse prog;
        prog.parse(argc, argv);
    }else if(cmd == "pileup"){
        ProgPileup prog;
        prog.parse(argc, argv);
    }
    return 0;
}
