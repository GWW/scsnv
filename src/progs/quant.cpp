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

#include "quant.hpp"
#include "../scmap/quant_worker.hpp"
#include "../scmap/reader.hpp"

using namespace gwsc;

int ProgQuant::run() {
    QuantBase base;
    unsigned int umi_len = 0;
    if(lib_type_ == "V2")      umi_len = Reader10X_V2::UMI_LEN;
    else if(lib_type_ == "V3") umi_len = Reader10X_V3::UMI_LEN;
    else {
        std::cerr << "unsupported library type\n";
        exit(1);
    }
    tout << "Loading the index\n";
    base.load_index(tx_idx_);
    tout << "Preparing the tag index(es)\n";
    base.prepare_tags(prefixes_);
    base.build_gene_groups(gene_groups_);
    base.run(umi_len, threads_, min_molecules_, bam_, full_cmd_);
    base.write_output(out_prefix_, min_molecules_);


    return EXIT_SUCCESS;
}
