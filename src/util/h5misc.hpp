#pragma once
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

#include <H5Cpp.h>
#include <vector>

template <typename T, typename F, typename P>
void write_h5_numeric(const std::string & name, const std::vector<T> & data, F & h5, P dtype){
    using namespace H5;
    //tout << "Writing samples\n";
    hsize_t dimsf[1];
    dimsf[0] = data.size();
    DataSpace dataspace(1, dimsf);
    //
    IntType datatype(dtype);
    datatype.setOrder( H5T_ORDER_LE );
    DSetCreatPropList ds_creatplist;  // create dataset creation prop list
    ds_creatplist.setChunk( 1, dimsf );  // then modify it for compression
    ds_creatplist.setDeflate( 6 );
    DataSet count_dataset = h5.createDataSet(name, datatype, dataspace, ds_creatplist);
    count_dataset.write(data.data(), dtype);
}

template <typename F>
void write_h5_string(const std::string & name, const std::vector<const char *> & data, F & h5){
    using namespace H5;
    StrType st(H5::PredType::C_S1, H5T_VARIABLE);
    st.setCset(H5T_CSET_UTF8);
    hsize_t dim[1] = {data.size()};
    H5::DataSpace ds = H5::DataSpace(1, dim);
    DSetCreatPropList ds_creatplist;  // create dataset creation prop list
    ds_creatplist.setChunk( 1, dim );  // then modify it for compression
    ds_creatplist.setDeflate( 6 );
    h5.createDataSet(name,st,ds, ds_creatplist).write(data.data(), st);
}

