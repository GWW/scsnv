# cython: c_string_type=str, c_string_encoding=ascii, language_level=3str

#Copyright (c) 2018-2020 Gavin W. Wilson
#
#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.

import h5py
import numpy as NP
import pandas as pd
import cython
cimport numpy as NP
from scipy import sparse
NP.import_array()

# Data for loading things from scmap
def load_barcode_rates(fname, full = False):
    with h5py.File(fname, 'r') as h5f:
        barcodes = NP.array(h5f.get('barcodes'))
        if full:
            brates = h5f['barcode_rates_full']
        else:
            brates = h5f['barcode_rates']

        data = {'barcode':NP.array(h5f.get('barcodes')), 'barcode_id':NP.array(h5f.get('barcode_ids'))}
        for k in brates['field_order']:
            if type(k) == bytes:
                k = k.decode('utf-8')
            if k not in brates:
                continue
            data[k] = NP.array(brates[k])
        df = pd.DataFrame(data)
        if type(df['barcode'].values[0]) == bytes:
            df['barcode'] = df['barcode'].str.decode('utf-8')
        return df

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def _make_matrix(NP.uint32_t[:] rows, NP.uint32_t[:] cols, NP.uint32_t[:] data, unsigned int NR, unsigned int NC):
    cdef NP.ndarray matrix = NP.zeros((NR, NC), dtype='int32')
    cdef Py_ssize_t i
    for i in range(0, len(cols)):
        matrix[rows[i], cols[i]] = data[i]
    return matrix

def build_cols(h5f, exonic, pids):
    cols = None
    if exonic:
        cols = NP.array(h5f.get('exonic/cols'))
    else:
        cols = NP.array(h5f.get('intronic/cols'))
    col_idx = NP.fromiter((i for i, c in enumerate(cols) if c in pids), dtype='uint')
    return cols[col_idx], col_idx

def build_matrix(fname, passed_ids, matrix_type = 'dense'):
    with h5py.File(fname, 'r') as h5f:
        barcode_ids = NP.array(h5f.get('barcode_ids'))
        gene_ids = NP.array(h5f.get('gene_ids'))
        gene_names = NP.array(h5f.get('gene_names'))

        pids = set(passed_ids)
        cidx = NP.fromiter((i for i,b in enumerate(barcode_ids) if b in pids), dtype='uint')
        pids = set(cidx)

        cols, col_idx = build_cols(h5f, True, pids)
        icols, icol_idx = build_cols(h5f, False, pids)

        rows = NP.array(h5f.get('exonic/rows'))
        rows = rows[col_idx]

        irows = NP.array(h5f.get('intronic/rows'))
        irows = irows[icol_idx]

        gids = NP.unique(NP.concatenate([NP.unique(rows), NP.unique(irows)]))

        gmap = NP.zeros(len(gene_ids), dtype='uint32')
        gmap[gids] = NP.arange(0, len(gids))
        cmap = NP.zeros(len(barcode_ids), dtype='uint32')
        cids = NP.unique(cols)
        cmap[cids] = NP.arange(0, len(cids))

        if matrix_type == 'dense':
            data = NP.array(h5f.get('exonic/data'))[col_idx]
            mat = _make_matrix(gmap[rows], cmap[cols], data, len(gids), len(pids))
            data = NP.array(h5f.get('intronic/data'))[icol_idx]
            imat = _make_matrix(gmap[irows], cmap[icols], data, len(gids), len(pids))
            print(mat.shape, imat.shape, len(gids))
        elif matrix_type == 'csr' or matrix_type == 'csc':
            data = NP.array(h5f.get('exonic/data'))[col_idx]
            ridx = gmap[rows]
            cidx = cmap[cols]
            if matrix_type == 'csr':
                mat = sparse.csr_matrix((data, (ridx, cidx)), shape=(len(gids), len(pids)), dtype='int32')
            else:
                mat = sparse.csc_matrix((data, (ridx, cidx)), shape=(len(gids), len(pids)), dtype='int32')

            data = NP.array(h5f.get('intronic/data'))[icol_idx]
            ridx = gmap[irows]
            cidx = cmap[icols]
            if matrix_type == 'csr':
                imat = sparse.csr_matrix((data, (ridx, cidx)), shape=(len(gids), len(pids)), dtype='int32')
            else:
                imat = sparse.csc_matrix((data, (ridx, cidx)), shape=(len(gids), len(pids)), dtype='int32')

        return gene_ids[gids], gene_names[gids], mat, imat
