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
from .data import build_cols
NP.import_array()

def count_genomes(fname, omap, passed_ids):
    with h5py.File(fname, 'r') as h5f:
        barcode_ids = NP.array(h5f.get('barcode_ids'))
        gene_ids = NP.array(h5f.get('gene_ids'))
        gene_names = NP.array(h5f.get('gene_names'))

        pids = set(passed_ids)
        cidx = NP.fromiter((i for i,b in enumerate(barcode_ids) if b in pids), dtype='uint')
        pids = set(cidx)

        cols, col_idx = build_cols(h5f, True, pids)
        icols, icol_idx = build_cols(h5f, False, pids)

        cmap = NP.zeros(len(barcode_ids), dtype='uint32')
        cids = NP.unique(cols)
        cmap[cids] = NP.arange(0, len(cids))

        rows = NP.array(h5f.get('exonic/rows'))
        rows = rows[col_idx]

        irows = NP.array(h5f.get('intronic/rows'))
        irows = irows[icol_idx]

        gids = NP.unique(NP.concatenate([NP.unique(rows), NP.unique(irows)]))

        gcounts = NP.zeros((4, len(passed_ids)), dtype='int')

        opos = NP.zeros(gids.max() + 1, dtype='int')
        for g in gids:
            opos[g] = omap[gene_ids[g]]

        data = NP.array(h5f.get('exonic/data'))[col_idx]
        for c, r, d in zip(cols, rows, data):
            gcounts[opos[r], cmap[c]] += d

        data = NP.array(h5f.get('intronic/data'))[icol_idx]
        for c, r, d in zip(icols, irows, data):
            gcounts[opos[r] + 2, cmap[c]] += d
    return gcounts

