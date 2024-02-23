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

import pandas as pd, numpy as NP, scipy
import scipy.io, h5py
from scipy.sparse import csc_matrix
import flammkuchen as fl
import os, sys
import tables, tarfile, argparse, time
import tarfile

def parse_data(fn, bset):
    with tables.open_file(fn, 'r') as f:
        gkey = next(iter(f.get_node('/')._v_children.keys()))
        print("Genome key ",gkey)
        mat_group = f.get_node(f.root, gkey)
        barcodes = NP.array([b.decode("utf-8") for b in f.get_node(mat_group, 'barcodes').read()])
        genes = f.get_node(mat_group, 'genes').read()
        names = f.get_node(mat_group, 'gene_names').read()
        data = getattr(mat_group, 'data').read()
        indices = getattr(mat_group, 'indices').read()
        indptr = getattr(mat_group, 'indptr').read()
        shape = getattr(mat_group, 'shape').read()

        cmatrix = csc_matrix((data, indices, indptr), shape=shape)
        matrix = cmatrix.astype('uint16')

        cidx = [i for i, b in enumerate(barcodes) if b.split('-')[0] in bset]
        matrix = matrix[:,cidx].toarray()
        return matrix, barcodes[cidx], genes, names, cmatrix[:,cidx], barcodes

def cellranger_cmd(args):
    parser = argparse.ArgumentParser(description='Replace a Cell Ranger filtered_gene_bc_matrices folder with a the barcodes from a text file')
    parser.add_argument('barcodes', help='A text file with a list of barcodes without "-" ie. the passed_barcodes.txt.gz file from the cells command')
    parser.add_argument('raw', help='Path to raw_gene_bc_matrices_h5.h5')
    parser.add_argument('out', help='Folder to write the matrices ie. sample/outs/filtered_gene_bc_matrices/genome')
    args = parser.parse_args(args[1:])

    keep = pd.read_csv(args.barcodes)
    bset = {x:i for i, x in enumerate(keep['barcode'].values)}

    #ts = time.strftime("%Y%m%d-%H%M%S")
    #tf = os.path.join(args.out, f'backup_{ts}.tar.gz')
    #print(f'Backing up {args.out} contents to {tf}')
    #with tarfile.open(tf, "w:gz") as tar:
    #    tar.add(args.out, arcname=os.path.basename(args.out))

    sm, sb, sg, sn, cm, all_barcodes = parse_data(args.raw, bset)

    fdir = args.out

    with open(os.path.join(fdir, 'all_barcodes.tsv'), 'w') as fp:
        for r in all_barcodes:
            b = r.split('-')[0]
            fp.write(f'{b}\n')

    with open(os.path.join(fdir, 'genes.tsv'), 'w') as fp:
        for g, n in zip(sg, sn):
            fp.write('{g}\t{n}\n')

    with open(os.path.join(fdir, 'barcodes.tsv'), 'w') as fp:
        for b in sb:
            fp.write(f'{b}\n')

    import scipy.io
    scipy.io.mmwrite(os.path.join(fdir, 'matrix.mtx'), cm)
    fl.save(os.path.join(fdir, 'expr.h5'), {'barcodes':sb, 'gene_ids':sg, 'gene_names':sn, 'matrix':cm, 'all_barcodes':all_barcodes})
    print("WARNING: Make sure to rename the analysis folder or this will cause problems with velocyto!!!")

