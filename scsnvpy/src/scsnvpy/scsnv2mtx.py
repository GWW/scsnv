#!/usr/bin/env python
#Copyright (c) 2018-2022 Gavin W. Wilson
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

def go():
    import numpy as NP
    import flammkuchen as fl
    import sys, os
    from scipy.io import mmwrite

    fd = sys.argv[1]
    gt_dir = sys.argv[2]

    dd = fl.load(fd)
    barcodes = dd['barcodes']
    if type(barcodes[0]) == bytes:
        barcodes = NP.array([x.decode('utf-8') for x in barcodes])
    snvs = dd['ann']
    strand_ref = dd['strand_ref_mat']
    strand_alt = dd['strand_alt_mat']


    print("Writing ref matrix")
    mmwrite(os.path.join(gt_dir, f'refs.mtx'), strand_ref)

    print("Writing alt matrix")
    mmwrite(os.path.join(gt_dir, f'alts.mtx'), strand_alt)

    print("Writing barcodes")
    fout = open(f'{gt_dir}/barcodes.txt', 'w')
    for b in barcodes:
        fout.write(b + "\n")
    fout.close()

    print('Writing the SNV data')
    snvs.to_csv(f'{gt_dir}/snvs.csv')
