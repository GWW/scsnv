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

import pandas as pd
import numpy as NP
import h5py, sys
import flammkuchen as flm
import argparse
from scipy.sparse import csr_matrix, coo_matrix
from .snvmats import align_matrices, merge_coo_dups, merge_coo_dups_to_csr
#from .snvcmp import snv_type_names, snv_type_funcs, snv_type_shorts

BMAP = {x:i for i, x in enumerate('ACGTN')}
RMAP = {x:y for x, y in zip('ACGTN', 'TGCAN')}
BORDER = ["coding", "non_coding", "utr3", "utr5", "intron", "intergenic"] 

from .annaux import AnnotateAux

def annotate_cmd(args):
    parser = argparse.ArgumentParser(description='Collate and annotate scSNV mutation data')
    parser.add_argument('-e', '--edits', help='RNA Editing database', type=str)
    parser.add_argument('-r', '--repeats', help='UCSC Repeat Masker Table', type=str)
    parser.add_argument('-d', '--g1000', help='Biallelic 1000 genomes processed text file (see manual)', type=str)
    parser.add_argument('-v', '--vcf', help='One or more vcf files to use for mutation annotation (space separated)', action='append', type=str)
    parser.add_argument('-n', '--name', help='Aliases for the vcf files. Must been in the same order as --vcf', action='append', type=str)
    parser.add_argument('-u', '--purity', help='Proportion of non-reference reads that must match a single allele', type=float, default=0.95)
    parser.add_argument('-s', '--strand-purity', help='Proportion of reads that must be derived from the same strand to call a stranded SNP', type=float, default=0.9)
    parser.add_argument('-c', '--capture', help='Capture regions if VCF files are whole exome', type=str)
    parser.add_argument('-f', '--flank', help='Capture regions flank', type=int, default=100)
    #parser.add_argument('-t', '--targets', help='Get ALTS from this target list', type=str)
    parser.add_argument('prefix', help='The pileup file prefix(es) generated by scsnv pileup', type=str, nargs='+')
    args = parser.parse_args(args[1:])

    aux = AnnotateAux()
    if args.vcf is not None and len(args.vcf) != len(args.name):
        print("Error there must be a name for each vcf file")
        sys.exit()

    for prefix in args.prefix:
        print(f'Processing {prefix}')
        ann = SNVAnnotate(prefix, args.purity, args.strand_purity, None)

        aux.annotate_vcfs(ann._snvs, args.vcf, args.name)
        aux.annotate_g1000(ann._snvs, args.g1000)
        aux.annotate_repeats(ann._snvs, args.repeats)
        aux.annotate_captured(ann._snvs, args.capture, args.flank)
        aux.annotate_edits(ann._snvs, args.edits)

        ann._snvs.to_csv(prefix + '_passed_snvs.txt.gz', sep="\t")
        flm.save(prefix + '_annotated.h5', {
            'barcodes':ann._barcodes, 'ann':ann._snvs, 
            'strand_ref_mat':ann.smats[0], 'strand_alt_mat':ann.smats[1],
            'filtering':{'kept':(ann._ppass & ann._spass).sum(), 'total':len(ann._ppass)}
        })


class SNVAnnotate(object):
    def __init__(self, prefix, purity, strand_purity, targets=None):
        df = pd.read_csv(prefix + ".txt.gz", dtype={'chrom':object}, sep='\t', index_col=[0, 1])
        #self._trates = pd.read_csv(prefix + '_base_rates.txt')
        h5f = h5py.File(prefix + '_barcode_matrices.h5', 'r')
        self._barcodes = NP.array(h5f['barcodes'])
        self._raw_snvs = df
        self._h5f = h5f
        self._prefix = prefix
        #self._check_cell_counts()
        if targets is not None:
            tmap = {}
            with open(targets) as fp:
                fp.readline()
                for l in fp:
                    t = l.rstrip().split("\t")
                    tmap[(t[0], int(t[1]))] = t[2]
            fixed = 0
            for i, idx in enumerate(df.index):
                na = tmap.get(idx)
                if na is not None:
                    fixed += 1
                    df['max_non_ref'].values[i] = na
                else:
                    print("Missing ", idx)
            print(f'Fixed {fixed} targets out of {len(tmap)}')


        self._call_raw_snvs(purity, strand_purity)
        self._build_matrices()
        self._filter_bad()

    def _filter_bad(self):
        keep = self._ppass & self._spass
        c1 = (self._ppass & ~self._spass)
        c2 = (~self._ppass & self._spass)
        c3 = (self._ppass & self._spass)
        fout = open(f'{self._prefix}_filtering.txt', 'w')
        fout.write(f'total\t{keep}\n')
        fout.write(f'biallelic_filtered\t{c1}\n')
        fout.write(f'strand_specific_filtered\t{c2}\n')
        fout.write(f'both_filtered\t{c3}\n')
        fout.write(f'total_passed\t{keep.sum()}\n')
        fout.close()
        print(f'Kept {keep.sum()} out of {len(keep)}')
        self._snvs = self._snvs[keep]
        #for i, m in enumerate(self.mats):
        #    self.mats[i] = csr_matrix(m[keep])

        #for i, m in enumerate(self.smats):
        #    self.smats[i] = csr_matrix(m[keep])
        print(keep, keep.shape, self.smats[0].shape, self.smats[1].shape)
        self.smats[0] = self.smats[0][keep]
        self.smats[1] = self.smats[1][keep]
        #indptr, indices, ref, alt = align_matrices(self.smats[0], self.smats[1])
        #self.smats[0] = csr_matrix((ref,indices,indptr), shape=(self.smats[0].shape[0],self.smats[0].shape[1]))
        #self.smats[1] = csr_matrix((alt,indices,indptr), shape=(self.smats[0].shape[0],self.smats[0].shape[1]))


    #Check to make sure the cell counts are correct
    def _check_cell_counts(self):
        for bi, b in enumerate("ACGT"):
            d = self._h5f[f'base_{b}']
            bids, sids, pbases, mbases = map(NP.array, (d['barcode_ids'], d['snps'], d['plus'], d['minus']))
            ctot = f'{b}_total'
            ctbc = f'{b}_total_barcodes'
            cp = f'{b}_plus_counts'
            cm = f'{b}_minus_counts'

            for i, r in enumerate(self._raw_snvs.itertuples()):
                vp = getattr(r, cp)
                vm = getattr(r, cm)
                vtot = getattr(r, ctot)
                vtbc = getattr(r, ctbc)

                ridx = NP.where(sids == i)[0]
                pc = pbases[ridx]
                mc = mbases[ridx]

                mp = pc.sum()
                mm = mc.sum()
                mtot = mp + mm
                mtbc = set(bids[ridx])
                if vp != mp or vm != mm or mtot != vtot or len(mtbc) != vtbc:
                    print("Error")

    def _call_raw_snvs(self, min_purity, min_spurity):
        snvs = self._raw_snvs
        plus_rcounts = NP.zeros(len(snvs), dtype='int32')
        plus_acounts = NP.zeros(len(snvs), dtype='int32')
        minus_rcounts = NP.zeros(len(snvs), dtype='int32')
        minus_acounts = NP.zeros(len(snvs), dtype='int32')
        for b in 'ACGT':
            idx = (snvs.ref == b).values
            plus_rcounts[idx] = snvs[f'{b}_plus_counts'].values[idx]
            minus_rcounts[idx] = snvs[f'{b}_minus_counts'].values[idx]

            idx = (snvs.max_non_ref == b).values
            plus_acounts[idx] = snvs[f'{b}_plus_counts'].values[idx]
            minus_acounts[idx] = snvs[f'{b}_minus_counts'].values[idx]


        tot = snvs[[f'{b}_total' for b in 'ACGT']].sum(axis=1).values
        rcounts = plus_rcounts + minus_rcounts
        acounts = plus_acounts + minus_acounts

        with NP.errstate(divide='ignore'):
            purity = (plus_rcounts + plus_acounts + minus_rcounts + minus_acounts) / tot
            purity[NP.isnan(purity)] = -1

        ppass = (purity >= min_purity)
        self._ppass = ppass

        with NP.errstate(divide='ignore'):
            spurity = (plus_rcounts + plus_acounts) / (rcounts + acounts)

        pstrand = ~pd.isnull(spurity) & (spurity >= min_spurity)
        mstrand = ~pd.isnull(spurity) & (spurity <= (1 - min_spurity ))
        spass = pstrand | mstrand
        self._spass = spass

        dfo = pd.DataFrame(index=snvs.index)
        dfo['ref'] = snvs.ref
        dfo['alt'] = '-'
        dfo.loc[ppass, 'alt'] = snvs.max_non_ref
        dfo['ref_count'] = rcounts
        dfo['ref_plus_count'] = plus_rcounts
        dfo['ref_minus_count'] = minus_rcounts

        dfo['alt_count'] = 0
        dfo.loc[ppass, 'alt_count'] = acounts[ppass]
        dfo.loc[ppass, 'alt_plus_count'] = plus_acounts[ppass]
        dfo.loc[ppass, 'alt_minus_count'] = minus_acounts[ppass]


        dfo['plus_base'] = snvs['plus_base']
        dfo['minus_base'] = snvs['minus_base']

        dfo['strand_ref'] = '-'
        dfo['strand_alt'] = '-'
        dfo['strand'] = '?'
        dfo.loc[pstrand, 'strand'] = '+'
        dfo.loc[mstrand, 'strand'] = '-'
        dfo.loc[pstrand, 'strand_ref'] = snvs.ref[pstrand]
        dfo.loc[pstrand, 'strand_alt'] = snvs.max_non_ref[pstrand]
        dfo.loc[mstrand, 'strand_ref'] = [RMAP[x] for x in snvs.ref[mstrand]]
        dfo.loc[mstrand, 'strand_alt'] = [RMAP[x] for x in snvs.max_non_ref[mstrand]]


        dfo.loc[pstrand, 'donor_dist'] = snvs.plus_donor_dist[pstrand]
        dfo.loc[pstrand, 'acceptor_dist'] = snvs.plus_acceptor_dist[pstrand]
        dfo.loc[mstrand, 'donor_dist'] = snvs.minus_donor_dist[mstrand]
        dfo.loc[mstrand, 'acceptor_dist'] = snvs.minus_acceptor_dist[mstrand]

        dfo['ref_strand_count'] = 0
        dfo['alt_strand_count'] = 0
        dfo.loc[pstrand, 'ref_strand_count'] = plus_rcounts[pstrand]
        dfo.loc[mstrand, 'ref_strand_count'] = minus_rcounts[mstrand]
        dfo.loc[pstrand, 'alt_strand_count'] = plus_acounts[pstrand]
        dfo.loc[mstrand, 'alt_strand_count'] = minus_acounts[mstrand]

        dfo['strand_AF'] = NP.nan
        dfo.loc[spass, 'strand_AF'] = 1.0 * dfo['alt_strand_count'] / (dfo['alt_strand_count'] + dfo['ref_strand_count'])
        dfo['strand_base'] = ''
        dfo.loc[pstrand, 'strand_base'] = dfo.loc[pstrand, 'plus_base']
        dfo.loc[mstrand, 'strand_base'] = dfo.loc[mstrand, 'minus_base']

        #print((purity >= min_purity).sum(), len(purity))
        #print(((spurity >= min_spurity) | (spurity <= (1 - min_spurity))).sum(), len(purity))
        self._snvs = dfo

    def _build_matrices(self):
        snvs = self._snvs
        all_barcodes = self._barcodes
        BC = len(all_barcodes)
        #strand_ref_mat = NP.zeros((len(snvs), BC), dtype='int32')
        #strand_alt_mat = NP.zeros((len(snvs), BC), dtype='int32')

        strand_ref_mat = None
        strand_alt_mat = None

        cmat = None
        test = (self._ppass & self._spass)
        dpos = NP.array(self._h5f['pos'])

        for b1 in 'ACGT':
            d = self._h5f[f'base_{b1}']
            bids1, sids1, pbases1, mbases1 = map(NP.array, (d['barcode_ids'], d['snps'], d['plus'], d['minus']))
            for b2 in 'ACGT':
                if b1 == b2:
                    continue
                d = self._h5f[f'base_{b2}']
                bids2, sids2, pbases2, mbases2 = map(NP.array, (d['barcode_ids'], d['snps'], d['plus'], d['minus']))
                base_idx = (snvs['ref'].values == b1) & (snvs['alt'].values == b2)

                for strand in '+-':
                    strand_idx = (base_idx) & (snvs['strand'].values == strand)
                    if ~NP.any(strand_idx):
                        continue
                    c1 = pbases1 if strand == '+' else mbases1
                    c2 = pbases2 if strand == '+' else mbases2
                    idx = NP.where(strand_idx)[0]
                    ridx1 = NP.where(NP.in1d(sids1, idx))[0]
                    ridx2 = NP.where(NP.in1d(sids2, idx))[0]

                    rmat = NP.zeros((len(ridx1) + len(ridx2), 3), dtype='int32')
                    rmat[:len(ridx1),0] = sids1[ridx1]
                    rmat[:len(ridx1),1] = bids1[ridx1]
                    rmat[:len(ridx1),2] = c1[ridx1]

                    amat = NP.zeros((len(ridx1) + len(ridx2), 3), dtype='int32')
                    amat[:len(ridx2),0] = sids2[ridx2]
                    amat[:len(ridx2),1] = bids2[ridx2]
                    amat[:len(ridx2),2] = c2[ridx2]



                    rmat[len(ridx1):,2] = 0
                    rmat[len(ridx1):,0] = amat[:len(ridx2),0]
                    rmat[len(ridx1):,1] = amat[:len(ridx2),1]

                    amat[len(ridx2):,0] = rmat[:len(ridx1),0]
                    amat[len(ridx2):,1] = rmat[:len(ridx1),1]
                    amat[len(ridx2):,2] = 0

                    rmat = rmat[NP.lexsort((rmat[:,2], rmat[:,1], rmat[:,0]))]
                    amat = amat[NP.lexsort((amat[:,2], amat[:,1], amat[:,0]))]

                    rmat = merge_coo_dups(rmat)
                    amat = merge_coo_dups(amat)

                    tmat = NP.zeros((amat.shape[0], 4), dtype='int32')
                    tmat[:,0] = amat[:,0]
                    tmat[:,1] = amat[:,1]
                    tmat[:,2] = rmat[:,2]
                    tmat[:,3] = amat[:,2]
                    tmat = tmat[tmat[:,2:].sum(axis=1) > 0]
                    if cmat is None:
                        cmat = tmat
                    else:
                        cmat = NP.vstack([cmat, tmat])

        self.smats = [
            csr_matrix((cmat[:, 2], (cmat[:, 0], cmat[:, 1])), shape=(len(snvs), BC)),
            csr_matrix((cmat[:, 3], (cmat[:, 0], cmat[:, 1])), shape=(len(snvs), BC))
        ]

        tb = NP.zeros(len(snvs), dtype='int')
        bb = NP.zeros(len(snvs), dtype='int')
        sa = NP.zeros(len(snvs), dtype='int')
        sr = NP.zeros(len(snvs), dtype='int')
        c_rmat, c_amat = self.smats

        for i, (s, e) in enumerate(zip(c_rmat.indptr, c_rmat.indptr[1:])):
            i1 = c_rmat.data[s:e] > 0
            i2 = c_amat.data[s:e] > 0
            tb[i] += (i1 | i2).sum()
            bb[i] += (i1 & i2).sum()
            sr[i] += (i1 & ~i2).sum()
            sa[i] += (i2 & ~i1).sum()

        self._snvs['total_barcodes'] = tb
        self._snvs['ref_strand_barcodes'] = sr
        self._snvs['alt_strand_barcodes'] = sa
        self._snvs['refalt_strand_barcodes'] = bb
