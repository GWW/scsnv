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

import numpy as NP
import pandas as pd
import flammkuchen as fl
from collections import Counter
from itertools import permutations
from flammkuchen.hdf5io import _HDFStoreWithHandle
from pandas.api.types import CategoricalDtype
import os, tables

border = ['5', 'C', 'E', '3', 'N', 'I']
bshort = ['utr5', 'coding', 'non_coding', 'utr3', 'intron', 'intergenic']
blong = ['5\'-UTR', 'Exon' ,'ncExon', '3\'-UTR', 'Intronic', 'Intergenic']
bcoding = [0, 1, 0, 0, 0, 0]
bmap = {s:n for s, n in zip(bshort, blong)}

#snv_type_names = ['Annotated', 'A-to-G Edits', 'Unannotated Transitions', 'Unannotated Transversions']
#snv_type_shorts = ['annotated', 'A_to_G_edits', 'unannotated_transitions', 'unannotated_transversions']
#snv_type_map = {s:n for n, s in zip(snv_type_names, snv_type_shorts)}

_edit = (lambda df: 
               (df.strand_ref == 'A') & (df.strand_alt == 'G') & \
               (df.REDIportal | df.rep_family.str.contains('Alu')))

def _is_transition(df):
    TM = set(('AG', 'GA', 'CT', 'TC'))
    return df.strand_change.isin(TM)

def _is_transversion(df):
    TV = set(('AC', 'CA', 'GT', 'TG', 'CG', 'GC', 'AT', 'TA'))
    return df.strand_change.isin(TV)

_snv_type_data = {
    'g1000':('1000 Genomes', (lambda df:df.g1000)),
    'A_to_G_edit':('A-to-G Edit', _edit),
    'germline':('Germline', (lambda df: (df.known.str.contains('germline') | df.known.str.contains('WGS')))),
    'somatic':('Somatic', (lambda df: df.known.str.contains('somatic'))),
    'captured':('Exome Captured', (lambda df:df.captured)),
    'raw':('Unfiltered Germline', (lambda df: df.known.str.contains('raw'))),
    'transition':('Transition', _is_transition),
    'transversion':('Transversion', _is_transversion),
}

class GroupBuilder(object):
    def __init__(self, captured=False):
        self._captured = captured
        self._groups = []
        self._names = []
        self._order = []
        self._vals = []

    def add_group(self, short, name, *groups):
        self._names.append(name)
        self._order.append(short)
        negs = []
        grps = []
        for g in groups:
            if g[0] == '~':
                grps.append('ann_' + g[1:])
                negs.append(False)
            else:
                grps.append('ann_' + g)
                negs.append(True)
        self._vals.append(negs)
        self._groups.append(grps)

    def apply(self, df):
        if self._captured:
            self._calc(df, True, 'group')
            self._calc(df, False, 'raw_group')
        else:
            self._calc(df, False, 'group')

    def _calc(self, df, captured, out):
        sidx = None
        if captured:
            sidx = df['ann_captured'].values
        else:
            sidx = NP.full(len(df), True, dtype='bool')

        for s, grps, vals in zip(self._order, self._groups, self._vals):
            idx = sidx.copy()
            for g, n in zip(grps, vals):
                res = df[g] == n
                idx &= (df[g] == n)
            df[s] = idx

        kvals = []
        for i, r in enumerate(df[self._order].itertuples(name=None)):
            idx = r[0]
            vals = r[1:]
            if NP.any(vals):
                kvals.append(','.join(self._order[i] for i, v in enumerate(vals) if v))
            else:
                kvals.append('')

        df[out] = kvals

    def __len__(self):
        return len(self._order)

    @property
    def shorts(self):
        return self._order

    @property
    def names(self):
        return self._names

snv_full_names = {k:v[0] for k, v in _snv_type_data.items()}

sub_order = sorted([c for c in permutations('ACGT', 2)])
sub_names = sorted([c[0] + c[1] for c in permutations('ACGT', 2)])

class SampleData(object):
    def __init__(self, fd, name, load_matrix = True, keep_barcodes = None, groups=None):
        self.name = name
        if load_matrix:
            dd = fl.load(fd)
            self.barcodes = dd['barcodes']
            self.snvs = dd['ann']
            self.strand_ref = dd['strand_ref_mat']
            self.strand_alt = dd['strand_alt_mat']
            if keep_barcodes is not None:
                cidx = [i for i, b in enumerate(self.barcodes) if b in keep_barcodes]
                print(f'Keeping {len(cidx)} barcodes out of {len(self.barcodes)}')
                self.strand_ref = self.strand_ref[:,cidx]
                self.strand_alt = self.strand_alt[:,cidx]
                self.barcodes = self.barcodes[cidx]

        else:
            with tables.open_file(fd, mode='r') as fp:
                hp = _HDFStoreWithHandle(fp)
                self.snvs = hp.get('ann')

        self.snvs['snv_idx'] = NP.arange(0, len(self.snvs))
        self.snvs['g1000_raw'] = self.snvs['g1000']
        self.snvs['strand_change'] = self.snvs['strand_ref'] + self.snvs['strand_alt']

        if 'known' not in self.snvs.columns:
            self.snvs['known'] = ''
    
        self._build_annotations(groups)

    @property
    def strand_AF(self):
        return self.strand_alt.toarray() / (self.strand_ref + self.strand_alt).toarray()

    def _build_annotations(self, groups):
        name_map = {}
        adf = pd.DataFrame(index=self.snvs.index)
        for k, v in _snv_type_data.items():
            idx = v[1](self.snvs)
            self.snvs['ann_'+k] = v[1](self.snvs)
        if groups is not None:
            groups.apply(self.snvs)



def build_snv_summaries_captured(groups, **kwargs):
    NG = len(groups) * 3
    summary_mat = NP.zeros((NG, len(kwargs) * 3), dtype='int')
    ssummary_mat = NP.zeros((NG * len(sub_order), len(kwargs) * 3), dtype='int')
    bsummary_mat = NP.zeros((NG * len(border), len(kwargs) * 3), dtype='int')


    for i, (name, mdata) in enumerate(kwargs.items()):
        snvs = mdata.merged
        for ci, captured in enumerate((1, 0, 2)):
            for ji, g in enumerate(groups.shorts):
                j = ji + len(groups) * ci
                if captured < 2:
                    df = snvs[snvs[g] & (snvs['captured'] == captured)]
                else:
                    df = snvs[snvs[g]]
                s1 = df['scSNV_alt_strand_count'] > 0
                s2 = df['CR_alt_strand_count'] > 0

                ts1, ts2, tov = s1.sum(), s2.sum(), (s1 & s2).sum()
                summary_mat[j, i * 3:(i + 1) * 3] = [ts1, ts2, tov]

                d = {}
                for idx, g in df.groupby('strand_change'):
                    s1 = (g['scSNV_alt_strand_count'] > 0).values
                    s2 = (g['CR_alt_strand_count'] > 0).values
                    d[idx] = (s1.sum(), s2.sum(), (s1 & s2).sum())

                for bi, b in enumerate(sub_names):
                    ssummary_mat[j * len(sub_order) + bi, i * 3:(i + 1) * 3] = d.get(b, 0)

                d = {}
                for idx, g in df.groupby('strand_base'):
                    s1 = (g['scSNV_alt_strand_count'] > 0).values
                    s2 = (g['CR_alt_strand_count'] > 0).values
                    d[idx] = (s1.sum(), s2.sum(), (s1 & s2).sum())

                for bi, b in enumerate(border):
                    bsummary_mat[j * len(border) + bi, i * 3:(i + 1) * 3] = d.get(b, 0)

    col_names = sum(([f'{k}_scSNV', f'{k}_CR', f'{k}_overlap'] for k in kwargs.keys()), [])
    gnames = [f'captured_{x}' for x in groups.shorts] + [f'outside_{x}' for x in groups.shorts] + [f'total_{x}' for x in groups.shorts]

    df_summary = pd.DataFrame(summary_mat, columns=col_names, index=gnames)
    df_summary.index.name = 'SNV_type'

    srow_names = sum(([s] * len(sub_names) for s in gnames), [])
    df_sub_summary = pd.DataFrame(ssummary_mat, columns=col_names, index=srow_names)
    df_sub_summary.insert(0, 'substitution', sub_names * NG)
    df_sub_summary.index.name = 'SNV_type'

    brow_names = sum(([s] * len(border) for s in gnames), [])
    df_base_summary = pd.DataFrame(bsummary_mat, columns=col_names, index=brow_names)
    df_base_summary.insert(0, 'base', bshort * NG)
    df_base_summary.index.name = 'SNV_type'

    return df_summary, df_sub_summary, df_base_summary

def build_snv_summaries(groups, **kwargs):
    summary_mat = NP.zeros((len(groups), len(kwargs) * 3), dtype='int')
    ssummary_mat = NP.zeros((len(groups) * len(sub_order), len(kwargs) * 3), dtype='int')
    bsummary_mat = NP.zeros((len(groups) * len(border), len(kwargs) * 3), dtype='int')


    for i, (name, mdata) in enumerate(kwargs.items()):
        snvs = mdata.merged
        for j, g in enumerate(groups.shorts):
            df = snvs[snvs[g]]
            s1 = df['scSNV_alt_strand_count'] > 0
            s2 = df['CR_alt_strand_count'] > 0

            ts1, ts2, tov = s1.sum(), s2.sum(), (s1 & s2).sum()
            summary_mat[j, i * 3:(i + 1) * 3] = [ts1, ts2, tov]

            d = {}
            for idx, g in df.groupby('strand_change'):
                s1 = (g['scSNV_alt_strand_count'] > 0).values
                s2 = (g['CR_alt_strand_count'] > 0).values
                d[idx] = (s1.sum(), s2.sum(), (s1 & s2).sum())

            for bi, b in enumerate(sub_names):
                ssummary_mat[j * len(sub_order) + bi, i * 3:(i + 1) * 3] = d.get(b, 0)

            d = {}
            for idx, g in df.groupby('strand_base'):
                s1 = (g['scSNV_alt_strand_count'] > 0).values
                s2 = (g['CR_alt_strand_count'] > 0).values
                d[idx] = (s1.sum(), s2.sum(), (s1 & s2).sum())

            for bi, b in enumerate(border):
                bsummary_mat[j * len(border) + bi, i * 3:(i + 1) * 3] = d.get(b, 0)

    col_names = sum(([f'{k}_scSNV', f'{k}_CR', f'{k}_overlap'] for k in kwargs.keys()), [])

    df_summary = pd.DataFrame(summary_mat, columns=col_names, index=groups.shorts)
    df_summary.index.name = 'SNV_type'

    srow_names = sum(([s] * len(sub_names) for s in groups.shorts), [])
    df_sub_summary = pd.DataFrame(ssummary_mat, columns=col_names, index=srow_names)
    df_sub_summary.insert(0, 'substitution', sub_names * len(groups))
    df_sub_summary.index.name = 'SNV_type'

    brow_names = sum(([s] * len(border) for s in groups.shorts), [])
    df_base_summary = pd.DataFrame(bsummary_mat, columns=col_names, index=brow_names)
    df_base_summary.insert(0, 'base', bshort * len(groups))
    df_base_summary.index.name = 'SNV_type'

    return df_summary, df_sub_summary, df_base_summary
