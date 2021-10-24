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

from collections import Counter
import numpy as NP
import pandas as pd
import re
from ncls import NCLS
from cyvcf2 import VCF

class AnnotateAux(object):
    def __init__(self):
        self.snps = None
        self.g1000 = None
        self.emap = None
        self.repeats = None
        self.captured = None

    def parse_vcfs(self, fnames, names):
        snps = {}
        for si, f in zip(names, fnames):
            keep = 0
            for r in VCF(f):
                if r.FILTER or r.is_indel or len(r.ALT) > 1:
                    continue
                pos = r.POS - 1
                ref = r.REF 
                alt = r.ALT[0]
                chrom = r.CHROM
                keep += 1
                snps.setdefault(chrom, {}).setdefault(pos, []).append((ref, alt, si))
            print(f'Loaded {keep:,} biallelic SNVs from vcf files')
        return snps

    def annotate_vcfs(self, snvs, fnames, names):
        if not fnames:
            snvs['known'] = ''
            snvs['known_af'] = ''
            snvs['known_het'] = ''
            return
        print("Parsing VCF files")
        if self.snps is None:
            self.snps = self.parse_vcfs(fnames, names)
        snps = self.snps
        refs = snvs['ref'].values
        alts = snvs['alt'].values
        print("Annotating with VCF files")
        found = 0
        known, known_af, known_het = [], [], []
        for i, k in enumerate(snvs.index):
            ref = snps.get(k[0], None)
            if ref is None:
                known.append('')
                continue
            xx = ref.get(k[1], [])
            st = '-'
            kk = []

            for x in xx:
                if refs[i] == x[0] and alts[i] == x[1]:
                    kk.append(x[2])
            if kk:
                found += 1

            known.append(','.join(kk))

        snvs['known'] = known
        print(f'Found in the VCF data {found}')

    def annotate_g1000(self, dfo, fg1000):
        if fg1000 is None:
            dfo['g1000'] = False
            return
        if self.g1000 is None:
            print("Parsing 1000 Genomes file")
            df = pd.read_csv(fg1000, sep="\t", dtype={'chrom':object}, names=['chrom', 'pos', 'ref', 'alt'], header=0)
            df.columns = ['chrom', 'pos', 'ref', 'alt']
            df['pos'] -= 1
            g1000 = {}
            for c, p, a in zip(df['chrom'].values, df['pos'].values, df['alt'].values):
                g1000.setdefault(c, set()).add((p, a))
            self.g1000 = g1000

        g1000 = self.g1000
        print("Annotating the SNVs")
        d = []
        for idx, a in zip(dfo.index, dfo['alt']):
            x = g1000.get(idx[0], None)
            f = False
            if x is not None and (idx[1], a) in x:
                f = True
            d.append(f)

        dfo['g1000'] = d
        ss = NP.sum(d)
        print(f'Found in 1000 Genomes: {ss}')

    def annotate_edits(self, dfo, fedits):
        #REDIPortal is pos + 1
        if fedits is None:
            dfo['REDIportal'] = False
            return

        if self.emap is None:
            print("Annotating edits")
            df = pd.read_csv(fedits, sep="\t")
            emap = set()
            for t in df.itertuples():
                tc = t.Region
                if (t.Ref == 'T' and t.Ed == 'C') or (t.Ref == 'A' and t.Ed == 'G'):
                    emap.add((tc, t.Position - 1))
            self.emap = emap
        emap = self.emap

        REDIportal = NP.zeros(len(dfo), dtype=bool)
        for i, (idx, r, a) in enumerate(zip(dfo.index, dfo['strand_ref'], dfo['strand_alt'])):
            if r == 'A' and a == 'G':
                REDIportal[i] = (idx[0], idx[1]) in emap
        print("REDIportal", REDIportal.sum())
        dfo["REDIportal"] = REDIportal


    def annotate_repeats(self, dfo, frep):
        #UCSC Repeat table is zero based
        if frep is None:
            dfo['rep_strand'] = ''
            dfo['rep_name'] = ''
            dfo['rep_class'] = ''
            dfo['rep_family'] = ''
            return

        print("Annotating repeats")
        if self.repeats is None:
            df = pd.read_csv(frep, sep="\t", skiprows=1, header=None, usecols=[5,6,7,9,10,11,12], names=['chrom', 'start', 'end', 'strand', 'rep_name', 'rep_class', 'rep_family'])
            df.sort_values(by=['chrom', 'start', 'end'], inplace=True)
            refs = {}
            rdata = []
            for t in df.itertuples():
                tc = t.chrom
                refs.setdefault(tc, []).append((t.start, t.end, len(rdata)))
                rdata.append((t.strand, t.rep_name, t.rep_class, t.rep_family))

            print("Building NCLS for each reference")
            for k, v in refs.items():
                starts = NP.array([x[0] for x in v])
                ends = NP.array([x[1] for x in v])
                ids = NP.array([x[2] for x in v])
                refs[k] = NCLS(starts, ends, ids)
            self.repeats = rdata, refs

        rdata, refs = self.repeats

        rep_name = [''] * len(dfo)
        rep_class = [''] * len(dfo)
        rep_family = [''] * len(dfo)
        rep_strand = [''] * len(dfo)

        print("Querying each SNV")
        sta, end, idx = NP.zeros(1, dtype='int'), NP.zeros(1, dtype='int'), NP.zeros(1, dtype='int')
        ocount = 0
        for i, x in enumerate(dfo.index):
            if x[0] not in refs:
                continue
            sta[0] = x[1]
            end[0] = x[1] + 1
            _, r_idxs = refs[x[0]].all_overlaps_both(sta, end, idx)
            overlaps = [rdata[j] for j in r_idxs]
            if len(overlaps) > 0:
                rep_strand[i] = ",".join(o[0] for o in overlaps)
                rep_name[i] = ",".join(o[1] for o in overlaps)
                rep_class[i] = ",".join(o[2] for o in overlaps)
                rep_family[i] = ",".join(o[3] for o in overlaps)
                ocount += 1

        dfo['rep_strand'] = NP.array(rep_strand)
        dfo['rep_name'] = NP.array(rep_name)
        dfo['rep_class'] = NP.array(rep_class)
        dfo['rep_family'] = NP.array(rep_family)
        print(f'SNVs with a repeat masker overlap: {ocount}')

    def annotate_captured(self, dfo, capture, flank):
        if capture is None:
            dfo['captured'] = False
            return
        print("Annotating capture regions")
        if self.captured is None:
            refs = {}
            with open(capture) as fp:
                h = fp.readline()
                for l in fp:
                    t = l.rstrip().split('\t', 3)
                    refs.setdefault(t[0], []).append((max(int(t[1]) - flank, 0), int(t[2]) + flank - 1))

            print("Building NCLS indexes")
            for k, v in refs.items():
                starts = NP.array([x[0] for x in v])
                ends = NP.array([x[1] for x in v])
                ids = NP.array(NP.arange(0 ,len(v)))
                refs[k] = NCLS(starts, ends, ids)
            self.captured = refs

        refs = self.captured
        captured = NP.zeros(len(dfo), dtype='bool')
        sta, end, idx = NP.zeros(1, dtype='int'), NP.zeros(1, dtype='int'), NP.zeros(1, dtype='int')
        for i, x in enumerate(dfo.index):
            if x[0] not in refs:
                continue

            sta[0] = x[1]
            end[0] = x[1] + 1
            _, overlaps = refs[x[0]].all_overlaps_both(sta, end, idx)
            captured[i] = len(overlaps) > 0

        dfo['captured'] = captured
        print(f'SNVs with capture region overlap: {captured.sum()}')

