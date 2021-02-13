import numpy as NP
import glob
import pandas as pd
from ncls import NCLS
from cyvcf2 import VCF
from collections import Counter
from .snvcmp import SampleData, _snv_type_data
import h5py, sys
from .annaux import AnnotateAux

def merge_samples(samples, groups):
    cols = ['ref', 'alt', 'known', 'g1000', 'g1000_raw', 'captured', 'strand', 'strand_ref', 'strand_alt', 'strand_base', 'strand_change', 'donor_dist', 'acceptor_dist',
            'REDIportal', 'rep_strand', 'rep_name', 'rep_class', 'rep_family']
    xcols = cols[:] + ['idx_tmp']

    mcols = [
       'ref_strand_count', 'alt_strand_count', 'strand_AF', 
       'total_barcodes', 
       'ref_strand_barcodes',
       'alt_strand_barcodes', 'snv_idx'
   ]

    print(f'Merging {len(samples)} Dataframes Annotation Data')
    index = {}
    invalid = set()
    merged = None
    for j, sd in enumerate(samples):
        s = sd.snvs
        s['idx_tmp'] = s.index
        if merged is None:
            merged = s[xcols].copy()
            print(f'  {sd.name} First Sample {len(merged)}')
        else:
            #merged = pd.merge(merged, s[['idx_tmp', 'strand_alt']], on=['idx_tmp', 'strand_alt'], how='outer')
            merged = pd.concat([merged,s[xcols]])
            X = len(merged)
            merged.drop_duplicates(subset=['idx_tmp', 'strand_alt'], inplace=True, keep='last')
            print(f'  {sd.name} {X} after {len(merged)}')

        print(f'    g1000 --> {merged.g1000_raw.sum():,}')
    merged.drop_duplicates(subset=['idx_tmp'], inplace=True, keep='last')
    print('After removing bad rows', len(merged))

    merged.set_index(pd.MultiIndex.from_tuples(merged.idx_tmp, names=['chrom', 'pos']), inplace=True)
    del merged['idx_tmp']
    print(f'Merged SNVs: {len(merged)}')
    for s in samples:
        merged = merged.merge(s.snvs[mcols].add_prefix(s.name + '_'), how='left', left_index=True, right_index=True)

    print("Cleaning up data")
    for f in merged.columns:
        if 'snv_idx' in f:
            merged[f].fillna(-1, inplace=True)
    merged.fillna(0, inplace=True, downcast='infer')

    print("Adding annotations")
    for k, v in _snv_type_data.items():
        idx = v[1](merged)
        merged['ann_' + k] = v[1](merged)
        ss =  merged['ann_' + k].sum()
        print(f'  {k} --> {ss:,}')

    groups.apply(merged)

    return merged

class SNVData(object):
    #df = pd.DataFrame(rows, columns=['chrom', 'pos', 'ref', 'alt', 'strand', 'genome'] + names)
    chrom = ''
    pos = ''
    ref = ''
    alt = ''
    strand = ''
    extra = None
    genome = 0
    valid = False

    def __init__(self, chrom, pos, ref, alt, strand):
        self.extra = {}
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.strand = strand
        self.genome = 0

    def key(self):
        return (self.chrom, self.pos)


def merge_with_genomes_cmd(args):
    import sys, argparse, os
    from .snvcmp import GroupBuilder
    parser = argparse.ArgumentParser(description='Merge SNVs from multiple VCFs and annotated SNV samples')
    parser.add_argument('-g', '--genome', help='Genome VCF file', action='append', type=str, required=True)
    parser.add_argument('-n', '--names', help='Key names for each dsc pileup', action='append', type=str, required=True)
    parser.add_argument('-c', '--capture', help='Use capture region annotations to limit annotated SNVs (Assumes VCF files are already filtered)', type=str)
    parser.add_argument('-f', '--flank', help='Capture regions flank', type=int, default=100)
    parser.add_argument('-o', '--output', help='Output text file (can be gz)', type=str, required=True)
    parser.add_argument('pileups', help='Annotated scsnv pileup prefix ie. sample/pileup for pileup_barcode_matrices.h5 and pileup_annotated.h5', type=str, nargs='+')
    args = parser.parse_args(args[1:])

    names = args.names
    pileups = args.pileups;

    groups = GroupBuilder(captured=False)
    groups.add_group('germline', 'Germline', 'germline')
    groups.add_group('g1000', 'Common', 'g1000', '~germline')
    groups.add_group('edit', 'A-to-G Edit', 'A_to_G_edit', '~germline', '~g1000')
    groups.add_group('transition', 'Transition', 'transition', '~A_to_G_edit', '~germline', '~g1000')
    groups.add_group('transversion', 'Transversion', 'transversion', '~A_to_G_edit', '~germline', '~g1000')

    if len(names) != len(pileups):
        print("number of names is not equal to the number of pileups")
        sys.exit()

    crefs = None
    if args.capture is not None:
        crefs = {}
        flank = args.flank
        with open(args.capture) as fp:
            h = fp.readline()
            for l in fp:
                t = l.rstrip().split('\t', 3)
                crefs.setdefault(t[0], []).append((max(int(t[1]) - flank, 0), int(t[2]) + flank - 1))
            print("Building NCLS indexes")
            for k, v in crefs.items():
                starts = NP.array([x[0] for x in v])
                ends = NP.array([x[1] for x in v])
                ids = NP.array(NP.arange(0 ,len(v)))
                crefs[k] = NCLS(starts, ends, ids)


    gsnps = {}
    srefs = {}
    sta, end, idx = NP.zeros(1, dtype='int'), NP.zeros(1, dtype='int'), NP.zeros(1, dtype='int')
    for f in args.genome:
        for r in VCF(f):
            if r.FILTER or r.is_indel or len(r.ALT) > 1:
                continue
            pos = r.POS - 1
            ref = r.REF 
            chrom = r.CHROM
            if crefs is not None:
                if chrom in crefs:
                    sta[0] = pos
                    end[0] = pos + 1
                    _, overlaps = crefs[chrom].all_overlaps_both(sta, end, idx)
                    if len(overlaps) == 0:
                        continue
                else:
                    continue
            alt = r.ALT[0]
            if (chrom, pos) not in gsnps:
                gsnps[(chrom, pos)] = (ref, alt)
                srefs.setdefault(chrom, []).append(pos)

    for k, v in srefs.items():
        srefs[k] = NP.sort(v)

    sdata = []
    print(srefs.keys())
    ckeep = set()
    for prefix, name in zip(pileups, names):
        print(f'Loading {name}: {prefix} srfs = {len(srefs)}')
        data = SampleData(f'{prefix}_annotated.h5', name, groups=groups, load_matrix=False)
        h5f = h5py.File(prefix + '_barcode_matrices.h5', 'r')
        d = h5f['coverage']
        cc, cp, ct = map(NP.array, (d['chroms'], d['pos'], d['tid']))
        cc = {c:i for i, c in enumerate(cc)}
        tot = 0
        for k, v in srefs.items():
            ci = cc[k]
            idx = NP.nonzero(ct == ci)[0]
            pp = cp[idx]
            r1, idx1, idx2 = NP.intersect1d(pp, v, assume_unique=True, return_indices=True)
            tot += len(idx1)
            for x in v[idx2]:
                ckeep.add((k, x))
        sdata.append(data)

    N = len(gsnps)
    gsnps = {k:v for k, v in gsnps.items() if k in ckeep}
    M = len(gsnps)
    del ckeep
    del srefs

    print(f'Kept {M:,} / {N:,} genome SNVs with coverage in at least one dscRNA-seq sample')

    print("Merging annotated samples")
    merged = merge_samples(sdata, groups=groups)
    if args.capture is not None:
        print("Merged annotated SNVs before capture filtering", len(merged), end=" ")
        merged = merged[merged.captured]
        print("After Captured Filtering", len(merged))

    print(f'Total merged SNVs {len(merged)}')
    all_snvs = {}
    for i, r in merged.iterrows():
        snv = SNVData(i[0], i[1], r.ref, r.alt, r.strand)
        for n in names:
            snv.extra[n] = (1 if r[f'{n}_alt_strand_count'] > 0 else 0)
        all_snvs[snv.key()] = snv

    bad = set()
    enames = []

    BEFORE = len(all_snvs)
    for k, (ref, alt) in gsnps.items():
        r = all_snvs.get(k)
        if r is not None:
            if r.ref == ref and r.alt == alt:
                r.genome = 1
            else:
                bad.add(k)
        else:
            snv = SNVData(k[0], k[1], ref, alt, '?')
            snv.genome = 1
            all_snvs[k] = snv

    print(f'SNVs from annotated samples: {BEFORE} after adding genome: {len(all_snvs)} and {len(bad)} SNVs will be filtered as they are not biallelic')
    rows = []
    for k, v in all_snvs.items():
        if k in bad:
            continue
        row = [v.chrom, v.pos, v.ref, v.alt, v.strand, v.genome]
        for n in names:
            row.append(v.extra.get(n, 0))
        for n in enames:
            row.append(v.extra.get(n, 0))
        rows.append(row)

    print(f'Writing {len(rows)} SNVs for accuracy calculations')
    df = pd.DataFrame(rows, columns=['chrom', 'pos', 'ref', 'alt', 'strand', 'genome'] + names)
    df.sort_values(['chrom', 'pos'], inplace=True)
    df.set_index(['chrom', 'pos'], inplace=True)
    df.to_csv(args.output + ".txt.gz", sep="\t")

    merged.sort_index(inplace=True)
    merged.to_csv(args.output + "_snvs.txt.gz", sep="\t")

def merge_results_cmd(args):
    import sys, argparse, os
    from .snvcmp import GroupBuilder
    parser = argparse.ArgumentParser(description='Merge SNVs from multiple VCFs and annotated SNV samples')
    parser.add_argument('-o', '--output', help='Output text file (can be gz)', type=str, required=True)
    parser.add_argument('-c', '--cols', help='Columns to add from merged data', type=str, required=False, default='strand,captured,germline,edit,g1000,transition,transversion,group')
    parser.add_argument('-a', '--aligners', help='Use these aligner prefixes to extract strand-specific count data ie. genome,scsnv,cellranger,starsolo', type=str, required=True)
    parser.add_argument('prefix', help='Result file prefix will merge data from all files', type=str)
    args = parser.parse_args(args[1:])
    

    mfile = ''
    sfile = '{args.prefix}.txt.gz'
    res = None
    aligners = list(args.aligners.split(','))
    for f in sorted(glob.glob(args.prefix + '_*.txt.gz')):
        if '_snvs' in f:
            mfile = f
            continue
        print(f)
        df = pd.read_csv(f, dtype={'chrom':object}, sep='\t', index_col=[0, 1])
        aln = None
        scell = False
        for a in aligners:
            if f'{a}_plus_ref' in df:
                aln = a
                scell = True
                break
            elif f'{a}_ref' in df:
                aln = a
                scell = False

        if aln is None:
            print('Aligner not found for df aligners = ',aligners,' columns = ', ','.join(df.columns))
            return 0


        sdata = NP.zeros((len(df), 6 if scell else 2), dtype='int')
        if scell:
            for st, sp in zip(('+', '-'), ('plus', 'minus')):
                idx = NP.where(df['strand'].values == st)[0]
                sdata[idx, 0] = df[f'{aln}_{sp}_ref'].values[idx]
                sdata[idx, 1] = df[f'{aln}_{sp}_alt'].values[idx]
                sdata[idx, 2] = df[f'{aln}_{sp}_ref_barcodes'].values[idx]
                sdata[idx, 3] = df[f'{aln}_{sp}_alt_barcodes'].values[idx]
                sdata[idx, 4] = df[f'{aln}_{sp}_altref_barcodes'].values[idx]
            sdata[:,5] = df[f'{aln}_total_barcodes'].values
            idx = NP.where(df['strand'] == '?')
            sdata[idx, 0] = df[f'{aln}_plus_ref'].values[idx] + df[f'{aln}_minus_ref'].values[idx]
            sdata[idx, 1] = df[f'{aln}_plus_alt'].values[idx] + df[f'{aln}_minus_alt'].values[idx]
            sdata[idx, 2] = df[f'{aln}_ref_barcodes'].values[idx]
            sdata[idx, 3] = df[f'{aln}_alt_barcodes'].values[idx]
            sdata[idx, 4] = df[f'{aln}_total_altref_barcodes'].values[idx]
        else:
            sdata[:, 0] = df[f'{aln}_ref'].values
            sdata[:, 1] = df[f'{aln}_alt' if f'{aln}_alt' in df else f'{aln}_alf'].values

        if res is None:
            idx = NP.where(df.columns == 'snv_idx_h5')[0][0]
            cols_to_use = df.columns[:idx + 1]
            res = df[cols_to_use].copy()
        if scell:
            names = ['ref', 'alt', 'ref_barcodes', 'alt_barcodes', 'ref_alt_barcodes', 'total_barcodes']
            names = [f'{aln}_{n}' for n in names]
        else:
            names = [f'{aln}_{n}' for n in ('ref', 'alt')]

        for j, n in enumerate(names):
            res[n] = sdata[:,j]


    groups = GroupBuilder(captured=False)
    groups.add_group('germline', 'Germline', 'germline')
    groups.add_group('g1000', 'Common', 'g1000', '~germline')
    groups.add_group('edit', 'A-to-G Edit', 'A_to_G_edit', '~germline', '~g1000')
    groups.add_group('transition', 'Transition', 'transition', '~A_to_G_edit', '~germline', '~g1000')
    groups.add_group('transversion', 'Transversion', 'transversion', '~A_to_G_edit', '~germline', '~g1000')
    print("Reading merged SNVs to fetch annotations")
    mdf = pd.read_csv(mfile, dtype={'chrom':object}, sep='\t', index_col=[0, 1])
    groups.apply(mdf)
    cols = list(args.cols.split(','))
    out = args.output

    print("Merging annotations")
    res = res.merge(mdf[cols], left_index=True,right_index=True, how='left')
    print("Writing data")
    res.fillna(False, inplace=True)
    res.to_csv(args.output, sep="\t")
