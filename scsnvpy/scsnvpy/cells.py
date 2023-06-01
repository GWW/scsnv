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
from .plot import rand_jitter
from .data import build_matrix, load_barcode_rates
from scipy.sparse import csr_matrix
from scipy.stats import gaussian_kde
import pandas as pd
import sys, os, argparse, gzip
import flammkuchen as fl
import matplotlib as mpl
import pylab as plt
from statsmodels.robust.scale import mad
import matplotlib.gridspec as gridspec

def dists_from_line(x, y):
    x = NP.log10(x)
    p1 = x[0], y[0]
    p2 = x[-1], y[-1]
    dists = []
    dist_p1p2 = NP.sqrt((p1[1] - p2[1])**2 + (p1[0] - p2[0]) ** 2)
    num = p2[0] * p1[1] - p1[0] * p2[1]
    y2y1 = p2[1] - p1[1]
    x2x1 = p2[0] - p1[0]
    maxd = 0
    maxi = -1
    for i in range(1, len(x) - 1):
        d = NP.abs(y2y1 * x[i] - x2x1 * y[i] + num) / dist_p1p2
        d = NP.abs(d)
        if d > maxd:
            maxd = d
            maxi = i
    return maxi

def find_start(y):
    x = NP.log10(NP.arange(1, len(y) + 1))
    #y = NP.log10(y)
    x2 = x[-1]
    y2 = y[-1]
    ups = []
    max_cnt = 0
    max_i = 0
    for i in range(1, len(x) - 1):
        slope = (y2 - y[i]) / (x2 - x[i])
        b = y2 - x2 * slope
        N = len(x) - i - 1
        cnt = ((x[i:len(x) - 1] * slope + b) <= y[i:len(x) - 1]).sum() / N
        if cnt > max_cnt:
            max_cnt = cnt
            max_i = i
        if cnt >= 1:
            break

    return max_i

def density_scatter(ax, x, y, s=10, cmap=plt.cm.inferno, logx=False, logy=False):
    if len(x) > 20000:
        ds_idx = NP.random.choice(NP.arange(0, len(x)), 20000, replace=False)
        x = x[ds_idx]
        y = y[ds_idx]
    lx, ly = x, y
    if logx:
        lx = NP.log10(x)
    if logy:
        ly = NP.log10(y)
    xy = NP.vstack([lx,ly])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]

    ret = ax.scatter(x, y, c=z, s=s, lw=0, cmap=cmap, marker='o')
    try:
        if logx:
            ax.set_xscale('log', subs=[2,3,4,5,6,7,8,9])
        if logy:
            ax.set_yscale('log', subs=[2,3,4,5,6,7,8,9])
    except:
        if logx:
            ax.set_xscale('log', subsx=[2,3,4,5,6,7,8,9])
        if logy:
            ax.set_yscale('log', subsy=[2,3,4,5,6,7,8,9])


def MT_cutoff(ax, df, key, mads, max_value, mt_perc):
    if key not in df.columns:
        print('Error {key} not in the barcode data columns')
        sys.exit()
    mols = df['molecules'].values
    pmt = 100.0 * df[key].values / mols
    med = NP.median(pmt)
    if mt_perc < 0.0:
        md = mad(pmt)
        co = min(mads * md + med, max_value)
    else:
        co = mt_perc

    ax.axhline(med, color=plt.cm.Set1.colors[0])
    ax.axhline(co, color=plt.cm.Set1.colors[0], ls='--')
    density_scatter(ax, mols, pmt, logx=True)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.grid(color='0.5', lw=0.5, ls='--')
    ax.set_xlabel('Molecules')
    ax.set_ylabel('MT (%)')
    passed = pmt <= co
    df['passed'] = df['passed'] & passed
    ax.set_title(f'Cutoff = {co:.2f}%\nPassed = {passed.sum():,} / {len(df):,}')

def dups_cutoff(ax, df, mads, max_mads):
    idx = df['passed'].copy()
    mols = df['molecules'].values[idx]
    pdup = 100.0 * df['pcr_dups'].values[idx] / df['reads'].values[idx]
    md = mad(pdup)
    med = NP.median(pdup)
    co = min(mads * md, max_mads)

    ax.axhline(med, color=plt.cm.Set1.colors[0])
    ax.axhline(med - co, color=plt.cm.Set1.colors[0], ls='--')
    density_scatter(ax, mols, pdup, logx=True)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.grid(color='0.5', lw=0.5, ls='--')
    ax.set_xlabel('Molecules')
    ax.set_ylabel('PCR Duplicates (%)')
    passed = pdup >= (med - co)
    df.loc[idx, 'passed'] &= passed
    ax.set_title(f'Cutoff = {med - co:.2f}%\nPassed = {passed.sum():,} / {idx.sum():,}')

def sat_cutoff(ax, df, mads, max_mads):
    idx = df['passed'].copy()
    mols = df['molecules'].values[idx]
    psat = 100 - 100.0 * mols / df['reads'].values[idx]
    md = mad(psat)
    med = NP.median(psat)
    co = min(mads * md, max_mads)

    ax.axhline(med, color=plt.cm.Set1.colors[0])
    ax.axhline(med - co, color=plt.cm.Set1.colors[0], ls='--')
    density_scatter(ax, mols, psat, logx=True)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.grid(color='0.5', lw=0.5, ls='--')
    ax.set_xlabel('Molecules')
    ax.set_ylabel('Saturation (%)')
    passed = psat >= (med - co)
    df.loc[idx, 'passed'] &= passed
    ax.set_title(f'Cutoff = {med - co:.2f}%\nPassed = {passed.sum():,} / {idx.sum():,}')

def barcode_cutoff(ax, df, mads, max_mads):
    idx = df['passed'].copy()
    mols = df['molecules'].values[idx]
    bcor = 100.0 * df['align_barcode_corrected'].values[idx] / df['align_total_reads'].values[idx]
    md = mad(bcor)
    med = NP.median(bcor)
    co = min(mads * md, max_mads)

    ax.axhline(med, color=plt.cm.Set1.colors[0])
    ax.axhline(med + co, color=plt.cm.Set1.colors[0], ls='--')
    density_scatter(ax, mols, bcor, logx=True)
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.grid(color='0.5', lw=0.5, ls='--')
    ax.set_xlabel('Molecules')
    ax.set_ylabel('Barcodes Corrected (%)')
    passed = bcor <= (med + co)
    df.loc[idx, 'passed'] &= passed
    ax.set_title(f'Cutoff = {med + co:.2f}%\nPassed = {passed.sum():,} / {idx.sum():,}')


def cells_cmd(cargs):
    parser = argparse.ArgumentParser(description='Find barcodes that represent cells')
    parser.add_argument('-u', '--min-umis', help='Minimum number of spliced UMIs for a cell', type=int, default=750)
    parser.add_argument('--skip-mt', help='Skip MT percent filtering', action='store_true')
    parser.add_argument('--mt-perc', help='Set the MT cutoff to this percent, set to -1 to enable the automatic cutoff below', type=float, default=25)
    parser.add_argument('--mt-mad', help='MT cutoff of median + X * (median absolute deviation)', type=float, default=5)
    parser.add_argument('--mt-max', help='Maximum MT Percent', type=float, default=25)
    parser.add_argument('--mt-key', help='Gene group name for MT reads', type=str, default='group_MT')
    parser.add_argument('--pcr-mad', help='PCR duplicate cutoff of median - X * (median absolute deviation)', type=float, default=5)
    parser.add_argument('--pcr-mad-max', help='Maximum PCR cutoff difference ie a cap for X * (median absolute deviation)', type=float, default=10)
    parser.add_argument('--saturation-mad', help='Saturation duplicate cutoff of median - X * (median absolute deviation)', type=float, default=3)
    parser.add_argument('--saturation-mad-max', help='Maximum Saturation cutoff difference ie a cap for X * (median absolute deviation)', type=float, default=10)
    parser.add_argument('--correct-mad', help='Barcode correction cutoff of median + X * (median absolute deviation)', type=float, default=10)
    parser.add_argument('--correct-mad-max', help='Barcode correction cutoff difference ie a cap for X * (median absolute deviation)', type=float, default=1.5)
    parser.add_argument('-o', '--output', help='Output prefix', type=str, required=True)
    parser.add_argument('-n', '--name', help='Sample Name instead of using summary file', type=str)
    parser.add_argument('summary', help='The summary.h5 file')
    args = parser.parse_args(cargs[1:])

    sample = args.summary
    brates = load_barcode_rates(sample)

    krates = brates[brates.molecules >= 1].copy()

    mols = krates['molecules'].values
    cc = NP.cumsum(mols)
    cc = cc / cc[-1]

    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['xtick.direction'] = 'out'
    mpl.rcParams['ytick.direction'] = 'out'
    mpl.rcParams['xtick.top'] = False
    mpl.rcParams['ytick.right'] = False

    fig = plt.figure(constrained_layout=True, figsize=(12, 12), dpi=250)
    gs = fig.add_gridspec(3, 3)

    gaxs = [gs[0, 0], gs[0, 1], gs[0, 2], gs[1, 0], gs[1, 1], gs[1, 2], gs[2,:]]
    axs = [fig.add_subplot(g) for g in gaxs]
    #fig.subplots_adjust(hspace=0.45, wspace=0.45)
    #axs = axs.flatten()

    colors = plt.cm.Set1.colors

    ax = axs[0]
    bidx = NP.arange(1, len(mols) + 1)
    ax.plot(bidx, 100.0 * cc, lw=0.5, color='k')
    ax.scatter(bidx, 100.0 * cc, s=5, color=colors[1], lw=0)
    start = find_start(cc)
    idx = dists_from_line(bidx[start:], cc[start:]) + start + 1
    ax.axvline(idx, lw=1, color=colors[0])
    ax.set_ylabel('Cumulative molecules (%)')
    ax.set_xlabel('Cumulative barcodes')
    sname = args.name if 'name' in args else sample
    ax.set_title(f'{sname}\nBarcodes Passed = {idx:,}')
    ax.set_axisbelow(True)
    ax.grid(color='0.5', lw=0.5, ls='--')
    try:
        ax.set_xscale('log', subs=[2,3,4,5,6,7,8,9])
    except:
        ax.set_xscale('log', subsx=[2,3,4,5,6,7,8,9])
    ax.minorticks_on()


    prates = krates.iloc[:idx].copy()
    prates['passed'] = True
    if not args.skip_mt:
        MT_cutoff(axs[1], prates, args.mt_key, args.mt_mad, args.mt_max, args.mt_perc)
    else:
        print('Skipping MT DNA filtering')
    dups_cutoff(axs[2], prates, args.pcr_mad, args.pcr_mad_max)
    sat_cutoff(axs[3], prates, args.saturation_mad, args.saturation_mad_max)
    barcode_cutoff(axs[4], prates, args.correct_mad, args.correct_mad_max)

    ax = axs[5]
    min_umi = args.min_umis
    passed = prates[prates.passed]['molecules'] >= min_umi
    prates['passed'] &= prates['molecules'] >= min_umi
    krates['passed'] = False
    krates['used'] = False
    krates.iloc[prates.index, krates.columns.get_loc('used')] = True
    krates.iloc[prates.index, krates.columns.get_loc('passed')] = prates['passed']
    ax.plot(bidx, mols, lw=0.5, color='k')
    ax.scatter(bidx[~krates['passed'].values], krates[~krates['passed']].molecules, s=10, color=colors[0], label='Low Quality', alpha=0.5, lw=0)
    idx = krates['used'] & krates['passed']
    ax.set_title(f'Min UMI Cutoff\nTotal passed = {passed.sum():,} / {len(passed)}')
    ax.scatter(bidx[idx], krates[idx].molecules, s=10, color=colors[2], label='Passed', alpha=0.5, lw=0)
    ax.set_ylabel('Total Molecules', fontsize=10)
    ax.set_xlabel('Cumulative barcodes', fontsize=10)
    ax.axhline(max(min_umi, krates[~krates['passed']].molecules.min()), lw=1, color=colors[2], ls='--')
    ax.set_axisbelow(True)
    ax.grid(color='0.5', lw=0.5, ls='--')
    try:
        ax.set_yscale('log', subs=[2,3,4,5,6,7,8,9])
        ax.set_xscale('log', subs=[2,3,4,5,6,7,8,9])
    except:
        ax.set_yscale('log', subsy=[2,3,4,5,6,7,8,9])
        ax.set_xscale('log', subsx=[2,3,4,5,6,7,8,9])
    ax.minorticks_on()






    ax = axs[6]
    box_keys = ['align_cdna', 'align_intronic', 'align_intergenic', 'align_multimapped', 'align_antisense', 'align_unmapped',
            'align_ambiguous', 'align_barcode_corrected', 'align_umi_qa_fail', 'align_tag_qa_fail']

    passed = krates[idx].copy()
    passed.reset_index(inplace=True, drop=True)

    for i, key in enumerate(box_keys):
        pvals = 100.0 * passed[key] / passed['align_total_reads']
        ax.scatter(rand_jitter(len(pvals), (i - 0.4, i + 0.4)), pvals, s=3, color=colors[1], lw=0, marker='o')
        ax.boxplot(pvals, positions=[i], showfliers=False, manage_ticks=False, widths=0.9, medianprops = dict(linestyle='-', linewidth=1.5, color='black'))

    ax.set_title('Passed Barcode Mapping Rates')
    ax.set_xticks(NP.arange(0, len(box_keys)))
    ax.set_xticklabels(['Spliced', 'Unspliced', 'Intergenic', 'Multimapped', 'Antisense', 'Unmapped', 'Ambiguous', 'Barcode Corrected', 'UMI Fail', 'Tag Fail'], rotation=45, ha='right')
    ax.set_ylabel('Reads (%)')


    fig.savefig(os.path.join(args.output, 'cells.png'), bbox_inches='tight')
    plt.close(fig)

    gene_ids, gene_names, mat, imat = build_matrix(sample, passed['barcode_id'].values, matrix_type='csr')

    from collections import Counter
    unq = Counter(gene_names)
    if len(unq) < len(gene_names):
        cc = Counter()
        for i, g in enumerate(gene_names):
            if unq[g] > 1:
                cc[g] += 1
                gene_names[i] = f'{g}-{cc[g]}'

    import anndata
    obs = {}
    for k in passed.columns:
        if k == 'barcode':
            obs[k] = NP.array(passed[k], dtype='str')
        elif k == 'index':
            continue
        else:
            obs[k] = NP.array(passed[k], dtype='int32')
    print(f'Total Passed Cells {len(passed)}')
    layers = {'unspliced':imat.T, 'spliced':mat.copy().T}
    adata = anndata.AnnData(mat.T, obs, {'var_names':gene_names, 'gene_ids':gene_ids}, layers=layers)
    adata.write(os.path.join(args.output, 'anndata.h5ad'), compression="lzf")

    genes = pd.DataFrame({'gene_id':gene_ids, 'gene_name':gene_names})
    print(mat.shape, imat.shape)
    out = {'genes':genes, 'cells':passed, 'spliced':mat, 'unspliced':imat}
    fl.save(os.path.join(args.output, 'cells.h5'), out)

    with gzip.open(os.path.join(args.output, 'passed_barcodes.txt.gz'), 'wt') as fout:
        fout.write("barcode\n")
        for b in passed['barcode'].values:
            fout.write(b + "\n")

