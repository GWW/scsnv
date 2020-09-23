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
import matplotlib
import pylab as plt

def format_axis(ax, hide_xticks = False, hide_yticks = False):
    ax.get_yaxis().tick_left()
    ax.get_xaxis().tick_bottom()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.get_yaxis().set_tick_params(direction='out')
    ax.get_xaxis().set_tick_params(direction='out')
    if hide_xticks:
        for t in ax.xaxis.get_ticklines(): 
            t.set_visible(False)
    if hide_yticks:
        for t in ax.yaxis.get_ticklines(): 
            t.set_visible(False)
    return ax

def format_subplots(r, c, bsize = None, dpi=150, **kwargs):
    if bsize is not None:
        kwargs['figsize'] = (c * bsize, r * bsize)
    fig, axs = plt.subplots(r, c, dpi=dpi, **kwargs)
    if kwargs.get('squeeze', True) and r == 1 and c == 1:
        format_axis(axs)
    else:
        for ax in axs.flat:
            format_axis(ax)
    return fig, axs

def add_ax(fig, lft, top, width, height, **kwargs):
    ax = fig.add_axes([lft, 1 - top, width, height], **kwargs)
    return format_axis(ax)

fsp = format_subplots
fax = format_axis


def rand_jitter(N, range, stdev = 0.120):
    mid = (range[1] - range[0]) / 2 + range[0]
    ret = NP.zeros(N, dtype='float') + NP.random.randn(N) * stdev + mid
    ret[ret > range[1]] = range[1]
    ret[ret < range[0]] = range[0]
    return ret

def get_colors():
    cb = plt.cm.tab20(NP.linspace(0, 1, 20))
    c1 = cb[NP.arange(0, 20, 2)]
    c2 = cb[NP.arange(1, 20, 2)]
    c3 = plt.cm.Dark2(NP.linspace(0, 1, 8))
    vcolors = NP.concatenate((c3, c1, c2))
    return vcolors

def density_scatter(ax, x, y, s=10, cmap=plt.cm.plasma):
    xy = NP.vstack([x,y])
    z = gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    return ax.scatter(x, y, c=z, s=s, edgecolor='', cmap=cmap, marker='o')
