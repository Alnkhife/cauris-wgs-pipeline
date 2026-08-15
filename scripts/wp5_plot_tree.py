import os
#!/usr/bin/env python3
"""Render the IQ-TREE consensus tree as a publication-quality annotated vector figure.
Tips coloured by resistance genotype; three lineages labelled.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO

import scipy.cluster.hierarchy as sch

TREE = os.path.expanduser('~/Projects/cauris-wgs-analysis/analysis/05_phylogeny/tree_freebayes.contree')
OUT_PDF = os.path.expanduser('~/Projects/cauris-wgs-analysis/output/Figure_4_phylogeny_freebayes.pdf')
OUT_PNG = os.path.expanduser('~/Projects/cauris-wgs-analysis/output/Figure_4_phylogeny_freebayes.png')

K143R = {'Sample01','Sample04','Sample05','Sample06','Sample07','Sample09','Sample10',
         'Sample13','Sample14','Sample15','Sample16','Sample19','Sample20'}
Y132F = {'Sample03','Sample08','Sample11','Sample12','Sample17','Sample18'}

# lightweight Newick reader -> dict of parent->children with branch lengths
def parse_newick(s):
    s = s.strip()
    # recursive parser
    pos = [0]
    def parse():
        # returns (name, children(list of (name, dist)), dist)
        if s[pos[0]] == '(':
            pos[0] += 1
            children = []
            while True:
                child = parse()
                children.append(child)
                if s[pos[0]] == ',':
                    pos[0] += 1
                    continue
                elif s[pos[0]] == ')':
                    pos[0] += 1
                    break
            name = ''
            dist = 0.0
            if pos[0] < len(s) and (s[pos[0]] == ':' or s[pos[0]].isalpha()):
                # internal node label (bootstrap)
                if s[pos[0]] != ':':
                    start = pos[0]
                    while pos[0] < len(s) and s[pos[0]] not in ':,()':
                        pos[0] += 1
                    name = s[start:pos[0]]
                if pos[0] < len(s) and s[pos[0]] == ':':
                    pos[0] += 1
                    start = pos[0]
                    while pos[0] < len(s) and s[pos[0]] not in ',()':
                        pos[0] += 1
                    dist = float(s[start:pos[0]])
            return ('', children, name, dist)
        else:
            start = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ':,()':
                pos[0] += 1
            name = s[start:pos[0]]
            dist = 0.0
            if pos[0] < len(s) and s[pos[0]] == ':':
                pos[0] += 1
                start = pos[0]
                while pos[0] < len(s) and s[pos[0]] not in ',()':
                    pos[0] += 1
                dist = float(s[start:pos[0]])
            return (name, None, '', dist)
    root = parse()
    return root

def collect(root, tips, internal, parent_dist):
    name, children, label, dist = root
    if children is None:
        tips.append((name, dist))
    else:
        internal.append((name, label, dist, parent_dist))
        for c in children:
            collect(c, tips, internal, dist)

root = parse_newick(open(TREE).read())
tips, internal = [], []
collect(root, tips, internal, 0.0)

# assign x positions by tree order (left-to-right tips)
y = {}
for i, (name, d) in enumerate(tips):
    y[name] = len(tips) - 1 - i

fig, ax = plt.subplots(figsize=(7, 8))
colors = {'Sample01': '#1a237e', 'Sample04': '#1a237e', 'Sample05': '#1a237e',
          'Sample14': '#1a237e', 'Sample15': '#1a237e', 'Sample16': '#1a237e',
          'Sample19': '#1a237e', 'Sample20': '#1a237e',
          'Sample06': '#00897b', 'Sample07': '#00897b', 'Sample09': '#00897b',
          'Sample10': '#00897b', 'Sample13': '#00897b',
          'Sample03': '#d62728', 'Sample08': '#d62728', 'Sample11': '#d62728',
          'Sample12': '#d62728', 'Sample17': '#d62728', 'Sample18': '#d62728'}

def draw(node, parent_x, parent_y):
    name, children, label, dist = node
    if children is None:
        x = parent_x + dist
        ax.plot([parent_x, x], [parent_y, parent_y], 'k-', lw=1)
        ax.plot([x, x], [parent_y, parent_y], 'k-', lw=1)
        col = colors.get(name, 'gray')
        ax.scatter([x], [parent_y], color=col, s=40, zorder=5)
        ax.text(x + 0.002, parent_y, name, va='center', ha='left', fontsize=8,
                color=col, fontweight='bold')
        return x
    x = parent_x + dist
    child_xs = []
    for c in children:
        cx = draw(c, x, y[c[0]] if c[1] is None else None)
        child_xs.append(cx)
    # vertical line
    ys = [y[tip[0]] for tip in tips if tip[0] in [t[0] for t in tips]]
    # simpler: draw vertical connector between min/max child y
    child_ys = []
    for c in children:
        if c[1] is None:
            child_ys.append(y[c[0]])
        else:
            # find descendant tips
            sub = []
            def sub_tips(n):
                if n[1] is None:
                    sub.append(n[0])
                else:
                    for cc in n[1]:
                        sub_tips(cc)
            sub_tips(c)
            child_ys.extend(y[t] for t in sub)
    ymin, ymax = min(child_ys), max(child_ys)
    ax.plot([x, x], [ymin, ymax], 'k-', lw=1)
    if label:
        ax.text(x + 0.001, (ymin + ymax) / 2, label, fontsize=7, color='0.4', va='center')
    return x

ax.set_xlim(0, 0.35)
ax.set_ylim(-1, len(tips) + 1)
ax.axis('off')
# root node
name, children, label, dist = root
for c in children:
    draw(c, 0, y[c[0]] if c[1] is None else None)

from matplotlib.lines import Line2D
legend = [Line2D([0], [0], marker='s', color='w', markerfacecolor='#1a237e', markersize=10, label='K143R/CDR1_V704L (clade B - navy)'),
          Line2D([0], [0], marker='s', color='w', markerfacecolor='#00897b', markersize=10, label='K143R/CDR1_V704L (clade A - teal)'),
          Line2D([0], [0], marker='s', color='w', markerfacecolor='#d62728', markersize=10, label='Y132F')]
ax.legend(handles=legend, loc='upper left', fontsize=8, frameon=False)
ax.set_title('Maximum-likelihood phylogeny of 19 C. auris isolates (IQ-TREE2, freebayes SNP set)', fontsize=9)
plt.tight_layout()
plt.savefig(OUT_PDF)
plt.savefig(OUT_PNG, dpi=200)
print('tree figure written:', OUT_PDF)
