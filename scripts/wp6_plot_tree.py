import os
#!/usr/bin/env python3
"""Render the 46-isolate combined tree (19 study + 27 Saudi context) as an
annotated publication figure: study isolates coloured by genotype, context
isolates in grey, major clades labelled."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re

TREE = os.path.expanduser('~/Projects/cauris-wgs-analysis/analysis/07_context/tree_combined.contree')
OUT_PDF = os.path.expanduser('~/Projects/cauris-wgs-analysis/output/Figure_WP6_combined_tree.pdf')
OUT_PNG = os.path.expanduser('~/Projects/cauris-wgs-analysis/output/Figure_WP6_combined_tree.png')

K143B = {'Sample01','Sample04','Sample05','Sample14','Sample15','Sample16','Sample19','Sample20'}
K143A = {'Sample06','Sample07','Sample09','Sample10','Sample13'}
Y132 = {'Sample03','Sample08','Sample11','Sample12','Sample17','Sample18'}
STUDY = K143B | K143A | Y132


def parse(s):
    pos = [0]
    def node():
        if s[pos[0]] == '(':
            pos[0] += 1
            children = []
            while True:
                children.append(node())
                if s[pos[0]] == ',':
                    pos[0] += 1
                    continue
                pos[0] += 1  # ')'
                break
            label = ''
            dist = 0.0
            if pos[0] < len(s) and s[pos[0]] != ')' and s[pos[0]] != ',':
                if s[pos[0]] == ':':
                    pos[0] += 1
                    st = pos[0]
                    while pos[0] < len(s) and s[pos[0]] not in ',)':
                        pos[0] += 1
                    dist = float(s[st:pos[0]])
                else:
                    st = pos[0]
                    while pos[0] < len(s) and s[pos[0]] not in ':,()':
                        pos[0] += 1
                    label = s[st:pos[0]]
                    if pos[0] < len(s) and s[pos[0]] == ':':
                        pos[0] += 1
                        st = pos[0]
                        while pos[0] < len(s) and s[pos[0]] not in ',)':
                            pos[0] += 1
                        dist = float(s[st:pos[0]])
            return ('', children, label, dist)
        st = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in ':,()':
            pos[0] += 1
        name = s[st:pos[0]]
        dist = 0.0
        if pos[0] < len(s) and s[pos[0]] == ':':
            pos[0] += 1
            st = pos[0]
            while pos[0] < len(s) and s[pos[0]] not in ',)':
                pos[0] += 1
            dist = float(s[st:pos[0]])
        return (name, None, '', dist)
    return node()


def tips(n, acc):
    name, ch, label, dist = n
    if ch is None:
        acc.append((name, dist))
    else:
        for c in ch:
            tips(c, acc)


root = parse(open(TREE).read())
tv = []
tips(root, tv)
n = len(tv)
y = {name: n - 1 - i for i, (name, d) in enumerate(tv)}

color = {}
for k in K143B: color[k] = '#1a237e'
for k in K143A: color[k] = '#00897b'
for k in Y132: color[k] = '#d62728'
for name, d in tv:
    if name not in color:
        color[name] = '#9e9e9e' if name.startswith(('ERR', 'SRR')) else '#000000'

fig, ax = plt.subplots(figsize=(11, 17))
ax.set_ylim(-1, n + 1)
ax.set_xlim(-0.005, 0.34)
ax.axis('off')


def draw(nod, px):
    name, ch, label, dist = nod
    if ch is None:
        py = y[name]
        x = px + dist
        ax.plot([px, x], [py, py], 'k-', lw=1.0)
        ax.scatter([x], [py], color=color.get(name, 'gray'), s=8, zorder=5)
        fs = 7.5 if name.startswith(('ERR', 'SRR')) else 9
        fw = 'normal' if name.startswith(('ERR', 'SRR')) else 'bold'
        ax.text(x + 0.001, py, name, va='center', fontsize=fs, color=color.get(name, 'gray'), fontweight=fw)
        return x, py
    x = px + dist
    pts = [draw(c, x) for c in ch]
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    ymin, ymax = min(ys), max(ys)
    ax.plot([x, x], [ymin, ymax], 'k-', lw=1.0)
    if label:
        ax.text(x + 0.001, (ymin + ymax) / 2, label, fontsize=7, color='0.4', fontweight='bold')
    return x, (ymin + ymax) / 2


for c in root[1]:
    draw(c, 0)

from matplotlib.lines import Line2D
leg = [Line2D([0], [0], marker='s', color='w', markerfacecolor='#1f77b4', markersize=8, label='Study: K143R/CDR1_V704L (clade B)'),
       Line2D([0], [0], marker='s', color='w', markerfacecolor='#00897b', markersize=8, label='Study: K143R/CDR1_V704L (clade A)'),
       Line2D([0], [0], marker='s', color='w', markerfacecolor='#d62728', markersize=8, label='Study: Y132F'),
       Line2D([0], [0], marker='s', color='w', markerfacecolor='#9e9e9e', markersize=8, label='Saudi context (Guan et al. 2025; Chow et al. 2020)')]
ax.legend(handles=leg, loc='upper left', fontsize=8, frameon=False)
ax.set_title('Combined phylogeny: 19 study isolates + 27 Saudi context isolates (IQ-TREE2, 1,301 core SNPs)', fontsize=9)
plt.tight_layout()
plt.savefig(OUT_PDF)
plt.savefig(OUT_PNG, dpi=200)
print('combined tree figure:', OUT_PDF)
