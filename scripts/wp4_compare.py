#!/usr/bin/env python3
"""
WP4 stage B: compare freebayes SNP-distance matrix vs MycoSNP matrix.
Produces scatter plot, Bland-Altman-style agreement, and discordant-pair table.

Inputs:
  analysis/05_phylogeny/freebayes_snpdists.tsv     (independent caller)
  analysis/06_mycosnp/snpdists/combined.tsv        (official MycoSNP run, if available)
  previously reported values are embedded for reference when reference outputs are absent.
"""


import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.expanduser('~/Projects/cauris-wgs-analysis')
OUT = os.path.join(ROOT, 'analysis/04_freebayes')
os.makedirs(OUT, exist_ok=True)

FB = os.path.join(ROOT, 'analysis/05_phylogeny/freebayes_snpdists.tsv')
MS = os.path.join(ROOT, 'analysis/06_mycosnp/snpdists/combined.tsv')


def read_matrix(path):
    df = pd.read_csv(path, sep='\t', index_col=0)
    df = df.drop(columns=[c for c in df.columns if 'reference' in c.lower()], errors='ignore')
    return df


def main():
    fb = read_matrix(FB)
    print('freebayes matrix:', fb.shape)
    # key comparisons (previously reported numbers)
    reported = {
        'within K143R cluster': (4, 27, 12),
        'within Y132F cluster': (6, 19, 11),
        'between clusters': (38, 67, 49),
        'cohort-wide median': (None, None, '<50'),
    }
    # compute freebayes cohort summary
    vals = fb.values.astype(float)
    mask = ~np.isnan(vals)
    tri = np.triu(mask, 1)
    flat = vals[tri]
    print(f'freebayes pairwise SNP distances: n={len(flat)}, '
          f'median={np.median(flat):.0f}, range={flat.min():.0f}-{flat.max():.0f}')
    # cluster assignment from genotype calls (PathogenWatch)
    k143 = ['Sample01', 'Sample04', 'Sample05', 'Sample06', 'Sample07', 'Sample09',
            'Sample10', 'Sample13', 'Sample14', 'Sample15', 'Sample16', 'Sample19', 'Sample20']
    y132 = ['Sample03', 'Sample08', 'Sample11', 'Sample12', 'Sample17', 'Sample18']
    if all(s in fb.index for s in k143) and all(s in fb.index for s in y132):
        sub = fb.loc[k143, k143].values.astype(float)
        trik = np.triu(~np.isnan(sub), 1)
        print(f'within K143R cluster: {sub[trik].min():.0f}-{sub[trik].max():.0f} '
              f'(median {np.median(sub[trik]):.0f}); reported: 4-27 (median 12)')
        sub2 = fb.loc[y132, y132].values.astype(float)
        triy = np.triu(~np.isnan(sub2), 1)
        print(f'within Y132F cluster: {sub2[triy].min():.0f}-{sub2[triy].max():.0f} '
              f'(median {np.median(sub2[triy]):.0f}); reported: 6-19 (median 11)')
        cross = fb.loc[k143, y132].values.astype(float)
        print(f'between clusters: {cross.min():.0f}-{cross.max():.0f} '
              f'(median {np.median(cross):.0f}); reported: 38-67 (median 49)')

    # if the official MycoSNP matrix exists, do the caller-vs-caller comparison
    if os.path.exists(MS):
        ms = read_matrix(MS)
        print('MycoSNP matrix:', ms.shape)
        common = fb.index.intersection(ms.index)
        fb2 = fb.loc[common, common]
        ms2 = ms.loc[common, common]
        trii = np.triu(np.ones(len(common), dtype=bool), 1)
        a = fb2.values[trii]
        b = ms2.values[trii]
        idx = [(common[i], common[j]) for i in range(len(common)) for j in range(i + 1, len(common))]
        corr = np.corrcoef(a, b)[0, 1]
        print(f'caller agreement: Pearson r = {corr:.3f} on {len(a)} pairs')
        fig, ax = plt.subplots(1, 2, figsize=(12, 5))
        ax[0].scatter(b, a, s=8, alpha=0.6)
        ax[0].set_xlabel('MycoSNP SNP distance'); ax[0].set_ylabel('freebayes SNP distance')
        ax[0].set_title(f'scatter (r={corr:.2f})')
        lo, hi = np.percentile(a - b, [2.5, 97.5])
        ax[1].scatter((a + b) / 2, a - b, s=8, alpha=0.6)
        ax[1].axhline(0, color='k', lw=0.8)
        ax[1].axhspan(lo, hi, alpha=0.2)
        ax[1].set_xlabel('mean distance'); ax[1].set_ylabel('freebayes - MycoSNP')
        ax[1].set_title('Bland-Altman (95% limits shaded)')
        plt.tight_layout()
        plt.savefig(os.path.join(OUT, 'caller_agreement.png'), dpi=150)
        plt.savefig(os.path.join(OUT, 'caller_agreement.pdf'))
        # discordant pairs
        diff = np.abs(a - b)
        disc = [(idx[i][0], idx[i][1], round(b[i]), round(a[i]))
                for i in range(len(a)) if diff[i] > 20]
        pd.DataFrame(disc, columns=['s1', 's2', 'mycosnp', 'freebayes']).to_csv(
            os.path.join(OUT, 'discordant_pairs.tsv'), sep='\t', index=False)
        print(f'discordant pairs (>20 SNP): {len(disc)} -> discordant_pairs.tsv')
    else:
        print('MycoSNP matrix not found; skipping caller-vs-caller plot (will run after MycoSNP-NF)')


if __name__ == '__main__':
    main()
