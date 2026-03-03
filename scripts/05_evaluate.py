#!/usr/bin/env python3
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2


def load_plink(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t', comment='#')
    if 'P' not in df.columns:
        raise ValueError(f'{label}: expected PLINK2 column P')
    out = pd.DataFrame({
        'method': label,
        'chr': df['#CHROM'] if '#CHROM' in df.columns else df['CHROM'],
        'pos': df['POS'],
        'snp': df['ID'],
        'p': df['P'],
    })
    return out.dropna(subset=['p'])


def load_gemma(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t')
    out = pd.DataFrame({
        'method': label,
        'chr': df['chr'],
        'pos': df['ps'],
        'snp': df['rs'],
        'p': df['p_wald'],
    })
    return out.dropna(subset=['p'])


def lambda_gc(pvals: pd.Series) -> float:
    pvals = pvals[(pvals > 0) & (pvals <= 1)]
    chisq = chi2.isf(pvals, 1)
    return float(np.median(chisq) / chi2.ppf(0.5, 1))


def qq_plot(df: pd.DataFrame, out_png: Path) -> None:
    plt.figure(figsize=(6, 6))
    for m, sub in df.groupby('method'):
        p = np.sort(sub['p'].values)
        p = p[(p > 0) & (p <= 1)]
        exp = -np.log10(np.arange(1, len(p) + 1) / (len(p) + 1))
        obs = -np.log10(p)
        plt.plot(exp, obs, '.', ms=2, label=m, alpha=0.7)
    maxv = plt.xlim()[1]
    plt.plot([0, maxv], [0, maxv], 'k--', lw=1)
    plt.xlabel('Expected -log10(p)')
    plt.ylabel('Observed -log10(p)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def manhattan_plot(df: pd.DataFrame, out_png: Path) -> None:
    methods = df['method'].unique()
    fig, axes = plt.subplots(len(methods), 1, figsize=(10, 3 * len(methods)), sharex=False)
    if len(methods) == 1:
        axes = [axes]
    for ax, m in zip(axes, methods):
        sub = df[df['method'] == m].copy()
        sub['mlog10p'] = -np.log10(np.clip(sub['p'].values, 1e-300, 1.0))
        ax.scatter(sub['pos'], sub['mlog10p'], s=4, alpha=0.7)
        ax.set_title(m)
        ax.set_xlabel('Position')
        ax.set_ylabel('-log10(p)')
    plt.tight_layout()
    plt.savefig(out_png, dpi=180)
    plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--lr', required=True)
    ap.add_argument('--lr-pcs', required=True)
    ap.add_argument('--lmm', required=True)
    ap.add_argument('--out-prefix', required=True)
    ap.add_argument('--top-n', type=int, default=100)
    args = ap.parse_args()

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    lr = load_plink(Path(args.lr), 'LR')
    lrp = load_plink(Path(args.lr_pcs), 'LR+PCs')
    lmm = load_gemma(Path(args.lmm), 'LMM')
    all_df = pd.concat([lr, lrp, lmm], ignore_index=True)

    qq_plot(all_df, out_prefix.with_suffix('.qq.png'))
    manhattan_plot(all_df, out_prefix.with_suffix('.manhattan.png'))

    lambdas = all_df.groupby('method')['p'].apply(lambda_gc).reset_index(name='lambda_gc')
    lambdas.to_csv(out_prefix.with_suffix('.lambda_gc.tsv'), sep='\t', index=False)

    top_sets = {}
    for m, sub in all_df.groupby('method'):
        top = sub.nsmallest(args.top_n, 'p')['snp']
        top_sets[m] = set(top)

    overlap_rows = []
    methods = list(top_sets.keys())
    for i in range(len(methods)):
        for j in range(i + 1, len(methods)):
            a, b = methods[i], methods[j]
            inter = len(top_sets[a] & top_sets[b])
            overlap_rows.append({'method_a': a, 'method_b': b, 'top_n': args.top_n, 'overlap': inter})
    pd.DataFrame(overlap_rows).to_csv(out_prefix.with_suffix('.top_overlap.tsv'), sep='\t', index=False)


if __name__ == '__main__':
    main()
