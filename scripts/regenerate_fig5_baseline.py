#!/usr/bin/env python3
"""Regenerate the n=14 random position/Fourier baseline and arithmetic points."""
from __future__ import annotations
import argparse, csv, sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from subset_states.core import entropy_from_state_vector, entropy_from_support, qft_state_from_support, summary_stats
from subset_states.experiments import almost_prime_union_supports, m_grid, random_subset_qft_samples

def write(path, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0].keys()));w.writeheader();w.writerows(rows)

def main():
    p=argparse.ArgumentParser();p.add_argument('--n',type=int,default=14);p.add_argument('--points',type=int,default=100);p.add_argument('--samples',type=int,default=100);p.add_argument('--seed',type=int,default=20250604);p.add_argument('--outdir',type=Path,default=ROOT/'data');a=p.parse_args()
    rows=[]
    for off,m in enumerate(m_grid(1,1<<a.n,a.points)):
        pos,four=random_subset_qft_samples(a.n,int(m),a.samples,seed=a.seed+off)
        ps,fs=summary_stats(pos),summary_stats(four)
        rows.append({'m':int(m),'sample_count':a.samples,'position_mean':ps.mean,'position_std':ps.std,'position_sem':ps.sem,'qft_mean':fs.mean,'qft_std':fs.std,'qft_sem':fs.sem,'delta_qft_minus_position':fs.mean-ps.mean})
        if off%10==0: print(f'baseline {off+1}/{a.points}')
    unions=[]
    for k,s in almost_prime_union_supports(a.n):
        unions.append({'k':k,'m':int(s.size),'position_entropy':entropy_from_support(a.n,s),'qft_entropy':entropy_from_state_vector(a.n,qft_state_from_support(a.n,s))})
    write(a.outdir/'fig5_random_qft_summary.csv',rows);write(a.outdir/'fig5_almost_prime_unions.csv',unions)
if __name__=='__main__': main()
