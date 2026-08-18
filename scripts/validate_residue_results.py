#!/usr/bin/env python3
from __future__ import annotations
import csv, json, math, sys
from itertools import combinations
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from subset_states.core import entropy_from_support, matrix_from_support
from subset_states.experiments import almost_prime_union_supports

def H(counts):
 c=np.asarray(counts,dtype=float);p=c[c>0]/c.sum();return float(-np.sum(p*np.log2(p))) if p.size else 0.0

def direct_entropy(n,support):
 omega=matrix_from_support(n,support);rho=omega.conj().T@omega
 w=np.linalg.eigvalsh(rho).real;w=w[w>1e-14];w=w/w.sum();return float(-np.sum(w*np.log2(w)))

def main():
 out=[];n=4;N=1<<n;m=n//2
 max_violation={1:-1e9,2:-1e9};max_direct=0.0;checked=0;single_class_violation=-1e9
 for mask in range(1,1<<N):
  s=np.fromiter((i for i in range(N) if mask>>i&1),dtype=np.int64)
  S=entropy_from_support(n,s);max_direct=max(max_direct,abs(S-direct_entropy(n,s)))
  for t in (1,2):
   q=1<<t;c=np.bincount(s%q,minlength=q);bound=m-t+H(c);max_violation[t]=max(max_violation[t],S-bound)
   if np.count_nonzero(c)==1:single_class_violation=max(single_class_violation,S-(m-t))
  checked+=1
 out.append(f'exhaustive_nonempty_supports_n4={checked}')
 out.append(f'max_direct_entropy_disagreement={max_direct:.3e}')
 for t in (1,2):out.append(f'max_residue_bound_violation_t{t}={max_violation[t]:.3e}')
 out.append(f'max_single_residue_factorization_bound_violation={single_class_violation:.3e}')
 # Prime-state parity consequence.
 prime=dict(almost_prime_union_supports(14))[1];M=len(prime);c=np.bincount(prime%2,minlength=2);bound=7-1+H(c);S=entropy_from_support(14,prime)
 out.append(f'prime_n14_M={M}')
 out.append(f'prime_parity_counts={c.tolist()}')
 out.append(f'prime_residue_entropy_bits={H(c):.12f}')
 out.append(f'prime_parity_ceiling_bits={bound:.12f}')
 out.append(f'prime_entropy_bits={S:.12f}')
 out.append(f'prime_forced_gap_from_balanced_max_bits={7-bound:.12f}')
 # Summary/raw consistency and exact count matching.
 summary=list(csv.DictReader(open(ROOT/'data/residue_matched_summary.csv')))
 raw=list(csv.DictReader(open(ROOT/'data/residue_matched_samples.csv')))
 max_mean_err=0.;max_std_err=0.;all_counts=True;all_below=True
 for r in summary:
  k=int(r['k']);t=int(r['matched_low_bits']);vals=[x for x in raw if int(x['k'])==k and int(x['matched_low_bits'])==t]
  pos=np.array([float(x['position_entropy']) for x in vals]);four=np.array([float(x['fourier_entropy']) for x in vals])
  max_mean_err=max(max_mean_err,abs(pos.mean()-float(r['null_position_mean'])),abs(four.mean()-float(r['null_fourier_mean'])))
  max_std_err=max(max_std_err,abs(pos.std(ddof=1)-float(r['null_position_std'])),abs(four.std(ddof=1)-float(r['null_fourier_std'])))
  all_below &= float(r['structured_position_entropy'])<pos.min() and float(r['structured_fourier_entropy'])<four.min()
  # regenerate residue counts from deterministic seed for a subset of 10 samples
  q=1<<t;target=np.array(json.loads(r['residue_counts']),dtype=int)
  # Raw file does not store supports, so verify sampler separately.
  rng=np.random.default_rng(int(r['seed']));N=1<<int(r['n'])
  for _ in range(10):
   parts=[]
   for residue,count in enumerate(target):
    if count:
     pool=np.arange(residue,N,q,dtype=np.int64);parts.append(rng.choice(pool,size=int(count),replace=False))
   s=np.concatenate(parts);all_counts &= np.array_equal(np.bincount(s%q,minlength=q),target)
 out.append(f'max_summary_mean_reconstruction_error={max_mean_err:.3e}')
 out.append(f'max_summary_std_reconstruction_error={max_std_err:.3e}')
 out.append(f'sampler_preserves_residue_counts={all_counts}')
 out.append(f'all_structured_values_below_all_1000_matched_samples={all_below}')
 # Explain fractions and basis stability.
 for k in (1,2,3):
  sub=[r for r in summary if int(r['k'])==k]
  best=max(sub,key=lambda r:float(r['cardinality_deficit_reduction_position']))
  out.append(f"k{k}_best_matched_bits={best['matched_low_bits']}")
  out.append(f"k{k}_position_deficit_reduction={float(best['cardinality_deficit_reduction_position']):.6f}")
  out.append(f"k{k}_fourier_deficit_reduction={float(best['cardinality_deficit_reduction_fourier']):.6f}")
  out.append(f"k{k}_position_residual_bits={float(best['position_residual_deficit']):.6f}")
  out.append(f"k{k}_fourier_residual_bits={float(best['fourier_residual_deficit']):.6f}")
 (ROOT/'RESIDUE_VALIDATION.txt').write_text('\n'.join(out)+'\n')
 print('\n'.join(out))
 if max(max_violation.values())>1e-10 or max_direct>1e-10 or single_class_violation>1e-10 or max_mean_err>1e-12 or max_std_err>1e-12 or not all_counts or not all_below: raise SystemExit(1)
if __name__=='__main__':main()
