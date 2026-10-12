# -*- coding: utf-8 -*-
"""
Precompute per-individual predictions over the parameter grid.

Grid: 13 STBR x 24 CABR x 16 KD_pl = 4,992 nodes.
The KD_pl axis is LOG-spaced; any integration over this grid must use
non-uniform quadrature weights (see quadrature.py).

Run: python grid_build.py <seconds>     # resumable; repeat until COMPLETE
"""
import numpy as np, pickle, os, time, sys
from core import *

d = load(); pop = pops(d)
rng = np.random.default_rng(2025)
idx = rng.permutation(len(d)); ntr = int(.8*len(d))
tr = d.iloc[idx[:ntr]].reset_index(drop=True)
te = d.iloc[idx[ntr:]].reset_index(drop=True)
cal = tr.iloc[rng.choice(len(tr), 200, replace=False)].reset_index(drop=True)

S_G = np.arange(16, 29, 1.0)                                   # 13 nodes, weeks FA
C_G = np.linspace(0.0, 1.2, 24)                                # 24 nodes, week^-1
K_G = np.logspace(np.log10(0.008), np.log10(1.2), 16)          # 16 nodes, L/day/dm2

F = 'grid.pkl'
if os.path.exists(F):
    G = pickle.load(open(F, 'rb'))
else:
    shape = (len(S_G), len(C_G), len(K_G), len(cal))
    G = dict(S=S_G, C=C_G, K=K_G,
             t2=np.full(shape, np.nan), cord=np.full(shape, np.nan),
             idx=idx, ntr=ntr, cal_i=cal.index.values)

budget = float(sys.argv[1]); t0 = time.time(); n = 0
for i, S in enumerate(S_G):
    for j, C in enumerate(C_G):
        for k, K in enumerate(K_G):
            if np.isfinite(G['t2'][i, j, k, 0]):
                continue
            o = run(cal, pop, S, C, K, dt=0.1)
            G['t2'][i, j, k] = o['t2']; G['cord'][i, j, k] = o['cord']; n += 1
            if time.time() - t0 > budget:
                pickle.dump(G, open(F, 'wb'))
                print('saved', n, 'new; remaining',
                      int(np.isnan(G['t2'][:, :, :, 0]).sum()), flush=True)
                sys.exit(0)
pickle.dump(G, open(F, 'wb'))
print('COMPLETE', n, 'new', flush=True)
