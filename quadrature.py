# -*- coding: utf-8 -*-
"""
Posterior evaluation by direct grid quadrature.

With only three structural parameters (STBR, CABR, KD_pl) the posterior can be
integrated exactly up to grid resolution, so no sampler is required and
convergence does not arise as a question.  The two observation-model standard
deviations (sigma_T2, sigma_cord) are marginalised numerically at every node.

Reads  : grid.pkl   (written by grid_build.py)
Writes : post.pkl   {'post': {name: (median, lo, hi)}, 'grid': ..., 'logpost': ...}

Run    : python quadrature.py
"""
import numpy as np, pickle
from core import load, GM

G = pickle.load(open('grid.pkl', 'rb'))
d = load()
tr = d.iloc[G['idx'][:G['ntr']]].reset_index(drop=True)
rng = np.random.default_rng(2025)
cal = tr.iloc[rng.choice(len(tr), 200, replace=False)].reset_index(drop=True)

o2 = np.log(cal.PB_C2.values)          # maternal T2 blood Pb, log scale
oc = np.log(cal.PB_C3.values)          # cord blood Pb, log scale
n  = len(o2)

S, C, K = G['S'], G['C'], G['K']
lt = np.log(np.clip(G['t2'],   1e-12, None))
lc = np.log(np.clip(G['cord'], 1e-12, None))

# ---- numerical marginalisation of the observation-model sigma -----------------
SIG = np.exp(np.linspace(np.log(0.15), np.log(3.0), 140))

def log_marginal(res, prior_scale):
    """log integral over sigma of a log-normal likelihood with half-normal prior."""
    SS = (res**2).sum(-1)[..., None]
    ll = (-n*np.log(SIG) - 0.5*SS/SIG**2
          - 0.5*(np.log(SIG)/prior_scale)**2 - np.log(SIG))   # last term: d sigma -> d log sigma
    m = ll.max(-1, keepdims=True)
    return m[..., 0] + np.log(np.trapz(np.exp(ll - m), SIG, axis=-1))

L = log_marginal(lt - o2, 0.7) + log_marginal(lc - oc, 0.9)

# ---- structural priors --------------------------------------------------------
#   STBR ~ N(21.0, 2.5^2)
#   CABR ~ half-normal(1.0), flat over the grid support above zero
#   KD_pl ~ log-uniform over the grid support  -> constant in log K
logpost = (L
           - 0.5*((S[:, None, None] - 21.0)/2.5)**2
           - 0.5*(C[None, :, None]/1.0)**2)
logpost -= logpost.max()
w = np.exp(logpost)

# ---- quadrature weights: trapezoid in S, C and **log K** ----------------------
def trap_w(x):
    dw = np.zeros_like(x, dtype=float)
    dw[1:-1] = (x[2:] - x[:-2])/2
    dw[0]    = (x[1] - x[0])/2
    dw[-1]   = (x[-1] - x[-2])/2
    return dw

W = (trap_w(S)[:, None, None]
     * trap_w(C)[None, :, None]
     * trap_w(np.log(K))[None, None, :])      # log spacing on the KD_pl axis

post = w*W
post /= post.sum()

# ---- marginals, medians, credible intervals -----------------------------------
def marginal(axis, vals):
    other = tuple(a for a in range(3) if a != axis)
    m = post.sum(axis=other); m = m/m.sum()
    cdf = np.cumsum(m)
    q = lambda p: float(np.interp(p, cdf, vals))
    return q(0.5), q(0.025), q(0.975), m

out = {}
print('Posterior marginals (grid quadrature)')
for ax, (nm, vals) in enumerate([('STBR', S), ('CABR', C), ('KDpl', K)]):
    md, lo, hi, m = marginal(ax, vals)
    out[nm] = (md, lo, hi)
    edge = m[-1]/m.max()
    print('  %-6s median %8.4f   95%% CrI [%.4f, %.4f]   density at upper grid edge %.3f of max'
          % (nm, md, lo, hi, edge))
    if edge > 0.01:
        print('          WARNING: upper tail is not negligible at the grid boundary; '
              'posterior mean and SD would be truncation-dependent. Widen the grid '
              'before reporting moments for this parameter.')

i, j, k = np.unravel_index(post.argmax(), post.shape)
print('Joint posterior mode: STBR=%.2f  CABR=%.4f  KDpl=%.4f' % (S[i], C[j], K[k]))
out['mode'] = (float(S[i]), float(C[j]), float(K[k]))

pickle.dump(dict(post=out, logpost=logpost, weights=W, grid=dict(S=S, C=C, K=K)),
            open('post.pkl', 'wb'))
print('written: post.pkl')
