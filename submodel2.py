# -*- coding: utf-8 -*-
"""
Sub-model 2: renal excretion and lactational transfer.

IMPLEMENTATION NOTE -- read before using.

This module does NOT integrate a postpartum ODE system.  It evaluates the
renal and lactational pathways from the Sub-model 1 solution:

  * the maternal plasma Pb trajectory is taken from Sub-model 1 (pbpk2.integrate);
  * the urinary transfer coefficient alpha_urine at each visit is the ratio of the
    observed creatinine-adjusted urinary Pb to the model-predicted plasma Pb at the
    corresponding gestational age;
  * the renal Pb clearance implied by that coefficient is compared against the
    clearance the physiological model predicts independently (GFR x CL_REN_COEF);
  * the milk-to-plasma partition coefficient P_c,milk is the ratio of the observed
    breast milk Pb geometric mean to the predicted maternal plasma Pb at delivery.

Renal elimination IS inside the Sub-model 1 mass balance (the R_ren term in
pbpk2.derivs).  Lactational elimination is NOT: milk Pb is obtained from maternal
plasma by partition after the fact, and no mass is removed from the maternal
compartments on account of lactation.  The calculation is therefore sequential
rather than jointly calibrated.  State this plainly in any text that describes it.

Reads  : grid.pkl, post.pkl
Writes : sub2.pkl
Run    : python submodel2.py
"""
import numpy as np, pickle
import pbpk2 as M
import core
from core import load, pops, build_P, GM

# ---- observed inputs (NOEMOC; see Table S2 of the article) -------------------
U_T1_CR   = 0.211      # ug/g Cr, creatinine-adjusted urinary Pb, T1  (n = 228)
U_T2_CR    = 0.322     # ug/g Cr, creatinine-adjusted urinary Pb, T2  (n = 163)
GA_T1, GA_T2 = 11.0, 25.5   # weeks, mean sampling age of the two urine visits
E_CREAT   = 1.2        # g creatinine/day, adult women
MILK_GM   = 0.1066     # ug/L, breast milk Pb geometric mean (n = 72)
V_MILK    = 0.45       # L/day, milk intake, postpartum weeks 1-4
BW_INFANT = 4.0        # kg
HCT_PP    = 0.40       # postpartum maternal haematocrit (non-pregnant value)

P = pickle.load(open('post.pkl', 'rb'))['post']
S, C, K = P['mode'] if 'mode' in P else (P['STBR'][0], P['CABR'][0], P['KDpl'][0])

d = load(); pop = pops(d)
core.run(d, pop, S, C, K, dt=0.05)                 # populates the steady-state cache
Pf = build_P(d, pop, S, C, K)
Yss = core._ss[(id(d), len(d))]
DI = d.PB.values/np.maximum(Yss[:, 0]/M.v_blood(0.0), 1e-12)
Pf['DI'] = DI
r = M.integrate(Yss*DI[:, None], Pf, 12.0, 40.0, 0.05)
t, Y = r['t'], r['Y']

def plasma_at(ga):
    i = int(np.argmin(abs(t - ga)))
    return GM(Y[i, :, 0]/M.v_blood(t[i])/M.rbp(M.hct_mat(t[i])))

# ---- renal pathway ------------------------------------------------------------
print('Sub-model 2 -- renal')
rows = []
for label, ga, u_cr in [('T1', GA_T1, U_T1_CR), ('T2', GA_T2, U_T2_CR)]:
    cp = plasma_at(ga)
    alpha = u_cr/cp
    cl_obs = alpha*E_CREAT                                    # L/day implied by the urine
    cl_mod = M.gfr_week(ga, pop['bun'], pop['bun'])/7.0*M.CL_REN_COEF   # L/day predicted
    print('  %s  plasma %.4f ug/L   alpha_urine %.2f   implied CL %.2f L/day   '
          'model CL %.2f L/day   obs/pred %.2f'
          % (label, cp, alpha, cl_obs, cl_mod, cl_obs/cl_mod))
    rows.append(dict(visit=label, GA=ga, plasma=cp, alpha=alpha,
                     CL_obs=cl_obs, CL_mod=cl_mod, ratio=cl_obs/cl_mod))

# ---- lactational pathway ------------------------------------------------------
bll40 = GM(Y[-1, :, 0]/M.v_blood(40.0))
cp_pp = bll40/M.rbp(HCT_PP)
pc_milk = MILK_GM/cp_pp
intake = MILK_GM*V_MILK

print('\nSub-model 2 -- lactation')
print('  maternal blood Pb at delivery   %.3f ug/L' % bll40)
print('  postpartum maternal plasma Pb   %.4f ug/L' % cp_pp)
print('  P_c,milk                        %.3f' % pc_milk)
print('  infant intake                   %.4f ug/day = %.4f ug/kg bw/day'
      % (intake, intake/BW_INFANT))
print('  milk Pb is %.0f%% below concurrent maternal whole-blood Pb'
      % ((1 - MILK_GM/bll40)*100))
print('\n  NOTE: the same infant intake follows directly from the observed milk Pb and '
      'the 0.45 L/day milk intake, independently of the model.')

pickle.dump(dict(renal=rows, bll40=bll40, cp_pp=cp_pp, pc_milk=pc_milk,
                 intake=intake, intake_per_kg=intake/BW_INFANT),
            open('sub2.pkl', 'wb'))
print('\nwritten: sub2.pkl')
