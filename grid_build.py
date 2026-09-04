import numpy as np, pickle, os, time, sys
from core import *
d=load(); pop=pops(d)
rng=np.random.default_rng(2025)
idx=rng.permutation(len(d)); ntr=int(.8*len(d))
tr=d.iloc[idx[:ntr]].reset_index(drop=True); te=d.iloc[idx[ntr:]].reset_index(drop=True)
cal=tr.iloc[rng.choice(len(tr),200,replace=False)].reset_index(drop=True)
S_G=np.arange(16,29,1.0)
C_G=np.array([0.0,0.005,0.01,0.02,0.035,0.05,0.07,0.10,0.15,0.22,0.32,0.45])
K_G=np.logspace(np.log10(0.02),np.log10(1.2),12)
F='grid.pkl'
if os.path.exists(F):
    G=pickle.load(open(F,'rb'))
else:
    G=dict(S=S_G,C=C_G,K=K_G,t2=np.full((len(S_G),len(C_G),len(K_G),len(cal)),np.nan),
           cord=np.full((len(S_G),len(C_G),len(K_G),len(cal)),np.nan),
           idx=idx,ntr=ntr,cal_i=cal.index.values)
budget=float(sys.argv[1]); t0=time.time(); n=0
for i,S in enumerate(S_G):
    for j,C in enumerate(C_G):
        for k,K in enumerate(K_G):
            if np.isfinite(G['t2'][i,j,k,0]): continue
            o=run(cal,pop,S,C,K,dt=0.1)
            G['t2'][i,j,k]=o['t2']; G['cord'][i,j,k]=o['cord']; n+=1
            if time.time()-t0>budget:
                pickle.dump(G,open(F,'wb')); print('saved',n,'new; remaining',int(np.isnan(G['t2'][:,:,:,0]).sum()),flush=True); sys.exit(0)
pickle.dump(G,open(F,'wb')); print('COMPLETE',n,'new',flush=True)
