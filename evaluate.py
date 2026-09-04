import numpy as np, pickle, pandas as pd
from core import *
G=pickle.load(open('grid.pkl','rb')); P=pickle.load(open('post.pkl','rb'))['post']
d=load(); pop=pops(d); idx=G['idx']; ntr=G['ntr']
tr=d.iloc[idx[:ntr]].reset_index(drop=True); te=d.iloc[idx[ntr:]].reset_index(drop=True)
S,C,K=P['STBR'][0],P['CABR'][0],P['KDpl'][0]
print('posterior medians: STBR=%.2f CABR=%.4f KDpl=%.4f'%(S,C,K))
res={}
for nm,df in [('Training (80%)',tr),('Hold-out test (20%)',te)]:
    o=run(df,pop,S,C,K,dt=0.05,mass=True)
    m2=metrics(o['t2'],df.PB_C2.values); mc=metrics(o['cord'],df.PB_C3.values)
    fm_p=GM(o['cord']/o['t2']); fm_o=GM(df.PB_C3/df.PB_C2)
    print('\n%s  n=%d  (mass-balance err %.1e)'%(nm,len(df),o['mberr']))
    print('  T2 BLL : MFE %.3f  MPE %+.1f%%  W2fold %.1f%%  r %.3f'%(m2['MFE'],m2['MPE'],m2['W2'],m2['r']))
    print('  Cord Pb: MFE %.3f  MPE %+.1f%%  W2fold %.1f%%  r %.3f'%(mc['MFE'],mc['MPE'],mc['W2'],mc['r']))
    print('  f/m GM : predicted %.3f  observed %.3f'%(fm_p,fm_o))
    print('  DI GM  : %.2f ug/day'%GM(o['DI']))
    res[nm]=dict(m2=m2,mc=mc,fm_p=fm_p,fm_o=fm_o,o=o,df=df)
# MEM-PBPK random effect + BLUP
o=res['Training (80%)']['o']; df=res['Training (80%)']['df']
eps=np.log(df.PB_C3.values/np.clip(o['cord'],1e-12,None))
sig=eps.std(ddof=1); lam=sig**2/(sig**2+0.20**2)
blup=np.clip(o['cord'],1e-12,None)*np.exp(lam*eps)
mb=metrics(blup,df.PB_C3.values)
print('\nMEM-PBPK: sigma_rand=%.3f  lambda=%.3f  BLUP MFE %.3f  W2fold %.1f%%'%(sig,lam,mb['MFE'],mb['W2']))
ot=res['Hold-out test (20%)']['o']; dt_=res['Hold-out test (20%)']['df']
et=np.log(dt_.PB_C3.values/np.clip(ot['cord'],1e-12,None))
bt=np.clip(ot['cord'],1e-12,None)*np.exp(lam*et)
mbt=metrics(bt,dt_.PB_C3.values)
print('          hold-out BLUP MFE %.3f  W2fold %.1f%%'%(mbt['MFE'],mbt['W2']))
pickle.dump(dict(res={k:{kk:vv for kk,vv in v.items() if kk!='df'} for k,v in res.items()},
                 sig=sig,lam=lam,mb=mb,mbt=mbt),open('eval.pkl','wb'))
