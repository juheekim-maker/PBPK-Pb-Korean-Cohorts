import numpy as np, pandas as pd, pbpk2 as M
KC='/mnt/user-data/uploads/FINAL_KoCHENS_DATASET_20260312_Cd.xlsx'
def load():
    kc=pd.read_excel(KC,na_values=['#N/A','99999','NA'])
    for c in ['PB','PB_C2','PB_C3','PB_C6','CA','ALB','BUN','HS_CRP','COTIN','D_PRE_WT','AGE','D_GWG_WT','BMI_PRE','D_B1_SEX']:
        kc[c]=pd.to_numeric(kc[c],errors='coerce')
    d=kc.dropna(subset=['PB','PB_C2','PB_C3','CA','ALB','BUN']).copy()
    d=d[(d.PB>0)&(d.PB_C2>0)&(d.PB_C3>0)].reset_index(drop=True)
    d['HS_CRP']=d['HS_CRP'].fillna(d['HS_CRP'].median())
    d['BW']=d['D_PRE_WT'].fillna(d['D_PRE_WT'].median())
    return d
GM=lambda x: float(np.exp(np.mean(np.log(np.clip(np.asarray(x,float),1e-12,None)))))
def pops(d): return dict(ca=d.CA.mean(),alb=d.ALB.mean(),bun=d.BUN.mean(),
                         crp=np.exp(np.log(d.HS_CRP.clip(0.01)).mean()))
def build_P(df,pop,STBR,CABR,KDpl):
    Vc,Vt=M.bone_volumes(df.BW.values)
    return dict(STBR=STBR,Vc=Vc,Vt=Vt,
        CABR=CABR*np.exp(0.024*(df.CA.values-pop['ca'])/0.34),
        KDpl=KDpl*np.exp(0.012*(np.log(df.HS_CRP.clip(0.01).values)-np.log(pop['crp']))),
        GFRw=M.gfr_week(28.0,df.BUN.values,pop['bun']),
        albf=np.exp(0.016*(df.ALB.values-pop['alb'])/0.25))
_ss={}
def run(df,pop,STBR,CABR,KDpl,dt=0.02,mass=False):
    P=build_P(df,pop,STBR,CABR,KDpl)
    key=(id(df),len(df))
    if key not in _ss:                       # steady state independent of STBR/CABR/KDpl
        Q=dict(P); Q['CABR']=np.zeros(len(df)); Q['DI']=np.ones(len(df))
        _ss[key]=M.steady_state_unit_intake(Q)
    Yss=_ss[key]
    Cb=Yss[:,0]/M.v_blood(0.0)
    DI=df.PB.values/np.maximum(Cb,1e-12)
    Y0=Yss*DI[:,None]; P['DI']=DI
    r=M.integrate(Y0,P,12.0,40.0,dt,track_mass=mass)
    t,Y=r['t'],r['Y']; i=np.argmin(abs(t-28))
    t2=Y[i,:,0]/M.v_blood(28.0)
    cord=Y[-1,:,3]/max(0.08*M.bw_fetal(38.0),1e-4)
    out=dict(t2=t2,cord=cord,DI=DI,Cp_T1=Y0[:,0]/M.v_blood(12.0)/M.rbp(M.hct_mat(12.0)))
    if mass:
        e=abs(Y0.sum(1)+r['cum_in']-r['cum_out']-Y[-1].sum(1))/np.maximum(Y0.sum(1),1e-12)
        out['mberr']=float(e.max())
    return out
def metrics(pred,obs):
    p=np.clip(np.asarray(pred,float),1e-12,None); o=np.asarray(obs,float)
    lr=np.log(p/o)
    return dict(MFE=float(np.exp(lr.mean())),MPE=float((np.exp(lr)-1).mean()*100),
                W2=float((np.abs(lr)<np.log(2)).mean()*100),
                r=float(np.corrcoef(np.log(p),np.log(o))[0,1]))
