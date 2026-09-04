import numpy as np, pickle
from scipy.interpolate import RegularGridInterpolator as RGI
from core import *
G=pickle.load(open('grid.pkl','rb'))
d=load(); pop=pops(d)
idx=G['idx']; ntr=G['ntr']
tr=d.iloc[idx[:ntr]].reset_index(drop=True); te=d.iloc[idx[ntr:]].reset_index(drop=True)
rng=np.random.default_rng(2025)
cal=tr.iloc[rng.choice(len(tr),200,replace=False)].reset_index(drop=True)
obs2=cal.PB_C2.values; obsc=cal.PB_C3.values
S,C,K=G['S'],G['C'],np.log(G['K'])
f2=RGI((S,C,K),np.log(np.clip(G['t2'],1e-12,None)),bounds_error=False,fill_value=None)
fc=RGI((S,C,K),np.log(np.clip(G['cord'],1e-12,None)),bounds_error=False,fill_value=None)
LO=[S[0],C[0],K[0]]; HI=[S[-1],C[-1],K[-1]]
def logpost(th):
    s,c,k,ls2,lsc=th
    if not(LO[0]<=s<=HI[0] and LO[1]<=c<=HI[1] and LO[2]<=k<=HI[2] and -3<ls2<1.5 and -3<lsc<1.5): return -np.inf
    p2=f2([[s,c,k]])[0]; pc=fc([[s,c,k]])[0]
    s2,sc=np.exp(ls2),np.exp(lsc)
    ll=np.sum(-0.5*((p2-np.log(obs2))/s2)**2-ls2)+np.sum(-0.5*((pc-np.log(obsc))/sc)**2-lsc)
    lp=-0.5*((s-21.)/2.5)**2-0.5*(c/1.0)**2-0.5*(ls2/0.7)**2-0.5*(lsc/0.9)**2
    return ll+lp
def chain(start,n,seed):
    r=np.random.default_rng(seed); th=np.array(start,float); lp=logpost(th)
    cov=np.diag([0.6,0.03,0.15,0.05,0.05])**2; sc=1.0; out=np.zeros((n,5)); acc=0
    for i in range(n):
        prop=th+r.multivariate_normal(np.zeros(5),cov*sc)
        lpn=logpost(prop)
        if np.log(r.random())<lpn-lp: th,lp=prop,lpn; acc+=1
        out[i]=th
        if (i+1)%100==0:
            a=acc/(i+1); sc*=np.exp((a-0.25)*1.2)
            if i>800: cov=np.cov(out[i-800:i+1].T)*2.38**2/5+1e-10*np.eye(5)
    return out,acc/n
NB,NK=4000,16000
starts=[[18,0.03,np.log(0.10),np.log(.45),np.log(.65)],
        [22,0.12,np.log(0.25),np.log(.55),np.log(.75)],
        [25,0.01,np.log(0.05),np.log(.40),np.log(.60)],
        [19,0.25,np.log(0.40),np.log(.60),np.log(.80)]]
ch=[];accs=[]
for k2,st in enumerate(starts):
    c,a=chain(st,NB+NK,200+k2); ch.append(c[NB:]); accs.append(a)
ch=np.array(ch)
ch[:,:,2]=np.exp(ch[:,:,2])
def rhat(x):
    m,n=x.shape; B=n*x.mean(1).var(ddof=1); W=x.var(1,ddof=1).mean()
    return np.sqrt(((n-1)/n*W+B/n)/W)
def ess(x):
    m,n=x.shape; v=x.var(); s=0
    for l in range(1,200):
        rl=np.mean([np.corrcoef(x[i,:-l],x[i,l:])[0,1] for i in range(m)])
        if rl<0.02: break
        s+=rl
    return m*n/(1+2*s)
names=['STBR','CABR','KDpl','log_sigma_T2','log_sigma_cord']
print('acceptance %.3f'%np.mean(accs))
post={}
for j,nm in enumerate(names):
    x=ch[:,:,j]; h=np.concatenate([x[:,:NK//2],x[:,NK//2:]],0)
    q=np.percentile(x,[2.5,50,97.5]); post[nm]=(q[1],q[0],q[2])
    print('%-15s median %8.4f  95%%CrI [%.4f, %.4f]  Rhat %.4f  ESS %.0f'%(nm,q[1],q[0],q[2],rhat(h),ess(x)))
pickle.dump(dict(ch=ch,post=post),open('post.pkl','wb'))
