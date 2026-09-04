# -*- coding: utf-8 -*-
"""
Extended perinatal Pb PBPK, Sub-model 1 -- corrected implementation (v2).

Defects repaired relative to the originally submitted code:
 (1) daily-intake back-calculation quadratic in T1 BLL -> linear steady-state reverse dosimetry
 (2) bone volumes fixed at a 60-kg linearisation       -> O'Flaherty (2000) allometric, 80/20
 (3) blood:plasma ratio two undocumented constants     -> single K_rbc applied to haematocrit
 (4) integration rtol 5e-3 / 300-step cap / silent fail-> fixed-step RK4, mass balance verified
 (5) GFR entered as 0.105 L/h (a per-minute value)     -> 6.30 L/h  [60x]
 (6) additional remodelling multiplied bone FORMATION, so uptake and release scaled together and
     no net skeletal mobilisation was possible         -> formation and resorption separated,
                                                          CABR acts on resorption only
States per individual: [A_blood, A_cort, A_trab, A_fetal_blood]  (ug); time in weeks of GA.
"""
import numpy as np

BSF_REF   = 1.73
K_RBC     = 30.0
HCT_NP, HCT_TERM, HCT_FET = 0.40, 0.33, 0.50
F_ABS     = 0.15
KF_BONE   = 0.001/168.0*168.0     # baseline fractional bone turnover, per week
K_BONE_PART = 15000.0
GFR_LPH   = 6.30                  # non-pregnant GFR, L/h  (0.105 L/min)
CL_REN_COEF = 0.0269              # fraction of GFR effective for Pb (anchored to NOEMOC T1 urine)
K_FET_LOSS  = 0.0002*168.0

def hct_mat(GA):
    f = np.clip(GA/40.0, 0.0, 1.0)
    return HCT_NP + (HCT_TERM - HCT_NP)*f

def rbp(hct):
    return (1.0 - hct) + hct*K_RBC

def bw_fetal(FA):
    FA = np.maximum(FA, 0.01)
    return 0.00137*np.exp(0.1974/0.01306*(1 - np.exp(-0.01306*FA)))

def s_sct(GA):
    return np.maximum(0.15 + 0.85*np.clip(GA/40.0, 0, 1), 0.15)

def t_sct(GA):
    return np.maximum(0.003 - 0.00291*np.clip(GA/40.0, 0, 1), 0.00009)

def bone_volumes(BW):
    Vt = 0.0168*np.power(BW, 1.188)
    return 0.80*Vt, 0.20*Vt

def v_blood(GA, bsf=BSF_REF):
    return 4.70*bsf*(1 + 0.45*np.clip(GA/40.0, 0, 1))

def gfr_week(GA, bun, bun_mean, bsf=BSF_REF):
    base = GFR_LPH*bsf*(1 + 0.50*np.clip(GA/40.0, 0, 1))*168.0
    return base*np.exp(0.017*(bun - bun_mean)/2.0)


def derivs(GA, Y, P):
    Ab, Ac, At, Afb = Y[:, 0], Y[:, 1], Y[:, 2], Y[:, 3]
    FA = np.maximum(GA - 2.0, 0.0)

    Vb  = v_blood(GA)
    Cp  = (Ab/Vb)/rbp(hct_mat(GA))

    # --- bone: formation and resorption are separate fluxes -------------------
    # baseline turnover is balanced; pregnancy-associated remodelling (BRRP)
    # accelerates RESORPTION, which is what mobilises stored skeletal Pb.
    BRRP  = np.where(FA >= P['STBR'], P['CABR']*(FA - P['STBR']), 0.0)
    BRRP  = np.maximum(BRRP, 0.0)
    BFR_c = KF_BONE*P['Vc']                       # formation, L/week
    BFR_t = KF_BONE*P['Vt']
    BRR_c = KF_BONE*P['Vc']*(1 + BRRP)            # resorption, L/week
    BRR_t = KF_BONE*P['Vt']*(1 + BRRP)

    R_up  = K_BONE_PART*Cp*(BFR_c + BFR_t)                       # plasma -> bone
    R_rel = (Ac/P['Vc'])*BRR_c + (At/P['Vt'])*BRR_t              # bone -> plasma

    R_ren = P['GFRw']*Cp*CL_REN_COEF*P['albf']
    R_abs = P['DI']*F_ABS*7.0

    bwf = bw_fetal(FA)
    Vfb = np.maximum(0.08*bwf, 1e-4)
    Cpf = (Afb/Vfb)/rbp(HCT_FET)
    Tr  = 0.00009/t_sct(GA)
    RD  = P['KDpl']*(Cp - Cpf)*Tr*s_sct(GA)*7.0
    RD  = np.where(FA > 0, RD, 0.0)

    bf  = np.where(FA >= 13, np.minimum((FA - 13)/13.0, 1.0), 0.0)
    Rfb = K_BONE_PART*Cpf*0.001*0.07*bwf*bf + K_FET_LOSS*Afb

    dAb  = R_abs + R_rel - R_up - R_ren - RD
    dAc  = K_BONE_PART*Cp*BFR_c - (Ac/P['Vc'])*BRR_c
    dAt  = K_BONE_PART*Cp*BFR_t - (At/P['Vt'])*BRR_t
    dAfb = RD - Rfb
    return np.stack([dAb, dAc, dAt, dAfb], axis=1), dict(R_abs=R_abs, R_ren=R_ren, Rfb=Rfb)


def _rk4(t, Y, P, dt):
    k1, f1 = derivs(t, Y, P)
    k2, _  = derivs(t + dt/2, Y + dt/2*k1, P)
    k3, _  = derivs(t + dt/2, Y + dt/2*k2, P)
    k4, _  = derivs(t + dt,   Y + dt*k3,  P)
    return np.maximum(Y + dt/6*(k1 + 2*k2 + 2*k3 + k4), 0.0), f1


def integrate(Y0, P, GA0, GA1, dt, track_mass=False):
    n = int(round((GA1 - GA0)/dt))
    Y = Y0.copy(); ts=[GA0]; out=[Y.copy()]
    cin = np.zeros(Y.shape[0]); cout = np.zeros(Y.shape[0])
    for i in range(n):
        t = GA0 + i*dt
        Y, f1 = _rk4(t, Y, P, dt)
        if track_mass:
            cin  += f1['R_abs']*dt
            cout += (f1['R_ren'] + f1['Rfb'])*dt
        ts.append(t+dt); out.append(Y.copy())
    r = dict(t=np.array(ts), Y=np.array(out))
    if track_mass: r['cum_in']=cin; r['cum_out']=cout
    return r


def steady_state_unit_intake(P, years=40, dt=0.25):
    """Pre-pregnancy steady state at DI = 1 ug/day (BRRP inactive at GA = 0)."""
    N = P['Vc'].shape[0]
    Y = np.zeros((N, 4))
    Q = dict(P); Q['DI'] = np.ones(N)
    for _ in range(int(years*52/dt)):
        Y, _ = _rk4(0.0, Y, Q, dt)
    return Y
