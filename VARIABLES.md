# VARIABLES.md — Variable dictionary

Extended perinatal lead (Pb) toxicokinetic model
Accompanies Kim JH, Seo SC, Toxicology and Applied Pharmacology (TAAP-D-26-01261).
Written against source version 1.1.0 (`pbpk2.py`, `core.py`, `grid_build.py`, `quadrature.py`, `mcmc2.py`, `submodel2.py`, `evaluate.py`).

---

## 0. How to read this file

Every symbol in the source is listed with its code name, units, defining expression and the
corresponding equation in the article. If a symbol appears in the article but not here, the
article is describing something the released code does not implement; see §9.

**Units.** Pb amounts in µg, concentrations in µg/L, volumes in L, areas in dm², lengths in dm.
**Internal time is in weeks of gestational age (GA)**, and every rate in `derivs` is therefore
per week. Daily quantities from the article (intake in µg/day, clearance in L/day) are
converted at the point of use by the factors `7.0` and `168.0` that appear in the source.

**Age convention.** `GA` = gestational age in weeks. `FA` = fetal age = `GA − 2.0`
(`pbpk2.derivs`, line 59). **STBR is in weeks FA**; the integration variable is GA. Confusing
the two shifts the bone-mobilisation window by two weeks and is the easiest error to make here.

---

## 1. State vector

`Y` has shape `(n_individuals, 4)`. Column order is fixed throughout the code.

| Index | Symbol | Code | Units | Meaning |
|---|---|---|---|---|
| `Y[:,0]` | A_blood | `Ab` | µg | Pb in maternal blood |
| `Y[:,1]` | A_cort | `Ac` | µg | Pb in maternal cortical bone |
| `Y[:,2]` | A_trab | `At` | µg | Pb in maternal trabecular bone |
| `Y[:,3]` | A_fet | `Afb` | µg | Pb in fetal blood |

Derived concentrations:

| Symbol | Code | Expression | Units |
|---|---|---|---|
| C_blood | `Ab/Vb` | maternal whole-blood Pb | µg/L |
| C_p,mat | `Cp` | `(Ab/v_blood(GA)) / rbp(hct_mat(GA))` | µg/L |
| C_p,fet | `Cpf` | `(Afb/Vfb) / rbp(HCT_FET)`, `Vfb = max(0.08·BW_fet, 1e-4)` | µg/L |

Reported endpoints (`core.run`): maternal T2 blood Pb is taken at the integration step nearest
**GA 28**; cord blood Pb is the final fetal blood concentration at **GA 40**, using a fetal blood
volume of `0.08 × bw_fetal(38)` (FA 38 = GA 40).

---

## 2. Fixed constants (`pbpk2.py`, lines 18–26)

| Symbol | Code | Value | Units | Meaning |
|---|---|---|---|---|
| BSF | `BSF_REF` | 1.73 | m² | Reference body surface area; multiplies both `v_blood` and `gfr_week` |
| K_rbc | `K_RBC` | 30.0 | — | Erythrocyte-to-plasma Pb partition constant, applied on both maternal and fetal sides |
| Hct_np | `HCT_NP` | 0.40 | — | Non-pregnant maternal haematocrit |
| Hct_term | `HCT_TERM` | 0.33 | — | Maternal haematocrit at term |
| Hct_fet | `HCT_FET` | 0.50 | — | Fetal haematocrit (constant) |
| f_abs | `F_ABS` | 0.15 | — | GI absorption fraction, constant across gestation |
| k_bone | `KF_BONE` | 0.001 | week⁻¹ | Baseline fractional bone turnover. *Written as `0.001/168.0*168.0`, which evaluates to 0.001; the round-trip conversion is a no-op and should be simplified.* |
| K_part,bone | `K_BONE_PART` | 15000.0 | µg/L per µg/L | Bone–plasma partition constant governing uptake into newly formed bone. *Not cited to a source in the code; supply the reference or state that it was set to reproduce the adult bone:blood ratio.* |
| GFR_ref | `GFR_LPH` | 6.30 | L/h | Non-pregnant GFR (ICRP 0.105 L/min) |
| k_GFR,Pb | `CL_REN_COEF` | 0.0269 | — | Fraction of GFR constituting renal Pb clearance; anchored to the NOEMOC T1 urine series and **not re-estimated** (article §4.6, limitation 7) |
| k_fet,loss | `K_FET_LOSS` | 0.0336 | week⁻¹ | First-order loss from fetal blood (written `0.0002*168.0`) |

---

## 3. Physiological functions (`pbpk2.py`, lines 28–54)

| Symbol | Code | Argument | Returns | Expression |
|---|---|---|---|---|
| Hct_mat(GA) | `hct_mat` | GA | — | `0.40 + (0.33 − 0.40)·clip(GA/40, 0, 1)` |
| R_bp(Hct) | `rbp` | Hct | — | `(1 − Hct) + Hct·K_RBC`. Maternal at term = **10.57**; fetal = **15.50**. The equilibrium feto-maternal ratio is bounded by the ratio of these two factors. |
| BW_fet(FA) | `bw_fetal` | FA | kg | Gompertz: `0.00137·exp((0.1974/0.01306)·(1 − exp(−0.01306·FA)))` |
| S_sct(GA) | `s_sct` | GA | dm² | `max(0.15 + 0.85·clip(GA/40,0,1), 0.15)` → 0.15 at conception, 1.00 at term |
| T_sct(GA) | `t_sct` | GA | dm | `max(0.003 − 0.00291·clip(GA/40,0,1), 0.00009)` → 0.003 early, 0.00009 at term |
| V_cort, V_trab | `bone_volumes` | BW (kg) | L | `V_total = 0.0168·BW^1.188` (O'Flaherty 2000); returns `0.80·V_total`, `0.20·V_total` |
| V_blood(GA) | `v_blood` | GA | L | `4.70·BSF_REF·(1 + 0.45·clip(GA/40,0,1))` → 8.13 L non-pregnant, 11.79 L at term |
| GFR(GA, BUN) | `gfr_week` | GA, BUN | L/week | `6.30·BSF_REF·(1 + 0.50·clip(GA/40,0,1))·168 · exp(0.017·(BUN − BUN̄)/2.0)` |

> **Check before release.** `GFR_LPH = 6.30 L/h` is already the ICRP value standardised to
> 1.73 m², and `v_blood` likewise takes 4.70 L as a per-1.73 m² figure. Both are then multiplied
> by `BSF_REF = 1.73` again. Confirm this is intended. The renal clearance that results
> (`GFR × 0.0269` = 8.1 L/day at GA 12 and 9.5 L/day at GA 28) matches the values reported in
> the article, so the calibrated `CL_REN_COEF` absorbs whatever scaling convention is in force;
> the point is that the convention should be stated rather than left implicit.

---

## 4. Right-hand side (`pbpk2.derivs`, lines 57–94)

| Symbol | Code | Units | Expression | Article |
|---|---|---|---|---|
| BRRP | `BRRP` | — | `CABR·(FA − STBR)` for FA ≥ STBR, else 0; floored at 0 | Eq. 2 |
| BFR_c, BFR_t | `BFR_c`, `BFR_t` | L/week | `k_bone·V_cort`, `k_bone·V_trab` — bone **formation** | Eq. 3 |
| BRR_c, BRR_t | `BRR_c`, `BRR_t` | L/week | `k_bone·V_cort·(1 + BRRP)`, same for trabecular — bone **resorption** | Eq. 3 |
| R_up | `R_up` | µg/week | `K_BONE_PART·C_p·(BFR_c + BFR_t)` — plasma → bone | — |
| R_rel | `R_rel` | µg/week | `(Ac/V_cort)·BRR_c + (At/V_trab)·BRR_t` — bone → plasma | — |
| R_ren | `R_ren` | µg/week | `GFR_week·C_p·CL_REN_COEF·albf` | §2.3.6 |
| R_abs | `R_abs` | µg/week | `DI·F_ABS·7.0` | Eq. 4 |
| RD_pl | `RD` | µg/week | `KD_pl·(C_p,mat − C_p,fet)·(0.00009/T_sct(GA))·S_sct(GA)·7.0`, zero before FA 0 | Eq. 1 |
| bf | `bf` | — | Fetal ossification ramp: 0 before FA 13, `(FA−13)/13` to FA 26, 1.0 thereafter | §2.3.2 |
| R_fb | `Rfb` | µg/week | `K_BONE_PART·C_p,fet·0.001·0.07·BW_fet·bf + K_FET_LOSS·Afb` — fetal bone deposition plus first-order loss | — |

Balance equations:

```
dA_blood/dGA = R_abs + R_rel − R_up − R_ren − RD
dA_cort /dGA = K_BONE_PART·C_p·BFR_c − (Ac/V_cort)·BRR_c
dA_trab /dGA = K_BONE_PART·C_p·BFR_t − (At/V_trab)·BRR_t
dA_fet  /dGA = RD − R_fb
```

**Sign convention.** Formation and resorption are separate fluxes and BRRP multiplies
**resorption only**, so raising CABR produces net mobilisation of skeletal Pb into maternal
plasma. An earlier version multiplied formation, which made increasing CABR *lower* maternal
blood Pb; see README defect 5.

---

## 5. Integration (`pbpk2.py`, lines 97–128)

| Item | Code | Setting |
|---|---|---|
| Integrator | `_rk4`, `integrate` | Fixed-step classical RK4; states floored at 0 after each step |
| Interval | `core.run` | GA 12 → 40 weeks |
| Step | `dt` | 0.05 weeks in `evaluate.py`; 0.1 in `grid_build.py`; 0.02 is the `core.run` default |
| Mass balance | `track_mass=True` | `cum_in = ∫R_abs`, `cum_out = ∫(R_ren + R_fb)`; relative error = \|Y₀ + cum_in − cum_out − Y_final\| / Y₀, reported as `mberr` |
| Pre-pregnancy steady state | `steady_state_unit_intake` | 40 years at `dt = 0.25` weeks with `DI = 1 µg/day` and BRRP inactive |

Placental transfer `RD` moves Pb between two tracked states and therefore does not appear in
`cum_out`; fetal loss `R_fb` does.

---

## 6. Individual parameterisation and reverse dosimetry (`core.py`)

`build_P` (lines 15–21) returns the per-individual parameter dictionary:

| Key | Expression | Meaning |
|---|---|---|
| `STBR` | scalar, not individualised | Onset of additional remodelling, weeks FA |
| `Vc`, `Vt` | `bone_volumes(BW)` | Cortical and trabecular bone volume, L |
| `CABR` | `CABR · exp(0.024·(Ca_i − Ca̅)/0.34)` | Serum-calcium scaling; 0.34 is the assumed SD of serum Ca |
| `KDpl` | `KD_pl · exp(0.012·(ln hs-CRP_i − ln hs-CRP_GM))` | Equivalent to the power-law form `(hs-CRP_i/GM)^0.012` in the article |
| `GFRw` | `gfr_week(28.0, BUN_i, BUN̄)` | **Evaluated once at GA 28 and held constant** for the whole simulation — the gestational GFR trajectory is not time-varying in this version |
| `albf` | `exp(0.016·(ALB_i − ALB̅)/0.25)` | Albumin scaling. **Applied to renal clearance only** (`R_ren`); it does not modify the free fraction available for placental diffusion |

`run` (lines 23–41) performs reverse dosimetry: the unit-intake steady state is computed once
per dataset and cached (`_ss`), the implied blood concentration `Cb` is read off, and

```
DI_i = PB_T1,i / Cb_unit          (µg/day, linear in T1 blood Pb)
Y0_i = Y_ss,unit · DI_i
```

so initial burdens in all four compartments scale with individual intake.

`metrics` (lines 42–47): `MFE = exp(mean(log(pred/obs)))`, `MPE = mean(pred/obs − 1)·100`,
`W2 = proportion with |log(pred/obs)| < log 2`, `r` = Pearson correlation on logs.

Smoking (`COTIN`) is loaded but **not used inside the model**, by design: its effect on intake
is already embedded in the T1 blood Pb from which `DI` is back-calculated.

---

## 7. Parameter grid and inference

`grid_build.py` precomputes predictions on a three-dimensional grid and caches them to
`grid.pkl`, resuming across runs under a wall-clock budget given as `argv[1]`.

| Axis | Code | Nodes | Values |
|---|---|---|---|
| STBR | `S_G` | 13 | `arange(16, 29, 1.0)` — 16 to 28 weeks FA |
| CABR | `C_G` | 24 | `linspace(0.0, 1.2, 24)` week⁻¹ |
| KD_pl | `K_G` | 16 | `logspace(log10(0.008), log10(1.2), 16)` — **non-uniformly spaced** |

Total 13 × 24 × 16 = **4,992 nodes**.

> **Any integration or averaging over this grid must use non-uniform trapezoidal weights.**
> A plain mean over nodes is incorrect because the KD_pl axis is logarithmic.

### `grid.pkl` keys

| Key | Contents |
|---|---|
| `S`, `C`, `K` | Node values for the three axes |
| `t2` | Predicted maternal GA-28 blood Pb, shape `(n_S, n_C, n_K, n_cal)` |
| `cord` | Predicted cord blood Pb, same shape |
| `idx` | Permutation of dataset rows; `idx[:ntr]` is the training partition |
| `ntr` | Size of the training partition (80%) |
| `cal_i` | Row indices of the 200-individual calibration subsample |

Partitioning is reproducible: `default_rng(2025)`, 80/20 split, calibration subsample of 200
drawn from the training partition.

`quadrature.py` evaluates the posterior on this grid. At every node the two observation-model
standard deviations are marginalised numerically over a 140-point log-spaced σ axis
(`log_marginal`, half-normal priors of scale 0.7 for σ_T2 and 0.9 for σ_cord), the structural
priors are applied, and the normalised posterior is integrated with `trap_w` weights that are
trapezoidal in STBR, CABR and **log** KD_pl. It prints marginal medians and 95% credible
intervals, reports the joint posterior mode used for deterministic prediction, and warns when a
marginal retains non-negligible density at the upper grid boundary — the condition under which a
posterior mean or SD would be an artefact of where the grid was truncated rather than a property
of the posterior. Output: `post.pkl`.

`mcmc2.py` is a cross-check only. It builds `RegularGridInterpolator` surrogates on `(S, C, log K)` of the log predictions
and runs four adaptive Metropolis–Hastings chains (4,000 burn-in + 16,000 retained) over
`(STBR, CABR, log KD_pl, log σ_T2, log σ_cord)`. Priors in `logpost` (line 21):

| Parameter | Prior in code |
|---|---|
| STBR | `N(21.0, 2.5²)` |
| CABR | half-normal, scale 1.0, restricted to the grid range |
| log KD_pl | uniform over the grid range (i.e. log-uniform in KD_pl) |
| log σ_T2 | `N(0, 0.7²)` |
| log σ_cord | `N(0, 0.9²)` |

Convergence diagnostics (`rhat`, `ess`) are computed and printed. The sampler reaches the same
posterior mode as the quadrature but mixes poorly along the STBR–CABR ridge, which is reported in
the article as evidence that those two parameters are only jointly identified.

`evaluate.py` runs the full training and hold-out partitions at the chosen parameter vector with
`dt = 0.05` and mass-balance tracking, then computes the MEM-PBPK random effect and BLUP:

| Symbol | Code | Expression |
|---|---|---|
| ε_i | `eps` | `log(Cord_obs,i / Cord_pred,i)` |
| σ_rand | `sig` | `eps.std(ddof=1)` |
| σ_meas | literal `0.20` | Assumed cord blood measurement error, log scale |
| λ | `lam` | `σ_rand² / (σ_rand² + 0.20²)` |
| BLUP | `blup` | `Cord_pred · exp(λ·ε)` — **in-sample by construction**; requires an observed cord value |

The hold-out BLUP reuses λ estimated on the training partition.

---

## 8. Dataset columns expected by `core.load`

Path is hard-coded in `core.KC`; change it to point at your copy. Rows are kept only where
`PB`, `PB_C2`, `PB_C3`, `CA`, `ALB`, `BUN` are all present and the three Pb values are positive.

| Column | Meaning |
|---|---|
| `PB` | Maternal blood Pb at T1, GA 9–14 wk (µg/L) |
| `PB_C2` | Maternal blood Pb at T2, GA 26–30 wk (µg/L) |
| `PB_C3` | Cord blood Pb (µg/L) |
| `PB_C6` | Infant blood Pb at 6 months (µg/L) — loaded, not used in Sub-model 1 |
| `CA`, `ALB`, `BUN`, `HS_CRP` | T1 clinical chemistry; `HS_CRP` missing values filled at the median |
| `COTIN` | Urinary cotinine; loaded, not used inside the model |
| `D_PRE_WT` | Pre-pregnancy weight (kg) → `BW`, missing filled at the median |
| `AGE`, `D_GWG_WT`, `BMI_PRE`, `D_B1_SEX` | Loaded for the pre-PBPK analyses |

`NA` sentinels handled on load: `#N/A`, `99999`, `NA`.

Cohort data are not redistributed. Ko-CHENS access is through NIER; NOEMOC access is by
request to the corresponding author.

---

## 9. Sub-model 2 (`submodel2.py`)

Sub-model 2 does **not** integrate a postpartum ODE system. It evaluates the renal and
lactational pathways from the Sub-model 1 solution.

| Symbol | Code | Expression | Units |
|---|---|---|---|
| C_p(GA) | `plasma_at` | Geometric mean plasma Pb from the Sub-model 1 trajectory at the given GA | µg/L |
| α_urine | `alpha` | `U_cr / C_p(GA_visit)`, one per visit | (µg/g Cr)/(µg/L) |
| CL_urine,implied | `cl_obs` | `α_urine · E_CREAT`, with `E_CREAT = 1.2` g/day | L/day |
| CL_urine,model | `cl_mod` | `gfr_week(GA)/7 · CL_REN_COEF` — predicted independently | L/day |
| BLL at delivery | `bll40` | `GM(A_blood/v_blood(40))` | µg/L |
| C_p,postpartum | `cp_pp` | `bll40 / rbp(0.40)` | µg/L |
| P_c,milk | `pc_milk` | `MILK_GM / cp_pp` | — |
| Infant intake | `intake` | `MILK_GM · 0.45` L/day, divided by 4 kg for the per-kg figure | µg/day |

Observed inputs are module constants at the top of the file: `U_T1_CR = 0.211`,
`U_T2_CR = 0.322` µg/g Cr at `GA_T1 = 11.0` and `GA_T2 = 25.5` weeks, and `MILK_GM = 0.1066` µg/L.

**Scope of the mass balance.** Renal elimination is inside the Sub-model 1 mass balance (the
`R_ren` term in `pbpk2.derivs`). Lactational elimination is not: milk Pb is obtained by partition
after the maternal trajectory has been solved, and no mass is removed from the maternal
compartments on account of lactation. The two sub-models are therefore evaluated sequentially,
and `P_c,milk` is conditional on the Sub-model 1 solution rather than estimated jointly with it.
This is stated in Sections 2.3.1, 2.3.6 and 4.8 of the article.

---

## 10. What this release does not implement

Listed so that a reader does not search for code that is not here.

- **Soft-tissue compartments.** The maternal model has four states. Liver/kidney, poorly
  perfused, well perfused and adipose/brain compartments, and the mature/juvenile bone
  subdivision of the parent frameworks, are not present.
- **Blood flows and cardiac output.** There is no `CO(t)` and no compartmental blood-flow
  vector; transport is not flow-limited. Maternal physiology enters through blood volume,
  haematocrit and GFR only.
- **A postpartum ODE integration.** See §9.
- **Time-varying GFR within a simulation.** `GFRw` is evaluated once at GA 28 (`core.build_P`).

---

## 11. File map

| File | Purpose |
|---|---|
| `pbpk2.py` | Sub-model 1 ODE system, physiology, RK4 integrator, pre-pregnancy steady state |
| `core.py` | Data loading, individual parameterisation, reverse dosimetry, performance metrics |
| `grid_build.py` | Precomputes predictions over the 4,992-node parameter grid (resumable) |
| `quadrature.py` | Posterior evaluation by grid quadrature; marginals, credible intervals, joint mode |
| `mcmc2.py` | Metropolis–Hastings cross-check on the grid surrogate |
| `submodel2.py` | Renal and lactational pathways evaluated from the Sub-model 1 solution |
| `evaluate.py` | Training / hold-out performance and MEM-PBPK BLUP |
| `VARIABLES.md` | This file |

---

*Licence: MIT. Questions and re-use: Ju Hee Kim, juheekim@khu.ac.kr, ORCID 0000-0001-6586-7244.*
