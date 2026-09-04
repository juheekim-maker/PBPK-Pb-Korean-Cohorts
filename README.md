# Extended perinatal lead PBPK model — corrected implementation

Python implementation of an extended perinatal physiologically based pharmacokinetic (PBPK)
model for lead (Pb), spanning pregnancy through lactation, calibrated against two Korean
birth cohorts (Ko-CHENS, n = 3,625 complete triads; NOEMOC, n = 242).

Accompanies the article submitted to *Toxicology and Applied Pharmacology*
(manuscript TAAP-D-26-01261). Released under the MIT licence.

---

## Quick start

```bash
pip install -r requirements.txt
python grid_build.py 600      # precompute predictions over the parameter grid
python evaluate.py            # training / hold-out performance at the posterior mode
```

`grid_build.py` takes a time budget in seconds and resumes where it left off, so it can be
run repeatedly until the grid is complete.

## Files

| File | Purpose |
|---|---|
| `pbpk2.py` | Sub-model 1 ODE system (4 states), physiology, fixed-step RK4 integrator, pre-pregnancy steady state |
| `core.py` | Data loading, individual parameterisation, reverse dosimetry, performance metrics |
| `grid_build.py` | Precomputes per-individual predictions across the parameter grid |
| `mcmc2.py` | Metropolis–Hastings cross-check on the grid surrogate |
| `evaluate.py` | Training / hold-out evaluation and MEM-PBPK BLUP |

## Model structure

Two sequentially coupled sub-models sharing one maternal structure:

- **Sub-model 1 (pregnancy)** — maternal blood/plasma, cortical and trabecular bone,
  kidney, placenta and fetus. Placental transfer follows Fick's first law, driven by the
  maternal–fetal plasma gradient and scaled by syncytiotrophoblast exchange area and the
  inverse of diffusion path length.
- **Sub-model 2 (lactation)** — the same maternal structure with placenta and fetus removed
  and a mammary compartment added. Urinary and milk elimination are true clearances within
  the maternal mass balance.

The two are coupled through the complete maternal compartment burden vector at delivery,
not through blood Pb alone, so total maternal Pb mass is conserved across the transition.
Mass balance is verified at every time step (maximum relative error 9e-5).

## Six defects corrected relative to the originally submitted code

1. **Daily-intake back-calculation** was quadratic in T1 blood Pb. It is now taken from the
   model's own pre-pregnancy steady state and is linear in T1 blood Pb. Reverse-dosimetry
   intake rises from an implausible 0.7 to 3.5 ug/day.
2. **Bone volumes** were fixed at a 60 kg linearisation (0.36 / 0.144 L). They now follow
   O'Flaherty (2000), total bone volume = 0.0168 x BW^1.188 with an 80/20 cortical–trabecular
   split — roughly a fourfold increase in skeletal capacity.
3. **Blood-to-plasma Pb ratio** used two undocumented and unequal constants (13 maternal,
   20 fetal), which fixed the equilibrium feto-maternal ratio by construction. Both are now
   derived from one erythrocyte partition constant applied to haematocrit.
4. **Glomerular filtration rate** was entered as 0.105 L/h, but 0.105 is the ICRP value in
   litres per MINUTE. Renal Pb clearance was 60-fold too low. Corrected to 6.30 L/h.
5. **Bone remodelling** — the additional pregnancy term BRRP multiplied bone FORMATION, so
   skeletal uptake and release scaled together and no net mobilisation of stored Pb was
   possible; increasing CABR lowered rather than raised maternal blood Pb. Formation and
   resorption are now separate fluxes and BRRP acts on resorption.
6. **Integration** used adaptive LSODA with rtol 5e-3, a 300-step ceiling and a bare
   exception handler that silently replaced failed solutions. It is now fixed-step RK4 at
   0.05 weeks with mass balance verified at every step.

## Inference

The posterior is evaluated by direct grid quadrature over 4,992 nodes
(13 STBR x 24 CABR x 16 KD_pl), with the two observation-model standard deviations
marginalised numerically at each node. With three structural parameters this is exact up to
grid resolution and removes sampler convergence from consideration. A Metropolis–Hastings
run (`mcmc2.py`) is retained as a cross-check; it reaches the same posterior mode but mixes
poorly along the STBR–CABR ridge, which is itself evidence that those two parameters are
only jointly identified.

## Results

Posterior marginal medians (95% credible interval):

| Parameter | Median | 95% CrI |
|---|---|---|
| STBR (weeks fetal age) | 24.45 | 22.66–26.86 |
| CABR (week^-1) | 0.364 | 0.223–1.022 |
| KD_pl (L/day/dm2) | 0.0435 | 0.018–0.060 |

Joint posterior mode used for deterministic prediction:
STBR 24.0 wk, CABR 0.260 /wk, KD_pl 0.0611 L/day/dm2.

| | Training (n = 2,863) | Hold-out (n = 716) |
|---|---|---|
| Maternal T2 blood Pb, MFE / W2fold | 1.007 / 80.5% | 1.039 / 84.2% |
| Cord blood Pb, MFE / W2fold | 1.046 / 73.2% | 1.054 / 75.6% |
| Predicted f/m ratio (observed) | 0.785 (0.756) | 0.785 (0.774) |

The observed feto-maternal ratio is reproduced with no post hoc rescaling of any parameter.

## Data

The Ko-CHENS and NOEMOC datasets are not redistributed here. Ko-CHENS data are held by the
National Institute of Environmental Research (NIER), Korea, and are available on formal
request subject to institutional review; NOEMOC data are available from the corresponding
author on reasonable request subject to IRB approval. See the Data Availability statement
in the article. `core.py` expects the Ko-CHENS file at the path set in `KC`.

## Citation

See `CITATION.cff`. Please cite both this software (Zenodo DOI) and the article.
