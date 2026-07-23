# GRLP literature notes

Source-grounded notes on the papers that define, correct, extend, and apply
**GRLP** (the gravel-river long-profile model). These notes exist to support the
documentation and to give a single place where the governing equations, their
one published correction, and the model's real-world applications are recorded
with their citations. They are working notes, not a substitute for the papers.

**The PDFs themselves are not tracked in git** (they are copyrighted and kept
locally under `docs/literature/`, ignored via `*.pdf`). This markdown file *is*
tracked.

## Licensing

All five papers are open access under **CC-BY** (the foundational Wickert &
Schildgen, 2019, and its corrigendum are CC-BY 4.0; the GRL, ESurf, and AGU
Advances papers are likewise CC-BY). CC-BY permits reuse and adaptation with
attribution, so equations and paraphrased results may be reproduced in the GRLP
documentation as long as each is attributed to its source. Even so, these notes
paraphrase and cite rather than copy prose at length.

## The papers at a glance

| Role for GRLP | Citation | DOI |
|---|---|---|
| **Foundation** — defines the model | Wickert, A. D. & Schildgen, T. F. (2019), Long-profile evolution of transport-limited gravel-bed rivers, *Earth Surf. Dynam.*, **7**, 17–43 | [10.5194/esurf-7-17-2019](https://doi.org/10.5194/esurf-7-17-2019) |
| **Correction** — fixes the valley-width term | Wickert & Schildgen (2020), Corrigendum to the above | [10.5194/esurf-7-17-2019-corrigendum](https://doi.org/10.5194/esurf-7-17-2019-corrigendum) |
| **Linearization / spectral response** | McNab, F., Schildgen, T. F., Turowski, J. M. & Wickert, A. D. (2023), Diverse responses of alluvial rivers to periodic environmental change, *Geophys. Res. Lett.*, **50**, e2023GL103075 | [10.1029/2023GL103075](https://doi.org/10.1029/2023GL103075) |
| **Network extension** | McNab, F., Schildgen, T. F., Turowski, J. M. & Wickert, A. D. (2025), Influence of network geometry on long-term morphodynamics of alluvial rivers, *Earth Surf. Dynam.*, **13**, 1059–1092 | [10.5194/esurf-13-1059-2025](https://doi.org/10.5194/esurf-13-1059-2025) |
| **Application / feature test** | Ruby, A., McNab, F., Schildgen, T. F., Wickert, A. D. & Fernandes, V. M. (2026), How sediment supply, sea-level, and glacial isostatic oscillations drive alluvial river long-profile evolution and terrace formation, *AGU Advances*, **7**, e2025AV002035 | [10.1029/2025AV002035](https://doi.org/10.1029/2025AV002035) |

---

## 1. Governing equations — the canonical (corrected) forms

These are the equations a correct GRLP implementation solves. They are the
Wickert & Schildgen (2019) equations **as amended by the 2020 corrigendum**. The
single substantive correction (see §2) is that the *valley width* `B` enters only
as a local prefactor `1/B` — never through its downstream derivative `∂B/∂x`.
The forms below already incorporate that fix.

**Symbols** (consistent throughout these notes and the code):

| Symbol | Meaning | GRLP code | Typical/constant |
|---|---|---|---|
| `z` | valley-floor (bed) elevation | `z` | — |
| `x` | down-valley distance | `x` | — |
| `t` | time | — | — |
| `Q_s` | bed-load sediment discharge | — | — |
| `Q` (`Q_w`) | channel-forming water discharge | `Q` | — |
| `B` | valley width | `B` | — |
| `U` | uplift / subsidence (source–sink) | via `set_uplift_rate` | — |
| `λ_p` | sediment porosity | `lambda_p` | **0.35** |
| `I` | intermittency (fraction of time at bank-full flow) | `intermittency` | e.g. 0.01–1 |
| `𝕊` | sinuosity (channel length / valley length) | `sinuosity` | ≥ 1 |
| `k_Qs` | lumped transport coefficient | `k_Qs` | **≈ 0.041** |
| `ε` | Parker excess-stress ratio at bank-full | `epsilon` | 0.2 |
| `τ*_c` | critical Shields stress | `tau_star_c` | 0.0495 |
| `φ` | Wong & Parker MPM coefficient | `phi` | 3.97 |
| `D` | median grain size D₅₀ | — | (cancels; see below) |

**Bed-load transport** — near-threshold, self-formed (Parker, 1978)
equilibrium-width channel (Eq. 18):

```
Q_s = -(k_Qs · I / 𝕊^(7/6)) · Q · (∂z/∂x) · |∂z/∂x|^(1/6)
```

A defining property: for an equilibrium-width channel the grain-size dependence
cancels (`q_s ∝ D^(3/2)` while channel width `b ∝ D^(-3/2)`), so **no grain size
appears in `Q_s`**. The lumped coefficient (Eq. 19) is

```
k_Qs = k_qs · k_b ≈ 0.041      (for ε = 0.2)
```

computed in the code from first principles in `bedload_lumped_constants()`.

**Mass conservation (Exner), corrected** (corrigendum Eq. 1 / B2):

```
∂z/∂t = -1/(B(1 - λ_p)) · ∂Q_s/∂x + U
```

`B` is a prefactor, **not** inside the derivative.

**Combined long-profile evolution equation** — a nonlinear diffusion equation
(corrigendum Eq. 20 / D1):

```
∂z/∂t = (k_Qs · I) / (𝕊^(7/6) (1 - λ_p)) · |∂z/∂x|^(1/6)
        · [ (7/6) · (Q/B) · ∂²z/∂x²  +  (1/B) · (∂Q/∂x) · (∂z/∂x) ]  +  U
```

The `(7/6)(Q/B) ∂²z/∂x²` group is the diffusive term; the `(1/B)(∂Q/∂x)(∂z/∂x)`
group is the advective term from downstream discharge growth. (There is **no**
`∂B/∂x` term — that was the corrigendum's fix.)

**Boundary conditions:**

- *Upstream* — Neumann (sediment-flux). A prescribed input `Q_s` sets the
  boundary slope (corrigendum Eq. D4):
  ```
  S_0 = -sgn(Q_s) · 𝕊 · (Q_s / (k_Qs · I · Q))^(6/7)
  ```
  Higher sediment-to-water supply → steeper boundary slope (Lane's balance).
  In GRLP this is `set_Qs_input_upstream(Q_s_0)`; `set_S0(S0)` is the equivalent
  slope-first entry point.
- *Downstream* — Dirichlet on elevation (base level); `set_z_bl(z_bl)`. Base-level
  rise/fall is applied through this boundary; `set_x_bl` additionally moves the
  outlet horizontally (used for shelf-coupled base level, §5).

**Numerical method:** semi-implicit finite difference. The diffusive part is
solved directly as a **tridiagonal** system; the weak nonlinearity from
`|∂z/∂x|^(1/6)` and any nonlinear `Q(x)`/`B(x)` is resolved by **Picard
iteration** (`set_niter`). The corrigendum also restored a `𝕊^(7/6)` sinuosity
factor that had been omitted in the original numerical equations, and corrected
the sign/direction of the upstream Neumann condition. The GRLP release that
brought the code in line with the corrected equations was **v1.4.1**.

**Analytical steady-state profile** (corrected Eq. 40) — the recommended
verification benchmark, between two known endpoints `(x_0, z_0)` and `(x_1, z_1)`:

```
z = (z_1 - z_0) · [ (x^a - x_0^a) / (x_1^a - x_0^a) ] + z_0 ,   a = 1 - (6/7)P_xQ
```

The exponent depends on the distance→discharge power `P_xQ` **only** (the
corrigendum removed the valley-width exponent `P_xB` that was previously present).

**Slope–area / concavity** (corrected Eq. 54):

```
S = S_0 (A_0 / A)^((6/7) P_AQ) ,   concavity θ = (6/7) P_AQ ,   k_s = S_0 A_0^θ
```

Channel concavity is set by the drainage-area→discharge exponent `P_AQ` **alone**.

---

## 2. What the corrigendum corrected (2020)

**Root cause:** in the original paper the valley width `B` was erroneously kept
*inside* the spatial derivative of the Exner equation (by analogy to the 1-D
single-channel case). That incorrectly implied a river could aggrade or incise
purely because its valley widens or narrows downstream.

**The fix, applied everywhere:** remove all terms involving the downstream
derivative of valley width. `B` becomes a local `1/B` prefactor.

**Consequences that matter for the model and its documentation:**

- A spurious *mechanism* is deleted: downstream valley-width change no longer
  drives vertical channel change.
- **Section 5.3 of the original paper is invalidated** — valley widening no
  longer produces concave-up long profiles.
- Concavity is redefined: `θ = (6/7) P_AQ` (no valley-width exponent).
- The numerical scheme gained a previously omitted `𝕊^(7/6)` sinuosity factor,
  and the upstream Neumann boundary condition was corrected for sign/direction.
- Published figures changed only slightly in practice: the authors report Figs.
  2, 3, 5, 6, 8 "remain qualitatively — and nearly quantitatively — identical"
  after correction, and Fig. 4 is now valid only along the *x* axis.

**Implication:** anchor any implementation or documentation to the corrected
discretized equations (corrigendum Eqs. D2/D3), *not* to the original Eq. (1)/(20).
The continuous corrected Eq. (20) is easy to misread as typeset; the discretized
D2/D3 are the unambiguous, load-bearing forms.

---

## 3. Linearized response — gain and lag (McNab et al., 2023)

Takes the GRLP governing equation and, under small sinusoidal forcing of
**upstream** water and/or sediment supply, linearizes it to a constant-coefficient
diffusion equation for the elevation perturbation `δz`:

```
∂δz/∂t = κ · ∂²δz/∂x²
```

with the **sediment-transport diffusivity** (Eq. 5)

```
κ = (7/6) · (k_Qs · I · Q̄_w) / (𝕊^(7/6) · B (1 - λ_p)) · |∂z̄/∂x|^(1/6)
```

and the **characteristic equilibration time** (Paola et al., 1992)

```
T_eq = L² / κ .
```

`T_eq` is the central control. A forcing of period `P` defines a response length
`L' = √(P κ)`. The response is governed by the dimensionless ratio `P / T_eq`:

- **Fast forcing (`P ≪ T_eq`):** strongly damped (gain `G_z → 0`), large lag —
  terraces unlikely to form.
- **Slow forcing (`P ≫ T_eq`):** `G_z → 6/7` uniformly, response nearly in phase
  — aggradation/incision and terrace formation favored.
- Near `P ≈ T_eq`: transitional; e.g. a quarter-cycle lag (`φ_z/P ≈ 0.25`) near
  the outlet.

**Gain / lag** quantify amplitude ratio and phase shift of the response relative
to the forcing. Valley-floor elevation response is independent of whether water
or sediment supply is varied; sediment-discharge response differs (varying `Q_w`
can *amplify* the sediment signal, `G_Qs > 1`, and even make it lead the forcing).
For natural catchments `T_eq` spans ~10 kyr–1 Myr, overlapping Milankovitch
periods — so the same climate forcing yields diverse responses across catchments.

These are the functions in `grlp.py` — `compute_diffusivity`,
`compute_equilibration_time`, `compute_z_gain`, `compute_z_lag`,
`compute_Qs_gain`, `compute_Qs_lag`, and the supporting series-term methods — and
the theory behind `examples/McNab_et_al_GRL/`.

---

## 4. Network extension (McNab et al., 2025)

Extends GRLP from a single long profile to a **drainage network** of converging
segments, each obeying the same Eq. (3). This is the science behind GRLP v3's
`Network` class.

**Junction rules** (the mathematics added at confluences):

1. **Sediment conservation:** a segment's upstream sediment supply is the *sum*
   of the sediment discharges of its upstream segments.
2. **Elevation continuity:** `z` is continuous across every junction.
3. **Discrete discharge accumulation:** water discharge steps up at each
   tributary junction (plus optional continuous along-stream addition via the `U`
   source term).

Solved with the same semi-implicit finite-difference GRLP scheme, generalized to
the branching graph, with supply imposed at channel-head (inlet) segments and
fixed base level at the outlet segment.

**Network geometry:** topologies are generated as **Shreve (1966, 1974) random
binary trees** (all topologically distinct trees equally likely) — the basis of
GRLP's `generate_random_network` / `Shreve_Random_Network`. Downstream discharge
growth follows a Hack-related power law `Q_w = k · (x + x_0)^p`. Horton–Strahler
and Schumm ratios (bifurcation `R_B`, length `R_L`, area `R_A`, discharge `R_Qw`),
Tokunaga metrics, Strahler orders, and topological lengths are all computed — the
analysis methods now living on the `Network` class (`compute_strahler_orders`,
`compute_horton_ratios`, `compute_tokunaga_metrics`, `compute_topological_lengths`,
`find_hack_parameters`, etc.).

**Central result:** integrated network response resembles a single segment, but
the effective response time scales with the network's **mean inlet-to-outlet
length `⟨L_1⟩`**, *not* its longest (main-stem) length:

```
L̂ = √(T̂_eq ⟨κ⟩) ,   T_eq ≈ k_L ⟨L_1⟩² / ⟨κ⟩ ,   with  L̂ / ⟨L_1⟩ ≈ 1.35–1.45 .
```

Spatial patterns of aggradation/incision within a network are complex and
structure-specific; short-period forcing at the outlet is controlled by the local
network structure near the outlet, not by distant headwaters.

> Citation note: this paper cites the network numerical scheme as "Wickert et
> al., 2025," distinct from the GRLP v2.0.0 code archive (Zenodo
> [10.5281/zenodo.17091246](https://doi.org/10.5281/zenodo.17091246)). Whether
> that scheme reference is the same code release or a separate methods paper is
> not fully disambiguated in the text — confirm before finalizing a bibliography.

---

## 5. Application and feature test — Río Santa Cruz (Ruby et al., 2026)

Applies GRLP to the ~230 km, tributary-free **Río Santa Cruz** (Patagonia,
Argentina), which drops ~180 m from Lago Argentino to the Atlantic and preserves
≥13 alluvial terraces up to ~400 m above the modern river. It asks which
driver(s) — oscillating `Q_s/Q_w`, sediment pulses, glacial isostatic adjustment
(GIA), or sea-level/base-level change — produced the terraces, given a response
time far longer than glacial cycles.

**GRLP features this exercises / validates** (directly relevant to v3):

- **Shelf-coupled base level:** sea-level change moves the shoreline *along the
  continental-shelf gradient*, so the outlet's `x` and `z` move together. Vertical
  base-level change drives incision/aggradation *only where the exposed shelf
  slope differs from the channel slope* — a steeper exposed shelf during
  regression launches an upstream-propagating incision wave; transgression drives
  aggradation. This is the physics behind GRLP's `set_x_bl` + `set_z_bl`
  combination and the `examples/one_dimensional/RioSantaCruz_set_S0_set_x_bl.py`
  example.
- **Slope (flux) upstream boundary** with constant `Q_w` (`set_S0` /
  `set_Qs_input_upstream`).
- **Spatially/temporally varying source–sink `U`** carrying regional uplift, GIA
  flexure, and downstream gravel fining (Sternberg-type loss).
- **Diffusivity / equilibration-time machinery** (`κ`, `T_eq`, `L' = √(Pκ)`).

**Río Santa Cruz parameters used** (report as-published; some are unstated):

| Quantity | Value |
|---|---|
| Valley width `B` | 10 km (constant) |
| Length `L` | ~230 km |
| Water discharge `Q_w` | 1300 m³/s |
| Intermittency `I` | 0.5 |
| Sediment supply `Q_s,0` | 0.0127 m³/s (bedload, back-derived from slope) |
| Transport coefficient `k_Qs` | ≈ 0.041 |
| Regional uplift `U` | ~0.1 mm/yr (0.3 mm/yr in one GIA-offset run for visibility) |
| Downstream fining | 1.63×10⁻⁶ /m (0.163 %/km) |
| Equilibration time `T_eq` | ~1.2 Myr (≫ forcing periods → permanently transient) |
| Shelf slopes | steep proximal ≈ 0.0016; shallow ≈ 0.0003; global avg 0.0009 (Kennett, 1982) |
| Sea-level forcing | one 120 m / 100 kyr cycle; sinusoidal runs at 20 & 30 m amplitude, 41 & 100 kyr |
| `Q_s/Q_w` forcing | amplitudes 33 % & 66 %; periods 23, 41, 100 kyr |

*Equilibrium slope `S_0` is not printed as a single number* — the profile is
described as near-linear (~180 m over ~200–230 km). The GRLP example uses the
implied `S_0 = 180 m / 250 km = 7.2×10⁻⁴`.

**Key result:** direction of signal propagation is diagnostic — upstream forcings
(`Q_s/Q_w`, pulses) propagate *downstream*; base-level fall propagates *upstream*.
Terrace formation can lag peak forcing by up to ~half a cycle, and GIA can even
produce terraces *before* peak glaciation (negative lag). No single driver
explains the full terrace sequence; base level and/or GIA are needed for the
downstream terraces, `Q_s/Q_w` and/or GIA for the upstream ones. The paper is a
strong utility test of GRLP's ability to handle a **moving base-level boundary
coupled to a shelf gradient** together with spatially varying uplift.

---

## How these map to citing GRLP

Following the tiered guidance in the repository README:

- **Always:** Wickert & Schildgen (2019) + the software release (`CITATION.cff`).
- **Linearization / spectral-response features** (gain, lag, diffusivity): also
  McNab et al. (2023).
- **Network features:** also Wickert et al. (GRLP software) and McNab et al. (2025).
