# Design note: second-order-in-time (BDF2, + adaptive Δt) for GRLP

**Status:** DESIGN / scope for review — not a commitment, no code yet.
**Motivates:** [MNiMORPH/GRLP#16](https://github.com/MNiMORPH/GRLP/issues/16).
**Benchmark:** `tests/test_time_accuracy.py`, `docs/accuracy.md` (current scheme is
first-order in time; final-profile error ∝ Δt).
**Template:** WTM's `benchmark/BDF2_ADAPTIVE_DESIGN.md` (recent, C++/PETSc) — same
family of problem (semi-implicit, Picard-frozen, dissipative), and its hard-won
lessons are inherited below.

---

## 1. Goal

Make GRLP's transient runs second-order in time so a target path accuracy is
reached with far fewer, larger steps (final error ∝ Δt² not Δt → for a fixed
tolerance, Δt roughly `1/√tol` larger). Keep the current backward-Euler scheme as
the default; gate the new path behind a flag. Then add adaptive Δt to size the
step to a user tolerance automatically.

## 2. Method: BDF2, not Crank–Nicolson

BDF2 fits a quadratic through the last three time levels:

$$\frac{3 z^{n+1} - 4 z^{n} + z^{n-1}}{2\,\Delta t} \;=\; \mathcal{L}(z^{n+1}),$$

local error O(Δt³), global O(Δt²). We choose it over Crank–Nicolson for the same
reason WTM did: **CN is A-stable but not L-stable** — its stiff-mode amplification
→ −1, so it *rings* on sharp fronts (e.g. an incision wave from a base-level drop).
**BDF2 is L-stable**, damping stiff modes toward zero. For a stiff diffusive
problem this is the correct 2nd-order method.

## 3. How it maps onto `grlp/solver.py` (small, local change)

The current step is backward Euler, `(z^{n+1} − z^{n})/Δt = 𝒟 z^{n+1} + src`,
assembled as `(I − Δt 𝒟) z^{n+1} = z^{n} + Δt·src` inside the Picard loop
(`assemble()` builds the stencil; `C0 ∝ Δt`; the RHS uses `seg.zold` = zⁿ). BDF2
rearranges to

$$\left(\tfrac{3}{2} I - \Delta t\,\mathcal{D}\right) z^{n+1}
   \;=\; \tfrac{1}{2}\,(4 z^{n} - z^{n-1}) + \Delta t\cdot src,$$

so only **two** things change, exactly as in the WTM note:

1. **Time diagonal** scaled `1 → 3/2` (the identity term in the stencil; `Δt·𝒟`
   and `Δt·src` are unchanged).
2. **RHS history** `z^{n} → (4 z^{n} − z^{n-1})/2` (replace the single `seg.zold`
   term with the two-level combination).

Implementation touches: store a second history level `seg.zold2`; carry a
`time_order` / θ-like switch into `assemble()` (and the junction rows) so the
diagonal identity coefficient and the RHS history term take BDF1 or BDF2 values;
the Picard freezing of the slope-dependent flux `𝒟` is untouched. **The SPD /
sparse-solve structure and unconditional stability are preserved.**

## 4. Why GRLP should get *clean* 2nd order (the key contrast with WTM)

WTM's achieved order collapsed to ~1 at fine Δt. It is worth being precise about
*why*, because it is **not** the intuitive answer. WTM's constitutive functions
are kinked (piecewise-C⁰ Fan `T`, `S`), and the natural hypothesis was that those
kinks cap the order — but that was **tested and disproved**: smoothing `T` to C∞
and widening the `S` transition band 100× left the order at ~0.9 and the errors
unchanged. The real cause was a **2-level secant effective storativity**
`(V(hⁿ⁺¹)−V(hⁿ))/Δh` — a *temporal-consistency* mismatch, its 2-level secant
driving BDF2's 3-level derivative — capping the order regardless of smoothness.
The proof: with **constant** storativity (`S ≡ porosity`, i.e. *linear* storage)
BDF2 recovered order ~2. The fix was "BDF2-on-V" (apply the 3-level difference to
the stored volume, not `S_eff·BDF2(h)`).

**GRLP does not have this trap (for now) — it is exactly WTM's constant-`S`,
order-recovering case.** GRLP's storage term is `∂z/∂t` with prefactor
`(1−λ_p)·B`; with **B constant in z**, the stored sediment is *linear* in z, so
the "secant" `ΔV/Δz` *is* the exact constant derivative — no 2-level/3-level
mismatch, `BDF2-on-z` *is* `BDF2-on-volume`, no order loss. GRLP's only
nonlinearity is the slope-dependent flux `𝒟` (a *coefficient*), Picard-frozen —
and WTM showed coefficient freezing does **not** cap the order. **Measured
(single-profile uplift benchmark): BDF2 achieves rate ~2.10, backward Euler
~1.03 — clean second order, as predicted.** (Kinks remain a *separate*,
milder flux-side risk at `S≈0`; see §8.)

> **Caveat, tied to #32 — and a design constraint to bake in there.** When valley
> width becomes dynamic (`B = B(z)` or `B(t)`), the stored sediment
> `V = (1−λ_p)∫B dz` is nonlinear in z and GRLP *would* inherit the WTM trap:
> writing #32 the natural `z`-first way (time term `≈ S_eff·(z−zold)/dt`,
> `S_eff = ΔV/Δz`) is the 2-level secant that caps BDF2 at first order and risks
> mass-conservation error. **The fix is a formulation choice made at #32's outset,
> not a later patch: conserve the sediment *volume* `V`, derive `z` from valley
> geometry** — then the time term is BDF2-on-`V`
> (`(3Vⁿ⁺¹−4Vⁿ+Vⁿ⁻¹)/(2Δt) = −∂Q_s/∂x + src`), conservation is exact, and BDF2
> stays second order. `V` is exactly the sediment budget's channel-storage integral
> `∫(1−λ_p)·B·z dx`, so this aligns with `compute_sediment_budget`. Much cheaper
> baked in than discovered empirically later (as WTM did).
>
> **DONE — implemented as issue #18 (2026-07).** The solver is now volume-first:
> `assemble()` row-scales the elevation system by `J = dV/dz` and carries the
> history in `V`-space, with a Jacobian-based linearization correction (not a
> secant). The geometry is a pluggable `Segment.valley_width(z)` primitive with
> `storage_volume`/`storage_jacobian` (rectangular default = constant `B`, exact
> reproduction to machine precision, 328 golden-master tests green). A synthetic
> nonlinear map confirms BDF2 holds 2nd order (2.10) — and, with the correction
> disabled, collapses to ~0.3 (the secant, live). So #32's dynamic width plugs
> into a solver that already conserves volume and keeps BDF2 second-order; it only
> needs to override `valley_width`/`storage_volume` with the real geometry.

## 5. Bootstrap and history validity

- **Step 0** has no `z^{n-1}` → take one backward-Euler step to seed the history
  (standard).
- **Discontinuous forcing/BC changes** (a step in uplift, a base-level jump)
  break the smooth history BDF2 assumes → **re-bootstrap** one BE step at the
  discontinuity. (GRLP forcing is usually smooth, so this is rare; but the
  accuracy benchmark's uplift *switch-on* is exactly such a point — the figure's
  left panel deliberately measures away from it.)

## 6. Adaptive Δt (phase 2)

Size Δt to a user **path tolerance** without a ground truth, by comparing two
orders at the same step:

- **Embedded estimate (preferred):** BDF2 carries a BDF1/predictor companion that
  shares the solve; `‖z_BDF2 − z_predictor‖` estimates the local error (~free).
- **Step-doubling (fallback):** one Δt step vs two Δt/2 steps (~3 solves/step).
- **Controller:** accept if `err < tol`; `Δt_new = Δt·(tol/err)^{1/(p+1)}` (PI
  controller for smoothness), clamped growth ratio; reject+retry if `err > tol`.
- **Why local control ⇒ global accuracy here:** the problem is **dissipative**
  (diffusion damps old truncation error), exactly as WTM argued — so local-error
  control is trustworthy (unlike hyperbolic/chaotic systems). GRLP shares this.
- **Variable-step BDF2** must use the non-uniform 3-level coefficients (they
  change with the step ratio ω = Δtₙ/Δtₙ₋₁) or the order silently drops to 1; cap
  the growth ratio.

## 7. Backward compatibility, tests, and the golden masters

> **SUPERSEDED (2026-07):** BDF2 + iterate-to-convergence is now the **default**
> (see CHANGELOG and `docs/accuracy.md`). This section planned BDF2 as opt-in with
> backward Euler default; that flipped. Golden-master sets are now kept for *both*
> schemes rather than pinning only the backward-Euler one.

- **Backward Euler stays the default.** BDF2/adaptive are opt-in (e.g.
  `Network`/`LongProfile` `set_time_integration("BDF2")` or a solver flag), so
  every existing result and the **golden-master characterization tests are
  unchanged** (they pin the default BE scheme). No silent re-baselining.
- **Acceptance:** `tests/test_time_accuracy.py` already measures the convergence
  rate and asserts ~1 for BE; add a BDF2 parametrization asserting ~2 on the same
  single-profile uplift benchmark (and, later, a small network). For adaptive:
  assert the achieved path error tracks the requested tolerance.

## 8. Risks / open questions (check before/while coding)

- **Non-smoothness in the flux capping the order.** The transport `∝ |S|^{7/6}`
  has a mild kink where slope `S → 0` (a locally flat / aggrading node). WTM found
  *coefficient* non-smoothness did **not** cap order — but confirm on GRLP by the
  self-convergence check; if a run spends time near `S≈0`, watch the achieved
  rate.
- **BDF2's order is contingent on Picard convergence.** The slope-flux coefficient
  is Picard-frozen and only becomes the fully-implicit value at `z^{n+1}` once the
  iteration converges; an *under-converged* solve reintroduces a coefficient lag
  that caps the order. Measured: BDF2 rate is 1.82 at `niter=1` but 2.10 at
  `niter≥2` (Picard converges superlinearly here, so 2–3 iterations suffice for
  the current weak `S^{1/6}` nonlinearity). Stiffer physics (dynamic `B`, #32) or
  very large steps could need more — which is the standing argument for
  **iterating Picard to a tolerance instead of a fixed count**
  ([#17](https://github.com/MNiMORPH/GRLP/issues/17)); that enhancement would make
  BDF2's second order unconditionally robust. The accuracy test guards this (the
  order drops if the solve is under-converged).
- **Junction rows.** The BDF2 diagonal/RHS change must be applied consistently in
  the multi-tributary junction cell (`_face_conductance` path), not just interior
  nodes.
- **Networks.** Verify order-2 on a small network, not only a single profile.

## 9. Phases (bounded, each independently verifiable)

1. **BDF2, fixed step**, behind a flag; `zold2` history + BE bootstrap; assert
   order-2 on the benchmark. (Small, self-contained — the highest-value slice.)
2. **Adaptive Δt** (embedded estimate + PI controller), behind a flag; assert the
   path error tracks tolerance; variable-step BDF2 coefficients.
3. (Optional) network + non-smooth-slope stress cases; docs update
   (`docs/accuracy.md` right panel bends 1 → 2).

Keep BE the default throughout; land phase 1 before phase 2.
