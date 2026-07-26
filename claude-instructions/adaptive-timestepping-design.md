# Design note: adaptive time stepping for GRLP (#16 Phase 2)

**Status:** DESIGN / scope for review — not started in code.
**Part of:** [MNiMORPH/GRLP#16](https://github.com/MNiMORPH/GRLP/issues/16) (part 3).
**Builds on (both DONE):** BDF2 second-order-in-time (#16 Phase 1) and the
volume-first solver (#18). Companion: `second-order-time-bdf2-design.md`.
**Template:** WTM `benchmark/BDF2_ADAPTIVE_DESIGN.md` §3.

## Goal

Size the time step automatically to hold a user **path tolerance** (e.g. "keep
the transient profile within X m"), taking small steps when the transient is
fast/nonlinear and growing them as it settles. Backward Euler / fixed-step BDF2
stay the defaults; adaptive is opt-in.

## Key architectural finding (from prototyping, 2026-07)

Adaptive Δt **cannot** be a wrapper around repeated `evolve(nt=1)`: the current
`solver.evolve` re-bootstraps the BDF2 history on the first step of *every* call
(`net._bdf2_step = bdf2_requested and ti > 0`), so each 1-step call is a
backward-Euler bootstrap (an embedded `‖z_BDF2 − z_BE‖` estimate came out exactly
0 for this reason). Adaptive stepping needs a **dedicated loop** that (a) carries
the `(zⁿ, zⁿ⁻¹)` history across steps, (b) estimates the per-step error, and
(c) resizes Δt each step.

## Plan (increments, each independently verifiable)

1. **Dedicated `evolve_adaptive` loop.** Carry `zold`/`zold2` across steps;
   bootstrap once at the start (and re-bootstrap only at genuine discontinuities,
   e.g. a step change in a boundary condition). Fixed Δt first, to confirm it
   reproduces `evolve` BDF2 (order ~2) before adding control.
2. **Local-error estimate — step-doubling (chosen).** Andy's call (2026-07),
   after prototyping both. The *embedded* `‖z_BDF2 − z_BE‖` estimate (1 extra
   solve) is cheap and always conservative, BUT it estimates the backward-Euler
   (O(Δt²)) error, not BDF2's (O(Δt³)); measured, its ratio to the true BDF2
   error **drifts from ~5 to ~140** as the solution smooths (it is a constant-in-Δt
   factor only at fixed order, and the two orders differ), so controlling on it
   throttles step growth ~10× in exactly the smooth regions adaptivity is for. So
   use **step-doubling**: one BDF2 step of Δt (`z_big`) vs two of Δt/2 (`z_small`),
   `est = ‖z_big − z_small‖`, advance with `z_small` (local extrapolation),
   ~3 solves/step (cheap here — the per-step sparse solve is a few-hundred nodes).
3. **PI step controller.** Accept if `est ≤ tol`; propose
   `Δt_new = Δt·(tol/est)^{1/(p+1)}` (p = 2), smoothed by a PI controller, clamped
   by a safety factor and max growth ratio; **reject and retry** with smaller Δt
   if `est > tol`. Never step past a forcing cadence if one is set.
4. **Variable-step BDF2 coefficients — DONE (prerequisite, pulled forward).**
   Was planned last; step-doubling's first half-step starts from a Δt-spaced
   history, so it *needs* non-uniform BDF2 before the estimator can exist. With
   ω = Δtₙ/Δtₙ₋₁ the 3-level weights are a=(1+2ω)/(1+ω), b=1+ω, c=ω²/(1+ω) (time
   term `[a zⁿ⁺¹ − b zⁿ + c zⁿ⁻¹]/Δtₙ`); ω=1 reduces to (3/2, 2, 1/2) bit-for-bit.
   Implemented in `assemble` (reads `net._bdf2_omega`, default 1); volume-first V
   history uses the same (b, c). Verified 2nd order on a random-ratio non-uniform
   sequence; uniform weights on it collapse to ~1.6.
5. **API + tests.** `LongProfile`/`Network.set_adaptive_timestep(tol)` (or an
   `evolve(..., tol=)` mode). Test: the achieved path error (vs a fine reference)
   tracks the requested tolerance across a range of tol; the step grows as the
   transient settles; results are Δt-history-independent to tolerance.

## Progress and empirical findings (2026-07, unpushed)

**Done and committed** (see solver.py / test_adaptive_timestepping.py):
- Increment 1: `evolve_adaptive` dedicated history-carrying loop; reproduces
  `evolve` BDF2 bit-for-bit at fixed Δt, order ~2.1.
- Increment 4 (variable-step BDF2 weights): as above.
- `_picard_step` / `_picard_config` extracted so `evolve` and `evolve_adaptive`
  share one step-taker.

**Increments 3 & 5 (controller + API) — DONE and committed.**
- `solver.evolve_adaptive(net, T, dt_init=None)` is now the adaptive integrator
  (advance total time T). `_trial_step` does the step-doubling; an I-controller
  `Δt·safety·(tol/est)^{1/(p+1)}` (p=2, p=1 bootstrap) with reject-and-retry sizes
  the step; advances with the finer half-step solution. Error-controlled BE
  bootstrap → result independent of `dt_init`.
- API: `set_adaptive_timestep(tol, dt_init, dt_min, dt_max, safety, max_grow,
  max_shrink)` on Network + LongProfile; `evolve_threshold_width_river_adaptive(T)`
  (LongProfile) / `..._network_adaptive` (Network).
- Verified (tests/test_adaptive_timestepping.py): achieved error falls
  monotonically with tol, dt_init-independent to ~1e-6 m over 1000×, lands exactly
  on T, denser stepping in the fast early transient, runs on a confluence network;
  golden masters untouched (340 passed).
- **Chose an I-controller, not PI** (design note said PI): the I-controller meets
  the acceptance test and avoids gain tuning; PI smoothing is a documented
  refinement, not a requirement. **Semantics:** `tol` is a *per-step local*
  tolerance (ODE-solver convention); achieved path error is comparable, not equal.

**Step-doubling estimator — characterized (prototype, informed the build):**
- The advanced solution's true error ≈ **0.83·est**, and this ratio is
  **constant across Δt** (measured 0.40 for `est/3`, flat over 500–8000 yr steps).
  Constant calibration is exactly what the controller needs — the estimator is
  well-behaved, unlike the embedded one. The factor is not the textbook 1/3
  because `z_small`'s first half-step is a non-uniform ω=0.5 BDF2 step (different
  error constant); it does not matter that it is 0.40 not 0.33, only that it is
  constant. **Controller should target `est ≤ tol`** (raw difference) → achieved
  error ≈ 0.83·tol, safe and ~1.2× conservative.
- **Prototyping gotcha (cost real time):** isolated local-error-vs-Δt scaling is
  *not* cleanly measurable — varying Δt either changes ω (if the warm-up spacing
  is held fixed → large-ω error blows up) or moves the start-state along the
  transient (if warm-up spacing = Δt → different C(t)). Both confound the scaling.
  Rely on the **end-to-end order-2 tests** for order; use the isolated runs only
  for the *calibration ratio* (which is robust to these confounds).

**Status: all increments (1, 3, 4, 5) DONE.** Adaptive time stepping is
implemented, tested, and documented (CHANGELOG + docs/accuracy.md). Possible
future refinements (not required): PI step controller (smoother step sequence);
per-unit-time (EPUS) tolerance so achieved *path* error ≈ tol directly;
re-bootstrap on boundary-condition discontinuities; expose adaptive on the public
`evolve_threshold_width_river(..., adaptive=True)` signature if a single entry
point is preferred over the `_adaptive` method.

## Risks / interactions

- **Variable-step BDF2 must use the correct non-uniform weights** — the single
  most likely place to silently lose the order. Guard with the self-convergence
  order check.
- **Picard convergence (#17) — DONE.** `set_picard_tolerance` landed; the
  controller should turn it on so the per-step solve it trusts is converged (a
  fixed `niter` under-converges at large adaptive steps and corrupts the estimate).
- **Volume-first (#18).** All of the above operates on V-space history; the
  estimator differences z (the observable), so keep the V↔z map consistent.
- **Bootstrap after discontinuities** (BC steps): re-bootstrap or the history is
  invalid.

## Recommendation

A dedicated loop + step-doubling estimator + PI controller + variable-step BDF2 is
a multi-increment effort; do it with the order-check after each increment. #17
(iterate-to-convergence) and the variable-step weights are **done**; the remaining
work is the estimator function, the controller, and the API. The one retired risk
worth restating: step-doubling's estimator is **well-calibrated (constant factor)**,
so a plain `est ≤ tol` controller is trustworthy.
