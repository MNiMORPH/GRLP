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
2. **Embedded local-error estimate.** At each step (past bootstrap) compute the
   BDF2 solution and the backward-Euler solution from the *same* Picard-frozen
   operator and history; `est = ‖z_BDF2 − z_BE‖` estimates the (lower-order) local
   error — ~one extra solve, not a re-integration. Valid because the problem is
   **dissipative** (diffusion damps old truncation error → local control ⇒ global
   accuracy; unlike hyperbolic/chaotic systems). Fallback: step-doubling (1 step
   of Δt vs 2 of Δt/2, ~3 solves/step).
3. **PI step controller.** Accept if `est < tol`; propose
   `Δt_new = Δt·(tol/est)^{1/(p+1)}` (p = 2), smoothed by a PI controller, clamped
   by a safety factor and max growth ratio; **reject and retry** with smaller Δt
   if `est > tol`. Never step past a forcing cadence if one is set.
4. **Variable-step BDF2 coefficients.** When consecutive steps differ
   (ω = Δtₙ/Δtₙ₋₁ ≠ 1) the 3-level weights are no longer (3/2, −2, 1/2); they must
   use the **non-uniform BDF2 formula** or the order silently drops to 1. This
   touches the volume-first time term: the diagonal factor and the V-space history
   combination `(… Vⁿ … Vⁿ⁻¹ …)` become ω-dependent. (A cruder first cut:
   piecewise-constant Δt in blocks — uniform BDF2 within a block, re-bootstrap on
   change — avoids the ω math but wastes a BE step per change; prefer proper
   variable-step.)
5. **API + tests.** `LongProfile`/`Network.set_adaptive_timestep(tol)` (or an
   `evolve(..., tol=)` mode). Test: the achieved path error (vs a fine reference)
   tracks the requested tolerance across a range of tol; the step grows as the
   transient settles; results are Δt-history-independent to tolerance.

## Risks / interactions

- **Variable-step BDF2 must use the correct non-uniform weights** — the single
  most likely place to silently lose the order. Guard with the self-convergence
  order check.
- **Picard convergence (#17).** The error estimate and the step assume a converged
  nonlinear solve each step; a fixed `niter` that under-converges at a large
  adaptive step corrupts both. Iterate-to-convergence (#17) makes this robust —
  arguably a prerequisite for trustworthy adaptivity.
- **Volume-first (#18).** All of the above operates on V-space history; the
  estimator differences z (the observable), so keep the V↔z map consistent.
- **Bootstrap after discontinuities** (BC steps): re-bootstrap or the history is
  invalid.

## Recommendation

A dedicated loop + embedded estimator + PI controller + variable-step BDF2 is a
multi-increment effort; do it with the order-check after each increment. Consider
landing #17 (iterate-to-convergence) first, so the per-step solve the controller
trusts is unconditionally converged.
