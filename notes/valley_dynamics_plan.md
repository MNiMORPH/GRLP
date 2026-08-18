# Valley-dynamics implementation plan

Implementation map for the transient valley-width prototype (branch
`valley-realism`). Physics is settled in `notes/valley_width.md`; the confined
transient widening follows Turowski et al. (2025). All valley-process rules and
state live on `Segment` (locality principle); the solver adds one orchestration
line; `Network` is untouched.

## Settled equations this implements

- **Widening** (deterministic Turowski 2025, Eqs. 9, 12, 17): exact-exponential
  update over a step,
  `B <- B_max - (B_max - B)*exp(-dt/tau)`, `tau = (B_max - b)*(1 - h/H_valley)/zetadot`.
  The maximum (unconfined) valley width `B_max` is the *prescribed* width from
  `set_B` (Fergus McNab's observed `B_br`), captured when B is set -- NOT a
  depth-scaled `k0*h + b` (that earlier form was unverified against the paper and
  contradicted by Fergus's code, which uses a measured constant; removed). Widening
  is opt-in via the `widen_by_migration` flag on `set_valley_dynamics`. Asymptotes
  to `B_max`, no overshoot, uncapped. Where the walls are no taller than the
  channel (`H_valley <= h`) the confined factor is undefined, so fall back to the
  unconfined `tau = (B_max - b)/zetadot` and warn loudly. Corrected wall factor is
  `(1 - h/H_valley)` (paper Eq. 17 prints a flipped `(1 - H_W/h)`; typo confirmed
  against Eqs. 15-18, 23).
- **Incision** (vertical walls): where `dz/dt < 0`, `B -> b` and
  `H_valley += |dz/dt|*dt`.
- **Aggradation deposit partition** (Wickert 2013 Eqs. 23, 31; `a_R=1, p_R=0`):
  where `dz/dt > 0`, `f_ch = B_c/B`,
  `B_c = b + (B - b)*(1 - exp(-zetadot*h / ((B - b)*etadot)))`,
  `etadot = dz/dt` (single uniform rate for channel and floodplain, for now);
  `f_ch = 1` elsewhere.
- **Migration rate from discharge** (deferred; prescribe `zetadot` for now):
  `zetadot_mig = Xi/(1 - lambda_p) * Q_s/(h*dx)`, `dx = 1 m` (Wickert 2013
  Eq. 29 with the mean channel width replaced by a unit downstream length).

## `Segment` (grlp/grlp.py)

Already built (committed): state `B`, `H_valley`, `f_ch`, `zetadot`, `B_max`;
setters `set_lateral_migration_rate`, `set_valley_wall_height`,
`set_valley_dynamics`. `B_max` (prescribed max valley width) is captured from
`set_B`, not a separate setter.

| method | status | role |
|---|---|---|
| `update_valley(dt)` | rework | between-step entry point the solver calls; no-op guard, then dispatches the three rules |
| `_narrow_by_incision(dt)` | new (from inline) | incision rule; gated by `narrow_by_incision` |
| `_widen_valley(dt)` | new (rework old block) | exact-exponential widening toward `B_max`; gated by `zetadot != 0` and `widen_by_migration` |
| `_partition_by_aggradation(dt)` | new | deposit fraction `f_ch = B_c/B`; gated by `partition_by_aggradation` |
| `storage_volume(z)`, `storage_jacobian(z)` | rework | fold in `f_ch`: `(1-lambda_p)*f_ch*B*z` and `(1-lambda_p)*f_ch*B`. Guard: 328 golden-master tests stay green at `f_ch = 1` |
| `lateral_migration_rate_from_Qs()` | deferred | `zetadot_mig` from `Q_s` (Fergus: `compute_q_L`) |

## `solver.py`

One line, twice: `for seg in net.segments: seg.update_valley(dt)` right after
`dz_dt` is published -- in `evolve` (~:420) and `evolve_adaptive` (~:521).

## `Network`

Nothing. It already produces `seg.Q_s` (grlp.py:1091); all valley rules are local.

## Build order

`_widen_valley` (rework) -> `_narrow_by_incision` (refactor) ->
`_partition_by_aggradation` -> wire `update_valley` -> storage fold
(run golden-master) -> solver orchestration -> tests. Granular commit per step.

## To reconcile with Fergus (co-author of Turowski et al. 2025)

1. **Where widening happens -- now a user option** (`set_valley_coupling`).
   `'between_step'` (default) updates the valley once per step after the solve:
   storage stays linear in `z`, BDF2 keeps second order, geometry lags one step.
   `'in_picard'` recomputes the valley from each Picard iterate inside the solve
   (Fergus's `notes/incision.md` direction): tighter coupling, but storage is
   nonlinear in `z`, so BDF2 is no longer strictly second order (a warning fires;
   backward Euler recommended). `update_valley` is idempotent (resets to the
   frozen `Bold`/`Hold` each call) so it can run once or every iteration.
   In-Picard is `evolve`-only for now (not the adaptive solver). **Resolved:**
   keep the in-Picard *end-point* form -- the geometry converges to `B(z^{n+1})`
   (backward-Euler in the geometry: L-stable, monotone) -- rather than Fergus's
   midpoint `B^{n+1/2}` (Crank-Nicolson-like). The midpoint is formally second
   order in the geometry but buys nothing (the time integration is first order
   for valley dynamics anyway, since the state-dependent storage rules out BDF2)
   and, being non-L-stable, can ring on the sharp valley features (`B -> b`
   narrowing; the f_ch feedback). Tested both ways on aggradation and
   incision+narrowing: they agree to ~1e-5, so the end-point form is kept as the
   strictly more robust default. L-stability over formal order, per ADW.
2. **Naming.** Private rule methods here match the flag names; Fergus used public
   `compute_*`. Can switch the shared helpers to his convention if desired.
3. **Eq. 17 typo** in the 2025 paper (`(1 - H_W/h)` should be `(1 - h/H_W)`) --
   **agreed**: ADW verified independently and Fergus (a co-author) agrees, so it
   is confirmed; needs only a paper correction (his call as co-author). The code
   already uses the corrected `(1 - h/H_valley)` in `_widen_valley`.
