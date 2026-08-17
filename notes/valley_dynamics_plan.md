# Valley-dynamics implementation plan

Implementation map for the transient valley-width prototype (branch
`valley-realism`). Physics is settled in `notes/valley_width.md`; the confined
transient widening follows Turowski et al. (2025). All valley-process rules and
state live on `Segment` (locality principle); the solver adds one orchestration
line; `Network` is untouched.

## Settled equations this implements

- **Widening** (deterministic Turowski 2025, Eqs. 9, 12, 17): exact-exponential
  update over a step,
  `B <- W0 - (W0 - B)*exp(-dt/tau)`, `tau = (W0 - b)*(1 - h/H_valley)/zetadot`,
  `W0 = k0*h + b`. Asymptotes to `W0`, no overshoot, uncapped. Where the walls
  are no taller than the channel (`H_valley <= h`) the confined factor is
  undefined, so fall back to the unconfined `tau = (W0 - b)/zetadot` and warn
  loudly. Corrected wall factor is `(1 - h/H_valley)` (paper Eq. 17 prints a
  flipped `(1 - H_W/h)`; typo confirmed against Eqs. 15-18, 23).
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

Already built (committed): state `B`, `H_valley`, `f_ch`, `zetadot`, `k0`;
setters `set_lateral_migration_rate`, `set_valley_wall_height`,
`set_channel_belt_coefficient`, `set_valley_dynamics`.

| method | status | role |
|---|---|---|
| `update_valley(dt)` | rework | between-step entry point the solver calls; no-op guard, then dispatches the three rules |
| `_narrow_by_incision(dt)` | new (from inline) | incision rule; gated by `narrow_by_incision` |
| `_widen_valley(dt)` | new (rework old block) | exact-exponential widening; gated by `zetadot != 0` and `k0` set |
| `_partition_by_aggradation(dt)` | new | deposit fraction `f_ch = B_c/B`; gated by `partition_by_aggradation` |
| `channel_belt_width()` | new | `W0 = k0*h + b` (Turowski Eq. 4; Fergus: `compute_B_0`) |
| `storage_volume(z)`, `storage_jacobian(z)` | rework | fold in `f_ch`: `(1-lambda_p)*f_ch*B*z` and `(1-lambda_p)*f_ch*B`. Guard: 328 golden-master tests stay green at `f_ch = 1` |
| `lateral_migration_rate_from_Qs()` | deferred | `zetadot_mig` from `Q_s` (Fergus: `compute_q_L`) |

## `solver.py`

One line, twice: `for seg in net.segments: seg.update_valley(dt)` right after
`dz_dt` is published -- in `evolve` (~:420) and `evolve_adaptive` (~:521).

## `Network`

Nothing. It already produces `seg.Q_s` (grlp.py:1091); all valley rules are local.

## Build order

`channel_belt_width` -> `_widen_valley` (rework) -> `_narrow_by_incision`
(refactor) -> `_partition_by_aggradation` -> wire `update_valley` -> storage fold
(run golden-master) -> solver orchestration -> tests. Granular commit per step.

## To reconcile with Fergus (co-author of Turowski et al. 2025)

1. **Where widening happens.** This design updates the valley *between* implicit
   steps (state frozen within the solve -> storage stays linear in `z`, BDF2
   clean). Fergus's `notes/incision.md` proposes solving `B^{n+1/2}` *inside* the
   Picard iteration. Real architectural difference, not cosmetic.
2. **Naming.** Private rule methods here match the flag names; Fergus used public
   `compute_*`. Can switch the shared helpers (`channel_belt_width -> compute_B_0`)
   to his convention.
3. **Eq. 17 typo** in the 2025 paper (`(1 - H_W/h)` should be `(1 - h/H_W)`) --
   flag to the authors.
