# Next step: remove x_ext / A_ext (make LongProfile fully node-based)

## Status
De-pad + solver unification are **done, green, on `master`** (also merged and
force-pushed to origin/master; backups: remote `archive/remote-master-pre-rebuild`,
local `depad-session-backup`). Last committed step: `eb6c5d7` removed the
standalone **z_ext**. Full suite green (304 passed, 1 xfailed = Sternberg gravel
loss, deliberately broken). Local master is ~3 commits ahead of origin (CHANGELOG
+ z_ext + this note if committed) — **not pushed** yet; push is gated on Andy.

## Goal (Andy: "fully node-based; we'll thank ourselves later")
Remove the remaining static padded arrays from `LongProfile`: `x_ext`, `dx_ext`,
`dx_ext_2cell`, the dead `dx_ext_2cell__left/cent/right`, **and** `A_ext` (padded
drainage area). The walker never uses them; they're only a setup/diagnostic
convenience. `Q_ext`/`z_ext` are already gone.

## Design decision to make (Andy's hint: scalar / per-head-array-lookup / generate-in-place)
Two kinds of boundary-ghost quantity:
1. **Ghost node POSITIONS** (`x_ext[0]`, `x_ext[-1]`): trivially
   `2*x[0]-x[1]` and `2*x[-1]-x[-2]`. Recommend **generate in place** — that's
   exactly what `assemble_by_walking` already does, so it stays consistent and
   needs no storage.
2. **Analytic power-law ghost VALUES** (`Q_ext[0]` from `k_xQ*x_ext**P_xQ`;
   likewise `A_ext`): already handled for Q as **scalar attributes**
   `self.Q_ghost_upstream/downstream` (set in `set_Q`) — the walker reads them
   with a linear fallback. Follow the same scalar-attribute pattern for A if the
   `q_R`/`A_R` path needs a ghost. `z_ghost_upstream` is the existing precedent.
   (Per-head arrays are a Network concern; a single LongProfile has one head, so
   scalars suffice.)

So: **positions -> generate in place; power-law ghost values -> scalar attrs.**

## Sites (grep `x_ext`, `dx_ext`, `A_ext` in grlp/grlp.py)
- `__init__` (~25-28): drop `x_ext/dx_ext/dx_ext_2cell` defaults (and `A_ext`?).
- `set_x` (~106-178): 3 branches (x / x_ext / dx,nx,x0) each build
  `x_ext,dx_ext,dx_ext_2cell,dx_ext_2cell__{left,cent,right}`. Keep `x,dx,dx_2cell`.
  For the `x_ext=`-given branch keep deriving `self.x = x_ext[1:-1]`.
  `dx_ext_2cell__*` are DEAD (only set here; used by the deleted build_matrices).
- `self.L` (line 178): `x_ext[-1]-x_ext[0]` -> `(2*x[-1]-x[-2]) - (2*x[0]-x[1])`
  (preserves value for uniform dx; L is only used self-referentially in the
  test_api equilibration/wavenumber tests, so exact value is not pinned).
- `set_z` guard (~198): `self.x_ext.any()` -> `self.x` check.
- `set_A` (~204-218): `A_ext` created 3 ways; the `k_xA` branch uses `x_ext`.
  Remove `A_ext`; keep `self.A`. Migrate `set_Q(q_R, A_R)` (lines ~251-252, uses
  `self.A_ext`) and the test `tests/test_setup.py:99` (`A_ext[1:-1]==A`).
- `set_Q` (~254-256): `Q_ext = k_xQ*x_ext**P_xQ`. Compute the two ghost Q values
  from the ghost positions inline: `k_xQ*(2*x[0]-x[1])**P_xQ`, etc., store as
  `Q_ghost_upstream/downstream` (already the pattern). Interior `Q` already uses
  `self.x`.
- `set_B` (~301): only a spurious `self.x_ext.any()` guard; body uses `self.x`.
  Change the guard to `self.x`.
- `set_x_bl` (~383): drop `self.x_ext[-1] = self.x_bl` (keep `self.x_bl = x_bl`).
- `analytical_threshold_width_uplift` (~570-571): `x0=x_ext[0]`,`x_bl=x_ext[-1]`
  -> `2*x[0]-x[1]`, `2*x[-1]-x[-2]`.
- `compute_Q_s` / `slope_area` (~600-670, my z_ext-removal code): currently use
  `self.dx_ext[0]` and `self.dx_ext_2cell`. Change `dx_ext[0]` -> `self.dx[0]`
  (equal: `x[1]-x[0]`), and build the 2-cell dx from a locally reconstructed
  `x_ext = hstack((2*x[0]-x[1], x, 2*x[-1]-x[-2]))` (mirror Network.compute_Q_s,
  which already builds `_x` locally).
- `compute_length` (Network, ~677): `for x in self.x_ext` — network segments have
  no x_ext now (would break if called); rewrite to use `seg.x` (x.max()-x.min()).

## Plan (staged; keep suite green + commit at each stage)
1. Diagnostics + analytical off x_ext/dx_ext (compute_Q_s, slope_area,
   analytical) — reconstruct ghosts inline. Test. Commit.
2. Setters: remove A_ext (set_A) and the x_ext-based ghost in set_Q; fix set_B
   guard; migrate set_Q(q_R) + the A_ext test. Test. Commit.
3. set_x / L / set_x_bl / set_z guard / __init__ / compute_length — remove
   x_ext/dx_ext creation. Test. Commit.

## Guards / gotchas
- Characterization goldens use `set_x(dx,nx,x0)`, where `x_ext[0]=x0-dx=2*x[0]-x[1]`,
  so generate-in-place keeps compute_Q_s goldens **bit-for-bit**. Verify.
- `dx_ext[0] == dx[0] == x[1]-x[0]` (uniform-ghost identity).
- `set_x(x=...)` currently builds `x_ext=[nan,x,nan]` (nan ghosts) — the x-only
  path already yields nan boundary diagnostics; generate-in-place is an
  improvement there.
- Examples using `set_Q(q_R,A_R)` / reading `x_ext`/`A_ext` are deferred with #25.
- Full suite ~2.3 min; run in background to a file.

Working method that mattered this whole effort: read + reuse Andy's functions
before writing; preserve his comments; investigate before proposing.
