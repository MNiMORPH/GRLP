# v3 example / template / library curation — plan of record

Goal: for the **v3.0.0** release, turn the current heap of example scripts into a
clean, limited, *runnable* set — and lift reusable machinery up into the library
where it belongs. This doc is the durable plan (survives compaction); edit freely.

## The four-way target (where code lives, and why)

Each script gets sorted in **two moves**: first *what machinery do I lift into
`grlp/`*, then *is the remaining driver an example or a template*.

1. **`grlp/`** — reusable **machinery** currently trapped inside scripts, promoted
   to importable submodules. If people re-type it, it belongs in the library.
   Candidates: Fergus's `build_synthetic_network` / `Shreve_Random_Network`; the
   DEM smooth / despike / coarsen pipeline; valley-width extraction.
2. **`examples/`** — teach-first, minimal, **must run on the v3 API**, wired into
   the docs / RTD. The clean small subset.
3. **`templates/`** — scientific **template code**: a full real-study workflow
   meant to be *copied and modified*, not read start-to-finish for a concept.
   Heavier, carries real data, application-specific. (Andy's framing: "these are
   often really template code for scientific implementation.")
4. **`deprecated/`** — frozen, not maintained, left for reference. Already exists.

The one hard promise this buys us, and the release gate we can actually test:
**everything in `examples/` runs on the v3 API.** `templates/` is explicitly
"adapt me, may need data / porting."

## The v3-feature coverage gap (write these new)

Nothing anywhere — in-repo or across the machine — demonstrates the new v3 API:

- **`set_S0` demo** — force the upstream boundary by slope (DEM-derived slope input).
- **`set_x_bl` demo** — move the outlet laterally (base-level mobility); optionally
  combined `x_bl` + `z_bl` (mobile outlet + subsidence/uplift).

These are `examples/` (teach-first, synthetic, small).

## Cross-computer survey (external scripts — curate UP into the repo)

Swept `/home/awickert` for GRLP user scripts: 46 external scripts across
`Downloads/`, `Dropbox/GRLP/`, `Papers/Completed/`, `Papers/InProgress/`, and the
`Stanford2024/GRLP-vardx_mod/` copy, deduplicated against library copies / JGR
submission archives / SRLP (a different model). Full detail in the session
scratchpad reports (`FINAL_REPORT.txt`, `SUMMARY_TABLE.txt`, `classify.txt`);
distilled here.

### Pull in (highest value / novel)

- **Fergus River network** — `Dropbox/Papers/InProgress/FergusRiverNetwork/`
  (11 variants → 1 canonical + README). The **only** script using
  `Shreve_Random_Network()`; real DEM-derived network (Fergus River, Argentina;
  FluvMORPH FM20240807). **OLD API — needs v3 port.** Canonical variant:
  `basic_FM20240807.py` or `along_stream_sources_FM20240807.py`.
  → machinery (`Shreve_Random_Network`, DEM routing) to **`grlp/`**; the driver to
  **`templates/`**. Highest priority.
- **Yanshui width** — `Downloads/width_aw.py` (OLD API) +
  `Dropbox/Papers/Completed/LandslidesOlivia/GRLP_baselevel.py` (modern API).
  Real DEM + coupled valley-width dynamics + landslide-driven base-level rise.
  Ties into the #32 valley-realism roadmap. → width machinery to **`grlp/`**;
  driver(s) to **`templates/`**.

### Reference-tier (keep as reference; not featured examples)

- **Rio Santa Cruz** — `Stanford2024/.../RioSantaCruz2023*.py` (2 variants; modern
  API; real field). Already have an in-repo copy under `one_dimensional/`. Check
  variant `b` for genuinely new forcing before deciding feature-vs-reference.
- **Tofelde experiments** — `Papers/Completed/TerraceAggradationIncision_Experiments/`
  (2 scripts; modern; lab-experiment validation).

### Skip (redundant / archival)

MPT_tests (7), stepchange (3), oscillatingQ (covered by in-repo
`oscillatingQ_2panel`), Corrigendum / LongProfileGravelBedRivers ModelRuns
duplicates, FluvTree NetworkX (superseded by Fergus), `Downloads/width.py`
(dup of `width_aw.py`), `GRLP_experiment_tests` (dups of completed).

## In-repo triage — FIRST PASS (proposal; correct me)

Current tree (tracked): `one_dimensional/` (15), `network/` (26),
`network_debugging/` (8), `McNab_et_al_GRL/` (5), `deprecated/` (13), plus
top-level `run_grlp.py`, `example_1d.ipynb`. API column = touches removed API
(`z_ext`/`x_ext`/`dx_ext`/`A_ext`) → needs a port if kept living.

### → `examples/` (teach-first; must run on v3)
- `run_grlp.py` — CLI/basic runner
- `example_1d.ipynb` — Eric Barefoot's intro notebook
- `one_dimensional/Basic1D.py`
- `one_dimensional/analytical_numerical.py` (+ `_uplift`; `_irregular_grid` [OLD])
- `one_dimensional/analytical_concavity.py`
- `one_dimensional/base_level_fall_transience.py`
- `one_dimensional/sudden_bl_fall_U0.py`
- `network/NewNetwork_1segment.py` / `_2segments` / `_3segments` [OLD] — network basics
- **NEW:** `set_S0` demo, `set_x_bl` demo (see gap section)

### → `templates/` (scientific; copy-and-adapt; may carry data)
- `one_dimensional/RioSantaCruz2023.py` (+ `b`) — real field, Argentina
- `one_dimensional/transient_figure_with_SA_10panel_*.py` [OLD] — paper figure
- `network/whitewater_2m_net*.py` (5) + `whitewater-*.json` — real DEM, Whitewater
- `network/perrot_net.py` / `perrot_newnet.py` + `Perrot-Network-Graph.svg` [OLD]
- `network/landslide5.py` / `landslideContinuous5` / `landslideContinuousUpper5` [OLD]
- `network/netBLfall*` / `netBLrise*` / `netLandslide*` [OLD] — network transients
- `McNab_et_al_GRL/*` (5) [OLD] — published-paper reproduction (McNab et al. 2025)
- (pull-ins) Fergus River driver, Yanshui width driver

### → `grlp/` (lift machinery into the library)
- `network/smooth_despike_coarsen_network*.py`, `smooth_despike_network.py` —
  DEM smooth/despike/coarsen pipeline → a preprocessing submodule
- `network/monotone.py` / `monotone2.py` — monotonic-profile enforcement (confirm
  this is machinery, not a one-off driver)
- (pull-ins) `Shreve_Random_Network` / DEM-network routing from Fergus; valley-
  width extraction from Yanshui

### → `deprecated/` or drop from repo
- `one_dimensional/Alluvial_Neumann_2018_2.py` [OLD], `Basic1D-var-dx.py` [OLD],
  `oscillatingQ_2panel*.py` [OLD], `test_uplift_equals_qs.py` [OLD] (test-like —
  candidate to fold into `tests/` instead)
- `network/Test3segs.py`, `TestQ.py` — dev checks
- `network_debugging/*` (8) — dev/debug scratch; keep as a dev-only dir or deprecate
- `deprecated/*` (13) — stays

## VERIFIED — `one_dimensional/` (15 scripts; ground-truthed against grlp.py)

Key correction: **`set_x(x_ext=…)` and `set_z(z_ext=…)` are SUPPORTED v3 API**
(the intended way to set explicit ghost positions — grlp.py `set_x` lines 137–149).
Only reading/writing `lp.z_ext` / `lp.A_ext` as **attributes** is broken. So the
port list is smaller than a raw token grep suggests.

Needs port (4) — all via removed *attribute* access:
- `test_uplift_equals_qs.py` — `lp.z_ext`, `lp.A_ext` (also test-like → consider `tests/`)
- `transient_figure_with_SA_10panel_2020_GRLP2update2025.py` — `lp.z_ext` ×5
- `oscillatingQ_2panel.py` — `lp.z_ext` slope read
- `oscillatingQ_2panel_outrange.py` — `lp.z_ext` slope read

v3-clean (no port needed):
- `analytical_numerical.py`, `_uplift`, `_irregular_grid` (uses supported
  `set_x(x_ext=…)`), `analytical_concavity.py`
- `Basic1D.py`, `Basic1D-var-dx.py` (uses supported `set_x(x_ext=…)`)
- `base_level_fall_transience.py`, `sudden_bl_fall_U0.py`
- `Alluvial_Neumann_2018_2.py` — token `z_ext` is a LOCAL var in its own
  finite-difference solver, NOT GRLP API (grep false-positive)
- `RioSantaCruz2023.py`

Drop: `RioSantaCruz2023b.py` — identical to `RioSantaCruz2023.py` but for one
blank line (line 67). Duplicate.

**DONE — 1D sort complete** (commits `5017db7`, `886e213`). Final placement:
- **`examples/one_dimensional/`** (10, all verified to run on v3):
  `analytical_numerical` + `_uplift` + `_irregular_grid`, `analytical_concavity`,
  `Basic1D`, `Basic1D-var-dx`, `base_level_fall_transience`, `sudden_bl_fall_U0`,
  `oscillatingQ_2panel`, `test_uplift_equals_qs` (kept as a validation demo).
- **`templates/one_dimensional/`** (4): `RioSantaCruz2023`,
  `transient_…10panel`, `Alluvial_Neumann_2018_2` (own FD solver),
  `oscillatingQ_2panel_outrange` (period-sweep animation).
- **Dropped**: `RioSantaCruz2023b` (exact duplicate).

`z_ext` fully removed from the library first: single-segment diagnostics
delegate to the walk (no padded array anywhere). All 1D examples ported off
`z_ext`/`A_ext`.

## VERIFIED — `network/` machinery candidates (SET B)

Confirmed **reusable machinery, not one-off drivers** — every script defines
generalized, parameterized, documented functions and then drives them under a
main block. Lift into the library:
- smooth/despike/coarsen family (`despike_1d`, `smooth_1d_fixed_ends`,
  `clean_edge_values`, `downsample_edge_along_s_new`) → e.g. `grlp/` smoothing
  submodule. The three files (`smooth_despike_network`, `_coarsen_network`,
  `_coarsen_network_3DEP`) are the **same functions** differing only by input
  JSON and whether the coarsen step runs → consolidate to ONE parameterized
  module + a thin example driver.
- monotone family (`monotone_envelope_edge`, `symmetric_running_mean_all_edges`,
  `super_smooth_edge_monotone`, …) → same or sibling submodule. `monotone.py`
  and `monotone2.py` are two generations (2 = 3-pass running-min-envelope) →
  keep the better one's machinery, consolidate.
All v3-clean (no removed API). Each loads a `grass-whitewater-*.json` network.

## Where we left off / next

- **`one_dimensional/` DONE** (sorted, ported, gated). `z_ext` removed library-wide.
- **`network/` sort + example porting DONE.** Placement:
  - `examples/network/` (8, all run on v3): `NewNetwork_1/2/3/5segments`,
    `netBLfall`, `landslide5`, `landslideContinuous5`, `netBLrise_animation` (new:
    FuncAnimation → gif, replaced the frame-dumper). Ported off removed API (shared
    outlet-connector plot `lp.x_ext[0][-2:]` → `[x[-1],x_ghost_downstream],[z[-1],z_bl]`;
    base level via `set_z_bl`/`set_x_bl`; landslide pulse `lp.z[2] += 30`).
  - `templates/network/`: 5 `whitewater_2m_net*` + 3 whitewater JSONs; `perrot_newnet`
    (+ Perrot SVG; needs external grass-perrot.json).
  - `examples/network_debugging/`: Test3segs, TestQ.
  - `examples/deprecated/network/`: landslideContinuousUpper5, perrot_net,
    NewNetwork_5segments_2-3-reverse, netBLfall_movieframes, netLandslide_movieframes.
  - **FLAG:** `landslideContinuous5` in-loop perturbation hits the LAST segment
    (loop-variable leak), not segment 2 where the initial pulse lands. Preserved
    as-is; ask Andy whether to fix to segment 2 for teaching clarity.
- **Next in `network/`: machinery lift** — `monotone`/`monotone2` +
  `smooth_despike_network`/`_coarsen_network`/`_coarsen_network_3DEP` (still in
  `examples/network/`) → consolidate into a `grlp/` submodule. Design-heavy new
  library API; its own step. Then external pull-ins (Fergus, Yanshui).
- **Still to write** (greenfield `examples/`): `set_S0` and `set_x_bl` demos — the
  v3-feature coverage gap.
- Open: does `deprecated/` (currently `examples/deprecated/`) move to repo root to
  match the new `examples/`/`templates/` siblings?

## Release gate (v3.0.0)

Everything landing in `examples/` must **run on the v3 API** (port the OLD-API
keepers). `templates/` may need data/porting and is documented as such. This
closes issue #25 (reframed: curate a clean subset, not port the heap).
