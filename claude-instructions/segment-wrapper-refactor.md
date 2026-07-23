# Segment / LongProfile-wrapper refactor (v3, pre-tag)

## Motivation

Today a single `LongProfile` plays two roles at once: it is a **member** that a
`Network` contains, *and* — through its solve methods (`evolve_threshold_width_river`,
`compute_Q_s`, `slope_area`) — a **self-contained whole** that builds
`Network([self])` around itself to get solved. A part instantiating its own
container is an inverted composition: it works and is simple internally, but it
reads as confusing *in use* ("the same object lies inside the network and is also
callable on its own and sticks itself into the network").

Fix: **split the two roles into two objects.**

- **`Segment`** — the obviously-embedded object. A pure network member: data +
  segment-level configuration + the per-node coefficients the solver consumes.
  It does **not** solve itself.
- **`LongProfile`** — a lightweight 1-D convenience wrapper. It *owns* a
  `Segment` (and the one-edge `Network` around it) and exposes the friendly
  standalone API by delegation. The wrapping becomes explicit and honest: a
  convenience whole that *composes* a part.

`Network` is built from `Segment`s. The solver (`grlp.solver`) is unchanged — it
already duck-types the network's members; only their class name changes.

User story, now unambiguous:
- **1-D** → `lp = grlp.LongProfile(...)` (never appears in a `Network`).
- **Network** → `grlp.Network([seg_a, seg_b, ...])` with `Segment`s; drive the
  `Network`. (`Segment` has no `evolve`, so it *cannot* be solved on its own.)

## Naming

- Class `LongProfile` (member) → **`Segment`**. Matches existing code
  (`upstream_segment_IDs`, `channel_head_segment_IDs`) and the fluvial term.
- New wrapper keeps the name **`LongProfile`** — a "long profile" *is* the
  single-river 1-D object, and `grlp.LongProfile()` stays the 1-D entry point.
- Variables: `lp` → **`seg`** wherever they name a segment. Tightens the
  `seg` / `upseg` / `downseg` neighbour notation for readability.
- `Network.list_of_LongProfile_objects` → **`Network.segments`** (keep a
  read-only `list_of_LongProfile_objects` property aliased for one release if any
  external code relies on it; otherwise drop — pre-tag, so likely drop).

## Method allocation

**On `Segment`** (everything that does *not* require a network):
- `__init__`, `set_ID`, `set_upstream_segment_IDs`, `set_downstream_segment_IDs`
- `basic_constants`, `bedload_lumped_constants`, `set_hydrologic_constants`
- `set_intermittency`, `set_x`, `set_z`, `set_A`, `set_Q`, `set_B`
- `set_uplift_rate`, `set_source_sink_distributed`, `set_Sternberg_gravel_loss`
- `set_niter`, `set_Qs_input_upstream`, `set_S0`, `set_z_bl`, `set_x_bl`, `set_bl`
- `build_LHS_coeff_C0` (the coefficient the solver reads)
- the analytical / spectral 1-D diagnostics that are pure math on the segment's
  own arrays: `analytical_threshold_width*`, `compute_channel_width`,
  `compute_flow_depth`, `compute_diffusivity`, `compute_length`,
  `compute_equilibration_time`, `compute_e_folding_time`, `compute_wavenumber`,
  `compute_*_series_terms`, `compute_*_gain`, `compute_*_lag`.
  (These don't self-wrap, so they may stay; the wrapper exposes them by
  forwarding. Some depend on `Q_s`/`S` being set first — see below.)

**Removed from `Segment`, provided by the wrapper / `Network`** (they require a
network walk):
- `evolve_threshold_width_river` — build/drive `Network([segment])`.
- `compute_Q_s`, `slope_area` — delegate to the owned `Network` (whose
  `compute_Q_s` walks the topology and sets each segment's `S`/`Q_s`).

A bare `Segment`'s `Q_s`/`S` are set by whatever network solves it (the wrapper's
network for 1-D, or the real `Network` when it's a member). So
`Segment.compute_channel_width` etc. work once a solve/`compute_Q_s` has run —
they just can't *originate* the walk.

## Delegation mechanism (the "lightweight" wrapper)

`LongProfile` holds `self.segment` (a `Segment`) and a cached
`self._network = Network([self.segment])` built lazily on first solve. It:
- forwards attribute access and configuration to the segment via `__getattr__`
  (so `lp.set_x(...)`, `lp.z`, `lp.analytical_threshold_width()` all reach the
  segment) — this keeps the wrapper genuinely lightweight, no ~25-method
  boilerplate;
- defines the three network operations explicitly (`evolve_threshold_width_river`,
  `compute_Q_s`, `slope_area`), each driving `self._network`.

The cached network holds the segment *by reference*, so later `set_*` calls
mutate the same segment the network sees — no rebuild needed each `evolve`
(an incidental win over today's rebuild-every-call). `__getattr__` only fires for
names the wrapper doesn't define, so the three explicit methods win; everything
else falls through to the segment.

Open question to settle in code: `__getattr__` forwarding (concise, slightly
magic) vs. a handful of explicit `set_*` delegations (verbose, explicit). Lean
`__getattr__` for "lightweight"; revisit if readability suffers.

## Execution plan (tests green at every commit)

1. **Introduce `Segment` behind an alias.** Rename `class LongProfile` →
   `class Segment`; add `LongProfile = Segment`. Pure rename — behaviour
   identical, all existing code works via the alias. Suite green. *(small, safe)*
2. **Split into pure `Segment` + real `LongProfile` wrapper.** Remove the three
   self-wrapping methods from `Segment`; replace the alias with the composition
   wrapper. Update network-construction sites (tests/examples that build a
   `Network`) to pass `Segment`s. 1-D sites using `grlp.LongProfile()` keep
   working through the wrapper. Suite green. *(the delicate step)*
3. **Readability renames.** `lp` → `seg`; `list_of_LongProfile_objects` →
   `segments`; tidy `seg`/`upseg`/`downseg`. Suite green.
4. **Docs / examples / changelog.** Update theory/api/tutorial narrative and the
   network-building examples to `Segment`; write the usage rule (1-D → LongProfile,
   network → Segment+Network); CHANGELOG entry (breaking: `Network` takes
   `Segment`s). Suite green.

Granular commits within each phase. The solver and the numerics do not change;
this is an interface/composition refactor, so the golden-master characterization
tests are the safety net that it stays behaviour-identical.

## Compatibility

Pre-v3.0.0 tag, so a breaking rename is in-window (this is exactly when to do it;
after release it would force a v4). The common **1-D** path (`grlp.LongProfile()`)
is preserved by name. The breaking change is **network construction**: callers
now build `Network([Segment, ...])` rather than `Network([LongProfile, ...])`.
