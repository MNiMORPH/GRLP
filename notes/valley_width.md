# Valley width

Design notes for adding valley-width dynamics to GRLP: a prototype that lets each
river track the width and wall height of the valley it occupies, narrow that
valley as it incises, widen it as it migrates laterally, and partition its
aggraded sediment between channel and overbank deposits.

The prototype is deliberately reduced. It carries valley width as a function of
downstream distance and time, $B(x,t)$, rather than the full cross-valley
geometry $B(x,z)$, and it is built so that a physically complete cross-section
model (TerraPIN) can later drop into the same interface. The reduction will cause
problems in general, but it gives a clean, testable first system.

## 1. Scope and simplifying assumptions

Three deliberate simplifications define the prototype:

1. **Width in time, not in elevation.** We remember $B(x,t)$ but not $B(x,z)$. At
   each downstream position the channel sees one valley width and one wall height;
   it cannot read a width that depends on how deep it has cut. This is the
   limiting assumption, and it is the piece TerraPIN restores. Numerically it is a
   gift: with $B$ frozen within an implicit step and updated between steps, the
   stored-sediment volume stays linear in $z$, so the BDF2 time integration keeps
   its second-order accuracy and the Picard iteration stays trivial.
2. **Vertical walls.** When the channel incises, the valley narrows to the channel
   width $b$; the abandoned floodplain stands as a vertical scarp, and its height
   adds to the wall height $H_v$. Unrealistic, simple, and a clean end-member for
   testing.
3. **A separate channel-deposit fraction.** Rather than fold overbank storage into
   an effective width, we carry the channel-deposit fraction $f_{ch}$ as its own
   variable. This keeps $B$ the true valley width – the quantity migration widens
   and incision narrows – and confines the channel-versus-overbank partition to a
   single, explicit factor, removing any risk of double-counting the stored
   sediment.

## 2. State variables

Each `Segment` carries, at every node, three valley-geometry fields that evolve in
time, plus the rate that drives widening:

| symbol | code | meaning | units |
|---|---|---|---|
| $B$ | `B` | valley width | m |
| $H_v$ | `H_valley` | valley-wall height | m |
| $f_{ch}$ | `f_ch` | channel-deposit fraction | – |
| $\dot\zeta$ | `zetadot` | lateral channel-migration rate | m s⁻¹ |

Here $b$ is the channel width and $h$ the flow depth, both already resolved by
GRLP.

## 3. Processes and their drivers

| process | driver | effect |
|---|---|---|
| incision | $\partial z/\partial t < 0$ | $B \rightarrow b$ (vertical walls); $H_v$ grows |
| lateral migration | $\dot\zeta$ (prescribed, or from $Q_s$) | $B$ widens |
| aggradation | $\partial z/\partial t > 0$ | channel-deposit fraction $f_{ch}$ set by the closure in §5 |

Channel mobility is the driver that ties the last two together: it sets both the
channel-deposit fraction and the frequency with which the river strikes the valley
walls, and it therefore feeds both the deposit partition (§5) and the widening
rate (§6).

## 4. Source models

Three papers supply the physics. All three share the same experimental lineage
(Wickert et al., 2013), so they compose rather than compete.

### 4.1 Channel mobility and deposit fraction – Wickert et al. (2013, *JGR-ES*)

Lateral channel mobility follows from a morphodynamic sediment balance
(their Eq. 30, sediment-flux term retained, avulsion term dropped):

$$\bar{P} = \mu \left[ \frac{1}{1-\lambda_p}\, \frac{Q_s}{b_N\, h_N} + \frac{\bar{P}_b\, H(0)}{h_N} \right],$$

where $\bar{P}$ is the mean lateral mobility, $\mu$ absorbs bank erodibility,
curvature, and vegetation, and $b_N$, $h_N$ are mean channel width and depth.

The channel-versus-overbank partition preserved in the stratigraphy (their
Eq. 31, the reservoir net-to-gross) is

$$f_{ch} = f_w + f_d\,(1 - p_R)\left[ 1 - \exp\!\left(-R\,\frac{h_N}{P_{ob}}\right)\right],
\qquad f_{ob} = 1 - f_{ch},$$

with $f_w$ the wetted (channel) fraction, $f_d = 1 - f_w$ the dry fraction, $p_R$
the reworking asymptote, $R$ the fluvial-surface reworking rate, and $P_{ob}$ the
overbank aggradation rate.

The measurements of mobility, planform overlap, and reworking are all linearly
related to one another as $t \rightarrow 0$ (their §5), which is the hook §5 uses
to close $R$.

### 4.2 Migration rate – Bufe et al. (2019, *ESPL*)

Dimensional analysis of the same experiments gives (their Eq. 7)

$$M_L = k\, \frac{Q_w^{1.13}\, Q_s^{0.13}}{H_b\, D},$$

with $M_L$ the lateral migration rate, $H_b$ the height of the actively reworked
valley walls, and $D$ the grain size. The headline result is that the reworked
sediment volume is constant, so $M_L \propto 1/H_b$: migration slows as the walls
it reworks grow taller. That bank height $H_b$ is our $H_v$.

### 4.3 Valley width – Turowski, Bufe & Tofelde (2024, *ESurf*)

The channel bevels its walls at lateral speed (their Eq. 3)

$$V = \frac{q_L}{H^+},$$

where $q_L$ is the lateral transport capacity per unit channel length and $H^+$ is
the bank height in the direction of motion – again our $H_v$. Treating migration
as a Poisson process with a valley-crossing timescale $\Delta t$, the valley width
is $W = \int_0^{\Delta t} V\,\mathrm{d}t + W_C$ (their Eq. 1), and the unconfined
channel-belt (maximum) width scales with flow depth (their Eq. 8),

$$W_0 = k_0\, h + W_C.$$

Incision raises the wall the channel must cut, so under uplift $U$ (their Eq. 9)

$$\frac{\mathrm{d}H^+}{\mathrm{d}t} = U,$$

which drives the valley below its maximum width (their Eq. 11). Lateral sediment
supply $q_H$ from the walls narrows it further (their Eq. 14),

$$W = W_0 - \frac{q_H}{q_L}\left(W_0 - W_C\right), \qquad q_H < q_L,$$

and the combined uplift-plus-supply width (their Eq. 16) is governed by the
dimensionless **mobility–uplift number** $M_U = q_L / (U\,W_0)$.

## 5. Unified deposit-fraction closure

Wickert's Eq. 31 leaves the reworking rate $R$ as an empirically measured
constant. Turowski's valley-crossing timescale closes it: the reworking rate *is*
the inverse of the time the river takes to cross the valley. Setting $R = 1/\tau$
turns $f_{ch}$ into a closed form with no free constant.

Write the partition as a width. With $f_w = b/B$ and $f_d = (B-b)/B$,

$$B_c \equiv f_{ch}\, B = b + (B - b)\left[1 - \exp\!\left(-R\,\frac{h}{\eta}\right)\right],$$

where $B$ is the valley width, $b$ the channel width, $\eta$ the aggradation rate
(Wickert's $P_{ob}$), and the asymptote $(1-p_R)$ is dropped for the prototype.
$B_c$ is the channel-belt width: the width that must fill with channel deposit per
unit aggradation.

Close $R$ through the crossing timescale and the sediment flux:

$$\tau = \frac{(B-b)\,h}{q_L}, \qquad R = \frac{1}{\tau}, \qquad q_L = \dot{s}\,h,
\qquad \dot{s}_f = \frac{\dot{s}}{B-b},$$

so that $R = \dot{s}_f / f_d = \dot{s}\,B/(B-b)^2$. Correcting the algebra
($B/(B-b)^2 \rightarrow 1/(B-b)$) gives $R = \dot{s}/(B-b)$, and therefore the
closed-form partition

$$\boxed{\; B_c = b + (B - b)\left[1 - \exp\!\left(-\frac{\dot{s}\,h}{\eta\,(B-b)}\right)\right], \qquad f_{ch} = \frac{B_c}{B}. \;}$$

The flux scale $\dot{s}$ is the sediment-flux velocity (§6). The exponent is
dimensionless, and the limits recover Wickert's net-to-gross behavior:

- $\eta \rightarrow 0$ (slow aggradation): $B_c \rightarrow B$, $f_{ch} \rightarrow 1$. The channel reworks the whole valley.
- $\eta \rightarrow \infty$ or $\dot{s} \rightarrow 0$ (burial faster than reworking): $B_c \rightarrow b$, $f_{ch} \rightarrow b/B$. Only the channel footprint is channel deposit.

This closure is derived for aggradation. Incision uses the separate narrowing rule
of §1, assumption 2; the deposit fraction is not the operative variable when the
river is removing material.

## 6. Migration and widening rate

Wickert's Eq. 30, sediment-flux term, is $\bar{P} = k\, Q_s/(b\,h)$ with
$k = \mu/(1-\lambda_p)$. Casting it in Turowski's form – lateral transport capacity
per unit downstream length, $q_L = Q_s/\mathrm{d}x$, and $V = q_L/H^+$ – gives a
single sediment-flux velocity

$$\dot{s} = \frac{Q_s}{h\, \mathrm{d}x},$$

which serves both the migration/widening rate ($V = q_L/H_v = \dot{s}\,h/H_v$) and
the deposit-fraction exponent of §5. One flux scale, not two.

**Caveat on $\mathrm{d}x$.** $Q_s$ is a physical flux and does not shrink under
mesh refinement, so $Q_s/\mathrm{d}x$ is grid-dependent if $\mathrm{d}x$ is the
numerical grid spacing. $\mathrm{d}x$ here must be a fixed physical reference
length (absorbed into the leading constant), not the solver's node spacing.

In the prototype $\dot\zeta$ is a prescribed input (`set_lateral_migration_rate`);
the $Q_s$-based closure above is the first external control to be added.

## 7. Dependency graph

```
 Q_s ──▶ mobility  (P̄, reworking rate R, reworked fraction f_R)
          │
          ├──▶ channel-deposit fraction  f_ch   (§5, Wickert + Turowski τ)
          │
          └──▶ frequency of striking the walls ──▶ valley widening  dB/dt   (§6, Turowski)

 H_v ──▶ widening-rate denominator; grows with incision (dH+/dt = U)
```

Mobility is the driver; the deposit fraction and the widening rate both hang off
it. Incision is a separate branch that narrows $B$ to $b$ and raises $H_v$.

## 8. Software architecture in GRLP

The split follows GRLP's existing **locality** principle: `Segment` holds anything
needing only its own nodes; `Network` handles anything that crosses a junction.

- **`Segment` owns state and rules.** The three fields $B$, $H_v$, $f_{ch}$ and the
  rate $\dot\zeta$; their setters; and the per-node process rules (incision,
  widening, deposit partition) composed in a single `update_valley_state(dt)`. The
  existing volume-storage hooks `valley_width`, `storage_volume`,
  `storage_jacobian` read the current $B$ and $f_{ch}$. Every input the rules need
  – $\partial z/\partial t$, $B$, $b$, $h$, $\mathrm{d}x$, $H_v$, $\dot\zeta$,
  $Q_s$ – is available on the segment.
- **The solver orchestrates.** The between-step update
  `for seg in net.segments: seg.update_valley_state(dt)` runs immediately after
  $\partial z/\partial t$ is published in `evolve` and `evolve_adaptive`.
- **`Network` stays nearly untouched.** It already produces $Q_s$ (stored on each
  segment by `compute_Q_s`), and it hosts whole-domain diagnostics (the sediment
  budget / mass-conservation check). No new network method is needed for the
  valley rules.

The vestigial elevation argument in `valley_width(z)` is the TerraPIN seam: the
prototype ignores it; the full cross-section model uses it.

## 9. Open questions

1. **Flux scale $\dot{s}$.** Fix its exact form and its relation to GRLP's $Q_s$:
   the Turowski form $Q_s/(h\,\mathrm{d}x)$, or the Bufe grain-size form
   $\propto Q/(h D)$.
2. **Incision and $f_{ch}$.** Confirm the deposit fraction is aggradation-only, so
   the incision branch uses the narrowing rule alone.
3. **Reference length $\mathrm{d}x$.** Pin it as a fixed physical length, not the
   grid spacing, and name it explicitly in code.
4. **Which height, where.** Flow depth $h$ in the reworking flux scale; wall height
   $H_v$ in the wall-bevelling denominator. Keep the two distinct.
5. **Valley continuity at confluences (deferred).** Each segment carries its own
   $B$; tributary valleys do not reconcile widths at a junction. A `Network`
   concern if ever wanted, and out of scope for the prototype.
