# Theory and governing equations

GRLP evolves the long profile of a transport-limited gravel-bed river. This page
gives the governing equations **as amended by the 2020 corrigendum** to Wickert &
Schildgen (2019) — the forms the code actually solves. The papers are open access
(CC-BY); see {doc}`references`. Fuller, source-linked notes on all of the
literature live in `docs/literature/README.md` in the repository.

## Symbols

| Symbol | Meaning | Code | Typical |
|---|---|---|---|
| $z$ | valley-floor (bed) elevation | `z` | — |
| $x$ | down-valley distance | `x` | — |
| $Q_s$ | bed-load sediment discharge | — | — |
| $Q$ | channel-forming water discharge | `Q` | — |
| $B$ | valley width | `B` | — |
| $U$ | uplift / subsidence (source–sink) | `set_uplift_rate` | — |
| $\lambda_p$ | sediment porosity | `lambda_p` | 0.35 |
| $I$ | intermittency | `intermittency` | 0.01–1 |
| $\mathbb{S}$ | sinuosity | `sinuosity` | $\ge 1$ |
| $k_{Q_s}$ | lumped transport coefficient | `k_Qs` | $\approx 0.041$ |

## Sediment transport

For a near-threshold, self-formed (Parker, 1978) **equilibrium-width** gravel
channel, the bed-load sediment discharge is

$$
Q_s = -\frac{k_{Q_s}\, I}{\mathbb{S}^{7/6}}\; Q\;
      \frac{\partial z}{\partial x}\,
      \left|\frac{\partial z}{\partial x}\right|^{1/6}.
$$

A defining feature of the equilibrium-width assumption is that **grain size
cancels**: the specific transport rate scales as $q_s \propto D^{3/2}$ while the
channel width scales as $b \propto D^{-3/2}$, so no grain size enters $Q_s$. The
lumped coefficient $k_{Q_s} = k_{q_s}\,k_b \approx 0.041$ is built from first
principles (Wong & Parker, 2006; Parker, 1978) in
`LongProfile.bedload_lumped_constants()`.

## Mass conservation (Exner)

$$
\frac{\partial z}{\partial t}
  = -\frac{1}{B(1-\lambda_p)}\frac{\partial Q_s}{\partial x} + U.
$$

Valley width $B$ enters **only as a local prefactor** $1/B$, not through its
downstream derivative. (Keeping $B$ inside the derivative was the error fixed by
the corrigendum — see below.)

## Long-profile evolution equation

Combining the two gives a nonlinear diffusion equation:

$$
\frac{\partial z}{\partial t}
  = \frac{k_{Q_s}\, I}{\mathbb{S}^{7/6}(1-\lambda_p)}
    \left|\frac{\partial z}{\partial x}\right|^{1/6}
    \left[\frac{7}{6}\,\frac{Q}{B}\,\frac{\partial^2 z}{\partial x^2}
          + \frac{1}{B}\,\frac{\partial Q}{\partial x}\,
            \frac{\partial z}{\partial x}\right] + U.
$$

The first bracketed term is diffusive; the second is an advective term arising
from downstream discharge growth. There is no valley-width-derivative term.

## Boundary conditions

Upstream (Neumann, sediment flux)
: A prescribed input sediment discharge sets the boundary slope,
  $S_0 = -\,\mathrm{sgn}(Q_s)\,\mathbb{S}\,
  \bigl(Q_s /(k_{Q_s} I Q)\bigr)^{6/7}$.
  A higher sediment-to-water supply ratio steepens the boundary (Lane's balance).
  Set with `set_Qs_input_upstream(Q_s_0)`, or slope-first with `set_S0(S0)`.

Downstream (Dirichlet, base level)
: Elevation is fixed at the outlet with `set_z_bl(z_bl)`. The outlet may also be
  moved horizontally with `set_x_bl` — for example to follow a shoreline
  migrating along a continental-shelf gradient (Ruby et al., 2026).

## Numerical method

The equation is solved by a **semi-implicit finite-difference** scheme. The
diffusive part is a tridiagonal system solved directly; the weak nonlinearity
from $|\partial z/\partial x|^{1/6}$ (and any nonlinear $Q(x)$, $B(x)$) is
resolved by **Picard iteration** (`set_niter`). On a network, one global sparse
system is assembled over all nodes by walking the channel topology
(`grlp.solver`).

## Steady-state and slope–area

The analytical steady-state profile between endpoints $(x_0, z_0)$ and
$(x_1, z_1)$, useful as a verification benchmark, is

$$
z = (z_1 - z_0)\,
    \frac{x^{a} - x_0^{a}}{x_1^{a} - x_0^{a}} + z_0,
    \qquad a = 1 - \tfrac{6}{7}P_{x,Q},
$$

with the exponent depending on the distance→discharge power $P_{x,Q}$ alone. The
corresponding slope–area relation is
$S = S_0 (A_0/A)^{(6/7)P_{A,Q}}$, so channel concavity
$\theta = \tfrac{6}{7}P_{A,Q}$ is set by the drainage-area→discharge exponent.

## Linearized response: gain and lag

Under small sinusoidal forcing, the equation linearizes to a constant-coefficient
diffusion equation for the elevation perturbation, with a **sediment-transport
diffusivity**

$$
\kappa = \frac{7}{6}\,
         \frac{k_{Q_s}\, I\, \overline{Q}}
              {\mathbb{S}^{7/6} B(1-\lambda_p)}\,
         \left|\frac{\partial\bar z}{\partial x}\right|^{1/6},
\qquad
T_{eq} = \frac{L^2}{\kappa}.
$$

The equilibration time $T_{eq}$ controls the response: fast forcing
($P \ll T_{eq}$) is strongly damped, slow forcing ($P \gg T_{eq}$) is nearly in
phase with gain $\to 6/7$. GRLP computes these directly
(`compute_diffusivity`, `compute_equilibration_time`, `compute_z_gain`,
`compute_z_lag`, `compute_Qs_gain`, `compute_Qs_lag`); the framework is from
McNab et al. (2023).

## The 2020 corrigendum

The original Wickert & Schildgen (2019) equations mistakenly kept valley width
$B$ inside the spatial derivative of the Exner equation, which incorrectly
implied that downstream valley widening or narrowing could by itself drive
aggradation or incision. The corrigendum removes every valley-width-derivative
term (so $B$ is a $1/B$ prefactor), restores an omitted $\mathbb{S}^{7/6}$
sinuosity factor in the numerical scheme, corrects the upstream boundary
condition's sign, and redefines concavity as $\theta = \tfrac{6}{7}P_{A,Q}$. The
GRLP release that brought the code into line with the corrected equations was
v1.4.1. **Use the corrected forms on this page.**

## Networks

A network is a set of long-profile segments joined at confluences, each obeying
the evolution equation above, coupled by three junction rules (McNab et al.,
2025): (1) a segment's upstream sediment supply is the **sum** of its upstream
segments' sediment discharges; (2) elevation is **continuous** across every
junction; (3) water discharge **steps up** at each tributary junction. Random
network topologies are generated as Shreve (1966, 1974) random binary trees
(`generate_random_network`). The integrated network response resembles a single
profile whose length is the network's **mean inlet-to-outlet distance**, not its
longest path.
