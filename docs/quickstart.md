# Quickstart

## A single long profile

Build a `LongProfile`, bring it to steady state, then perturb its upstream
sediment supply and watch the profile respond. This mirrors
`examples/one_dimensional/Basic1D.py`.

```python
import numpy as np
import grlp

# --- Physical setup (Río Santa Cruz-like numbers) --------------------------
Q  = 700.          # m3/s, channel-forming water discharge
z0 = 180.          # m, elevation at the upstream end
z1 = 0.            # m, base level
L  = 250e3         # m, valley length
S0 = (z0 - z1) / L # equilibrium slope
B  = 10000.        # m, valley width
U  = 1e-4 / 3.15e7 # m/s, uplift rate

# --- Build the profile -----------------------------------------------------
lp = grlp.LongProfile()
lp.basic_constants()
lp.bedload_lumped_constants()
lp.set_hydrologic_constants()

lp.set_x(dx=1000, nx=251, x0=0)   # downstream distance grid
lp.set_z(S0=S0, z1=z1)            # initial (linear) profile
lp.set_Q(Q)                       # water discharge
lp.set_B(B)                       # valley width
lp.set_niter(3)                   # Picard iterations per step
lp.set_z_bl(z1)                   # downstream base level (Dirichlet BC)
lp.set_uplift_rate(U)

# Upstream sediment-supply boundary (Neumann): the supply consistent with S0
Qs0 = lp.k_Qs * lp.Q[0] * S0 ** (7 / 6.)
lp.set_Qs_input_upstream(Qs0)

# --- Evolve to steady state, then perturb ----------------------------------
lp.evolve_threshold_width_river(10, 1e14)   # long steps -> equilibrium
z_eq = lp.z.copy()

lp.set_Qs_input_upstream(Qs0 * 4)           # quadruple the sediment supply
for _ in range(50):
    lp.evolve_threshold_width_river(1, 1000 * 3.15e7)  # 1000-yr steps

# lp.z now holds the transiently aggraded profile; compare with z_eq.
```

`evolve_threshold_width_river(nt, dt)` advances `nt` steps of size `dt` seconds.
The upstream boundary can be forced by sediment supply (`set_Qs_input_upstream`)
or directly by slope (`set_S0`); the downstream boundary is base level
(`set_z_bl`), which may also be moved horizontally with `set_x_bl` (e.g. to track
a shoreline along a continental shelf — see the `RioSantaCruz_set_S0_set_x_bl`
example).

## A river network

Every GRLP solution is a network solution; a single long profile is the trivial
one-edge case. Generate a random Shreve network, evolve it, and plot it:

```python
import grlp

# generate_random_network returns (network, topology).
#   magnitude     : number of channel heads (network size)
#   max_length    : length of the longest source-to-outlet path [m]
#   mean_discharge: sets segment discharges from drainage area
net, topo = grlp.generate_random_network(
    magnitude=8,
    max_length=2.0e4,
    mean_discharge=10.,
)
net.set_niter(3)
net.evolve_threshold_width_river_network(nt=200, dt=3.15e11)
net.plot()   # branching planform
```

See {doc}`examples` for point-source (landslide) forcing, base-level change,
and network analysis (Hack parameters, Strahler orders, Horton ratios).
