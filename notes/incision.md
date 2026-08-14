# Incision

From Turowski et al. (2025, ESurf), the rate of valley widening is given by:
    
```math
\frac{\text{d}B}{\text{d}t} = PV \tag{T9}
```

where $B$ is valley width, $t$ is time, $P$ is the fraction of time the river spends at the outer walls, and $V$ is the lateral migration rate at the valley walls.

The fraction of time at the valley walls is in turn given by

```math
P = \frac{B_0 - B}{(B_0 - b)(1 - \frac{H}{h})} \tag{T17}
```

while the migration velocity is

```math
V = \frac{q_L}{H} \tag{T1}
```

Here, $B_0$ is the unconfined channel-belt width (a theoretical width for the channel belt a river with this hydrology would create without any existing valley walls), $b$ is channel width, $H$ is the height of the valley walls, $h$ is the channel depth and $q_L$ is the lateral transport capacity.

Furthermore, the unconfined channel-belt width is given by

```math
B_0 = \frac{ch}{k_T} + b \tag{T4}
```

where $c = 2.2285$ and $k$ is a dimensionless constant.

We want to solve T19 for $B$. To avoid making edits to the main matrix to allow for changing $B$ through time (maybe would make sense long term), I propose a separate centred difference:
    
```math
\frac{B^{n+1} - B^n}{\delta t} = [PV]^{n+\frac{1}{2}}
```

which we can solve as part of the Picard iteration, with a first guess of $B^{n+1} = B^n$. We can then set `segment.B` to $B^{n+\frac{1}{2}}$. More explicitly we can write

```math
B^{n+\frac{1}{2}} = \frac{3}{2}B^n + \delta t [PV]^{n+\frac{1}{2}} =  \frac{3}{2}B^n + \delta t \frac{ [PV]^n + [PV]^{n+1} }{2}
```

**Important downside is that `segment.B` will lag behind other variables by half a timestep.** Maybe it would be an option to create a new `B` that goes into the solver to allow `segment.B` to keep up. For non-variable width cases these would just the same.

We will need a set of functions along the lines of (some will also be relevant for the aggradation case):
    
- `compute_q_L` (modified version of Aaron's equation?)
- `compute_B_0` (T4)
- `compute_widening_time_fraction` (T17)
- `compute_lateral_migration_velocity` (T1)

and to keep track of the values at the previous and new timestep within the iteration.

We will also need code to recognise the onset of incision, reset $B$ to $b$, and start tracking with amount of incision $H$.