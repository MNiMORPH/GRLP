# Network solver: global assembly by walking the topology

GRLP evolves a gravel-bed river long profile by a semi-implicit (Picard-iterated)
solve of a diffusion-like equation for bed elevation `z`. On a single segment
this is a tridiagonal system. On a **network** of segments joined at confluences,
GRLP assembles and solves **one global sparse linear system** over all nodes,
built by walking the channel topology. This note describes that assembly
(`Network.assemble_by_walking` / `Network._evolve_by_walking`).

## What changed, and why

The earlier network solver was **segment-centric**: it built a tridiagonal block
for each segment, arranged those blocks down the diagonal of a large matrix, and
then *separately* wrote off-diagonal entries to couple a segment's boundary node
to its neighbour segment across each junction. That works, but it treats the
segment as the primitive and the junction as a special case bolted on.

The current solver is **node-centric**: every node in the network is a row in one
global matrix, and each node's coupling to its neighbours — *including across a
junction* — comes from the same stencil. A confluence is then just a node with
more than one upstream neighbour; there is no separate block-assembly step. This
is both simpler and closer to how the physics actually connects (elevation
diffuses between adjacent nodes regardless of which "segment" they belong to).

## Global node numbering

Segments are concatenated in segment-ID order into a single global node vector of
length `n = Σ lengths`. Each segment's offset is the cumulative sum of the
preceding segment lengths:

```python
lengths = list(self.list_of_segment_lengths)
starts  = np.cumsum([0] + lengths)[:-1]      # offset of each segment
# node i of segment s has global index:
g = starts[s] + i
```

## Assembly: local stencil → global sparse matrix

The matrix is accumulated as COO triplets (`rows`, `cols`, `vals`) and built once:

```python
rows = []; cols = []; vals = []; RHS = np.zeros(n)
for lp in segments:
    for i in range(len(lp.z)):
        g = starts[lp.ID] + i
        # ... find neighbours by walking; compute center/left/right ...
        rows.append(g); cols.append(g);    vals.append(center)   # diagonal
        rows.append(g); cols.append(up_g); vals.append(left)     # upstream
        rows.append(g); cols.append(dn_g); vals.append(right)    # downstream
LHSmatrix = sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))
```

`up_g` and `dn_g` are the **global indices of the node's real neighbours**, found
by walking the topology:

- **interior node:** the adjacent node in the same segment (`g ± 1`);
- **channel head:** a ghost above the domain, elevation `z[0] + S0·dx` and
  discharge `2·Q[0] − Q[1]` (linear extrapolation);
- **river mouth:** a ghost below the domain, elevation `z_bl` (base level) and
  discharge `2·Q[-1] − Q[-2]`;
- **at a junction:** the neighbour is the *edge node of the adjacent segment*
  (e.g. a confluence node's upstream neighbour is a tributary's last node), whose
  global index is `starts[neighbour_segment] + edge_index`.

Because the neighbour's global index is used directly, the ordinary `left`/`right`
coefficient lands at the cross-segment position `(g, neighbour_global_index)`.
**Junction coupling is therefore an intrinsic off-diagonal entry, not a patch.**

## The stencil is the single-segment stencil

The per-node coefficients (`C1`, the `7/3` terms, the Neumann channel-head and
Dirichlet river-mouth boundary modifications) are **identical** to the
single-segment `LongProfile.build_matrices`. Only the *neighbour lookup* differs
(walk the topology vs. index a padded ghost array). For a one-segment network the
global matrix is exactly the standalone tridiagonal system, bit-for-bit.

## Confluences: a conservative junction cell

A multi-tributary confluence node uses a conservative three-node junction cell:
each junction face carries a single shared sediment-flux conductance used by both
adjacent nodes, so the flux across the face is single-valued and sediment is
conserved to machine precision. It is exact for uniform discharge per segment and
first-order (convergent) for within-segment-varying discharge. This cell writes
its own coupling entries into the same triplet lists.

## Structure of the resulting matrix

The sparsity pattern *is* the channel adjacency graph:

- diagonal: each node's self-coupling;
- off-diagonals: one entry per neighbour (`z_up`, `z_dn`).

For a lone chain the matrix is tridiagonal. For a branching network it is a
mostly-banded matrix (a tridiagonal diagonal block per segment) plus a few
off-diagonal entries wherever an edge node couples across a junction. It is
generally **non-symmetric** (the advective `dQ/dx` term breaks symmetry) but has
the connectivity of a weighted graph Laplacian.

## Relation to standard methods

This is the standard **global assembly** used in finite-element and
finite-volume methods on unstructured grids: compute local (element/node)
contributions and scatter them into a global sparse matrix indexed by global
degrees of freedom, using COO/triplet accumulation. Casting the domain as a graph
whose matrix sparsity equals the node adjacency — diagonal = total coupling,
off-diagonals = coupling to each neighbour — is the structure of a weighted graph
Laplacian, and is exactly how conduction/diffusion on a network is assembled in
circuit *nodal analysis* (Kirchhoff's laws), pipe- and groundwater-flow models,
and heat conduction on graphs. In that framing a river confluence is simply a
graph node of degree > 2, which needs no special treatment for the *matrix
structure* — only for *conservation* at the junction (the flux cell above).

## Time stepping

`Network._evolve_by_walking` runs the semi-implicit scheme: for each step it
freezes the right-hand side at the start-of-step elevation (`zold`) and
Picard-iterates, re-assembling and re-solving `LHSmatrix z = RHS` with `spsolve`
while the coefficients relinearize on the current iterate.
