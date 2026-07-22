#! /usr/bin/env python3
"""
Clean a DEM-extracted river-network graph into a GRLP input.

Reads a GRASS-style node-link JSON (a network already extracted from a DEM),
cleans the along-edge elevation profiles, and writes a GRLP-ready node-link
JSON. Run once per dataset and commit the output: that stored file is the
exact, inspectable input to the model.

The default pass is the monotone pipeline (monotone envelope -> symmetric
running mean -> coarsen). For the alternative despike + diffusive-smooth pass,
import network_preprocessing.clean_edge_values directly.

Example -- reproduce the Whitewater input shipped in templates/network/
(supply the raw GRASS export yourself; it is not committed):

    python3 clean_network.py \\
        grass-whitewater-3DEP-sand.json \\
        whitewater-3DEP-sand-smoothed-coarsened.json \\
        --outlet-node 0 --outlet-source-edge 20 0 --outlet-slope 4e-4 \\
        --smooth-window 101 --keep-fraction 0.1
"""

import argparse

import network_preprocessing as npp


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("infile", help="input GRASS node-link JSON")
    p.add_argument("outfile", help="output GRLP-ready node-link JSON")
    p.add_argument("--node-attr", default="z", help="node elevation attribute")
    p.add_argument("--edge-attr", default="z", help="edge elevation attribute")
    p.add_argument("--smooth-window", type=int, default=101,
                   help="symmetric running-mean window (odd)")
    p.add_argument("--keep-fraction", type=float, default=0.1,
                   help="fraction of along-edge samples to keep when coarsening")
    p.add_argument("--no-coarsen", action="store_true",
                   help="skip the coarsening step")
    p.add_argument("--enforce-child", action="store_true",
                   help="affine-rescale each edge so its end exactly hits the child node")
    # Give the outlet node a real elevation so nan does not propagate upstream:
    # set it from the terminal elevation of a source edge, minus a small slope.
    p.add_argument("--outlet-node", type=int, default=None,
                   help="outlet node id to assign a real elevation")
    p.add_argument("--outlet-source-edge", type=int, nargs=2, default=None,
                   metavar=("U", "V"),
                   help="edge (U V) whose last elevation seeds the outlet")
    p.add_argument("--outlet-slope", type=float, default=4e-4,
                   help="slope offset subtracted from the seed elevation")
    args = p.parse_args()

    G = npp.load_network_json(args.infile)

    if args.outlet_node is not None and args.outlet_source_edge is not None:
        u, v = args.outlet_source_edge
        seed = G.edges[u, v][args.edge_attr][-1]
        G.nodes[args.outlet_node][args.node_attr] = [seed - args.outlet_slope]

    H = npp.preprocess_network(
        G,
        node_attr=args.node_attr,
        edge_attr=args.edge_attr,
        enforce_child=args.enforce_child,
        smooth_window=args.smooth_window,
        coarsen=not args.no_coarsen,
        keep_fraction=args.keep_fraction,
    )

    npp.save_network_json(H, args.outfile)
    print("Wrote %s" % args.outfile)


if __name__ == "__main__":
    main()
