# Unify single-thread and network modes

## Step 0: Consider a minor release

Before making this larger change, it might be appropriate to make a minor release. This may involve updating machinery for PyPI integration and/or other things that are newer than my latest release code. This will be tagged, released, and thereby pushed to PyPI and picked up by Zenodo.

## Step 1: Generalize the network wtih single-thread as special case

Unify the single-thread and network modes of GRLP by making Network the default and single-thread a special case of this (1-segment = 1-edge network). This might end up being a breaking change depending on what the end user sees, and would be a significant refactor of the code.

The code was originally written for single-thread rivers, wtih the network capability added on. This will require this post-hoc addition to become primary. The class structure might need to change, but I think it might be okay.

## Step 2: NetworkX

Look into `examples/network/perrot_net.py`. This does more than I am looking for here, but I do want us to use NetworkX for the networks. This can replace our self-rolled network approach. Another useful resource could be $HOME/dataanalysis/r.fluvial/v.stream.network. This will expand capacity and reduce dependency.

## Step 3: instead of padded arrays, reach to neighboring arrays

Right now, the solutions rely on padded arrays that are used to build the matrices for the finite-difference solver. This is fine for a single segment. But for multiple segments, information is duplicated. Therefore, incorrect matrix padding can lead to a downstream segment seeing a different set of points upstream than the upstream one does downsteam. These should be linked.

I propose to no longer include padded matrices, and instead compute values for the solution matrices directly by walking up and down the network an pulling those from the arrays from adjacent river segments.

## Step 4: Abstract the processes (Not now, but to set up thinking for later)

Once the network mode is set up, we have the basis to build a generalized network approach to modeling rivers of any sort. This is the goal for the repository at https://github.com/MNiMORPH/FluvTree. This includes sand-bed rivers (very similar to gravel-bed rivers, but a slightly different structure) or stream-power-law bedrock rivers (more like FastScape). We should think about how to organize the softwware architecture, classes, and variable sharing in order to be able to instantiate multiple river types across a network. They might have boundaries between them that will move, and they might all exist on the full network but turn "on" and "off" based on conditions.
 
