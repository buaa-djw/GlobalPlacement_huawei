# Paper clustering conventions

## A. Paper-specified behavior

- Implements Figure 2 Modified Best-Choice control flow with a global maximum priority queue.
- Implements Equation (19) for clustering score.
- Uses the maximum score candidate for each merge, with lazy invalid/valid recomputation.
- Each clustering process targets an approximately fivefold block reduction.
- The placement hierarchy stops at the `current^2` movable-cluster threshold.
- Merging dynamically updates the cluster hypergraph/netlist and removes internal nets.

## B. Paper-underspecified implementation conventions

These items are engineering conventions, not claimed as fully specified by the TCAD paper.

- `A(a)=1/a`; the paper only says `A(a)` is inversely proportional to area, and a global proportional constant is ranking-equivalent within one PQ.
- `D_jk`/`E_j` use a paper-referenced implementation convention; the TCAD paper itself does not fully specify the connectivity expansion. For each active net with `p_e` unique active clusters, a shared pair receives `1/(p_e-1)`, and `E_j` is the sum over active neighbors.
- Score ties are deterministic: larger score first, then smaller canonical first cluster ID, then smaller canonical second cluster ID.
- Fixed objects are singleton clusters, remain in nets, affect connectivity, and never participate in movable merge candidates.
- A level target is `ceil(N/5.0)`.
- Cluster IDs are monotonic: original IDs first, merged clusters appended.
- Coarse geometry is deferred; clusters store only area as the sum of member areas.
