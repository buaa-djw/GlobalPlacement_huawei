# GlobalPlacement Huawei Baseline

This repository is currently a minimal C++17 placement-analysis baseline for reproducing the paper-defined metrics and exact density model from **“Nonsmooth Optimization Method for VLSI Global Placement”**:

- Bookshelf parser front end through the existing Limbo adapter, with a compatibility parser fallback.
- `PlacementDB` with continuous zero-based cell/net/pin/row IDs.
- Exact HPWL evaluation from the input `.pl` placement.
- Equations (3)--(8) exact weighted L1 wirelength model whose value matches exact HPWL.
- Equation (13) exact wirelength subgradient with previous-iteration tie handling.
- Equation (12) nonsmooth objective composition `f = W + lambda P`.
- Equations (9)--(12) exact density model based on exact cell-bin overlap lengths and quadratic overflow penalty.
- Equation (18) OFR.
- Section III-B exact density subgradient with paper-defined random subderivatives at nondifferentiable boundaries.

The current stage intentionally does **not** include initial placement, clustering, multilevel placement, Polak--Ribiere optimization, global optimization, legalization, or detailed placement. Running the program does not modify any cell coordinates.


## Paper Clustering Stage

This repository now includes the Section IV-A / Figure 2 clustering stage from the TCAD paper:

- Figure 2 Modified Best-Choice clustering with a global maximum priority queue.
- Equation (19) clustering score evaluation.
- Lazy invalid/valid recomputation when neighboring connectivity changes.
- Dynamic cluster hypergraph updates after every merge.
- Each hierarchy level targets approximately a fivefold reduction using `ceil(N/5)`.
- The hierarchy stops when movable clusters reach the `current^2` threshold, or earlier when no positive-score candidate remains.
- Fixed objects remain singleton clusters and do not participate in merges.
- Quadratic initialization is still not implemented.
- Polak--Ribiere optimization is still not implemented.
- `global_placer` still prints only HPWL and OFR and does not modify coordinates.

Implementation conventions for paper-underspecified details are documented in `docs/paper_clustering_conventions.md`.

## Build

```bash
rm -rf build
cmake -S . -B build \
  -DLIMBO_ROOT=/home/djw/Desktop/GlobalPlacement_huawei/thrid_party/Limbo
cmake --build build -j$(nproc)
```

ZLIB is required. If Limbo headers/libraries are unavailable, CMake keeps building with the compatibility Bookshelf reader.

## Run

```bash
./build/global_placer --aux <benchmark.aux> --bins 64 64 --target-density 0.9
```

The default result output is restricted to the paper placement metrics:

```text
==================[Placement Result]===================
HPWL <value>
OFR <value>
```

## Tests

```bash
ctest --test-dir build --output-on-failure
```

A benchmark smoke helper is available:

```bash
BENCHMARK_ROOT=testbench/ispd2005 scripts/run_adaptec2_baseline.sh
```

## Coordinates

Bookshelf database coordinates in `PlacementDB` store each cell lower-left corner:

```text
cell.x, cell.y = lower-left corner
```

The paper density formulas use cell centers, so the density code converts before applying equations (9)--(12):

```text
x_i = cell.x + 0.5 * cell.width
y_i = cell.y + 0.5 * cell.height
x_b = 0.5 * (bin.lx + bin.ux)
y_b = 0.5 * (bin.ly + bin.uy)
```

Pin positions are currently evaluated as N-orientation semantics:

```text
pin_x = cell.x + 0.5 * cell.width  + offset_x
pin_y = cell.y + 0.5 * cell.height + offset_y
```

Non-`N` orientations are warned about by the parser path when relevant, but no orientation transform or coordinate modification is performed in this stage.

## Exact Weighted L1 Wirelength

The paper equations (3)--(8) are implemented as an exact weighted L1 wirelength subgradient. The value equals exact HPWL but is evaluated through the paper pair construction rather than by an extrema-only gradient shortcut. Two-pin nets use one unit-weight pair, three-pin nets use all three unordered pairs with weight `0.5`, and each `p > 3` net uses `2p - 3` weighted pairs per axis with weight `1 / (p - 1)`.

For `p > 3`, each axis independently sorts pins by `(coordinate, pin_id)`, selects exactly one minimum boundary representative and one maximum boundary representative, treats all other pins as inner pins, and generates the min-to-all plus max-to-inner pairs. This deterministic tie convention resolves the paper's unspecified tied-extrema representative choice while preserving exact HPWL.

Equation (13) uses actual pin coordinates from `PlacementDB::pinPosition`, not only cell centers. If a canonical pair has equal current coordinates, the caller-provided previous pin-coordinate snapshot selects the Section III-A direction rule: previous greater samples `theta1 ~ Uniform[0, pi/3]`, previous less samples `theta1 ~ Uniform[2pi/3, pi]`, and no previous/equal previous uses `theta1 = 0`. The caller supplies `std::mt19937_64`, so a fixed seed and fixed previous snapshot reproduce random tie subgradients. Fixed/non-movable objects still participate in wirelength value and pair classification, but their gradient entries are zero.

## Nonsmooth Objective

Equation (12) is available as an internal evaluator:

```text
f(x,y) = W(x,y) + lambda * P(x,y)
```

`W` is the exact weighted L1 wirelength value, `P` is the unscaled quadratic density penalty, and the returned subgradient is the direct mathematical sum of wirelength subgradient plus `lambda` times density subgradient. The objective evaluator calls ExactWirelengthSubgradient before ExactDensitySubgradient when a shared `std::mt19937_64` is supplied, so a fixed seed is reproducible. Objective and density penalty remain internal optimization quantities; the default executable still prints only HPWL and OFR.

## Exact Density Definition

For each axis, equations (9) and (10) are implemented as the exact overlap length between a cell and a bin using center coordinates. For each bin:

```text
D_b = sum_i overlap_x(i,b) * overlap_y(i,b)
capacity_b = target_density * bin_width * bin_height
overflow_b = max(0, D_b - capacity_b)
P = sum_b overflow_b^2
```

The OFR metric follows equation (18):

```text
OFR = sum_b overflow_b / sum_i (cell.width_i * cell.height_i)
```

The denominator is the raw area of all cells before clipping and includes movable standard cells, fixed standard cells, `terminal`, `terminal_NI`, and cells completely outside the placement core.

Fixed and non-movable objects participate in `D_b`, overflow, penalty, and the OFR denominator. Their density subgradient entries are always zero; only movable standard cells receive equation (12) penalty gradients.

At nondifferentiable overlap boundaries, Section III-B is followed directly: positive-side boundaries use `theta2 ~ Uniform[pi/2, pi]`, negative-side boundaries use `theta3 ~ Uniform[0, pi/2]`, and the equal-size coincident case uses `theta ~ Uniform[0, pi]`. The caller supplies the `std::mt19937_64` seed so repeated runs with the same seed are reproducible.
