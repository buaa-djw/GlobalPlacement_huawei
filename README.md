# GlobalPlacement_huawei

## Project Goal

`GlobalPlacement_huawei` is a C++ prototype for reading ISPD2005 Bookshelf placement benchmarks and preparing data for later global-placement research. The current executable validates the data layer by building an in-memory placement database and printing a summary.

## Current Scope: PlacementDB

This revision implements only the PlacementDB layer and the Bookshelf parser adapter. It does **not** implement a placement optimizer, HPWL optimization, density optimization, legalization, routing, or visualization.

## Build Environment

The project is intended for Ubuntu 22.04 with a C++17 compiler and CMake 3.16 or newer.

## Dependencies

- Ubuntu 22.04
- C++17
- CMake
- Limbo Bookshelf parser
- ISPD2005 Bookshelf benchmark, for example `adaptec1`

Limbo is optional at configure time only so that the repository can still build in minimal CI containers. For production use, pass `-DLIMBO_ROOT=/path/to/limbo`; when the Limbo header is found, `LimboBookshelfAdapter` calls `BookshelfParser::read` and receives nodes, nets, pins, placement, and rows through Limbo callbacks. Without Limbo headers, the adapter prints a warning and uses a small compatibility reader for local smoke testing.

## Directory Structure

```text
GlobalPlacement_huawei/
├── CMakeLists.txt
├── README.md
├── main.cpp
├── include/
│   ├── db/
│   │   ├── Cell.h
│   │   ├── Net.h
│   │   ├── Pin.h
│   │   ├── PlacementDB.h
│   │   └── Row.h
│   └── parser/
│       └── LimboBookshelfAdapter.h
├── src/
│   ├── db/
│   │   └── PlacementDB.cpp
│   └── parser/
│       └── LimboBookshelfAdapter.cpp
├── tests/
│   └── test_placement_db.cpp
└── scripts/
    └── run_adaptec1.sh
```

## Data Structures

All geometric values are stored as `double` in Bookshelf database units. All IDs are zero-based contiguous indices owned by `PlacementDB`.

### Cell

`Cell` represents a Bookshelf node from `.nodes` and its placement from `.pl`.

- `id`: contiguous cell ID in `PlacementDB`.
- `name`: original Bookshelf node name.
- `width`, `height`: node dimensions.
- `x`, `y`: lower-left placement coordinate.
- `is_terminal`: true for Bookshelf terminals.
- `is_fixed`: true when the node is fixed by terminal status or `.pl` fixed flags.
- `pin_ids`: IDs of pins attached to this cell.

### Pin

`Pin` represents one pin entry from `.nets`.

- `id`: contiguous pin ID in `PlacementDB`.
- `cell_id`: owning cell ID.
- `net_id`: owning net ID.
- `cell_name`: original Bookshelf cell name copied for diagnostics.
- `net_name`: original Bookshelf net name copied for diagnostics.
- `offset_x`, `offset_y`: Bookshelf/Limbo pin offset values; the adapter stores the callback values without reinterpretation.
- `direction`: pin direction string such as `I`, `O`, or `B`; it may be empty if the source does not provide one.

### Net

`Net` represents one Bookshelf net from `.nets`.

- `id`: contiguous net ID in `PlacementDB`.
- `name`: original Bookshelf net name.
- `pin_ids`: IDs of all pins connected to the net.

### Row

`Row` represents one placement row from `.scl`.

- `id`: contiguous row ID in `PlacementDB`.
- `y`: row y coordinate.
- `height`: row height.
- `site_width`: width of one placement site.
- `site_spacing`: site pitch.
- `x_start`, `x_end`: horizontal row interval.
- `num_sites`: number of placement sites in the row.

### PlacementDB

`PlacementDB` owns vectors of cells, nets, pins, and rows plus name-to-ID maps. It validates IDs, rejects duplicate cells/nets, reports unknown netlist cells as errors, warns on unknown `.pl` cells, and guarantees that every pin is linked from both its cell and its net.

## Main Interfaces

Core APIs are declared in `include/db/PlacementDB.h`:

- `int addCell(const std::string& name, double width, double height, bool is_terminal)`: create a cell and return its ID.
- `int addNet(const std::string& name)`: create a net and return its ID.
- `int addPin(int cell_id, int net_id, double offset_x, double offset_y, const std::string& direction)`: create a pin, append it to the owning cell and net, and return its ID.
- `int addRow(const Row& row)`: append a placement row and assign its ID.
- `bool hasCell(const std::string& name) const`, `bool hasNet(const std::string& name) const`: query whether names exist.
- `int getCellId(const std::string& name) const`, `int getNetId(const std::string& name) const`: map Bookshelf names to DB IDs or throw on unknown names.
- `Cell& cell(int id)`, `Net& net(int id)`, `Pin& pin(int id)`: mutable object access with ID checking.
- `const std::vector<Cell>& cells() const`, `const std::vector<Net>& nets() const`, `const std::vector<Pin>& pins() const`, `const std::vector<Row>& rows() const`: read-only vector access for future optimizers.
- `void setCellLocation(const std::string& name, double x, double y, bool fixed)`: update placement and fixed status from `.pl`.
- `void printSummary(std::ostream& os) const`: print counts, region bounds, and sample cells/nets.

Future modules should use IDs and the const vector accessors for fast traversal. They should not depend on Limbo callback types.

## Parser Flow

```text
ISPD2005 Bookshelf .aux
→ Limbo Bookshelf parser
→ LimboBookshelfAdapter callbacks
→ PlacementDB cells/nets/pins/rows
→ main.cpp summary output
```

`main.cpp` only parses command-line arguments, creates `PlacementDB`, invokes `LimboBookshelfAdapter::read(aux_path)`, and prints the summary.

## How to Build

```bash
rm -rf build
mkdir build
cd build
cmake .. -DLIMBO_ROOT=/path/to/limbo
make -j$(nproc)
```

If Limbo is installed in a nonstandard location, `LIMBO_ROOT` should point to the directory that contains the `limbo/` include tree or its immediate source include root.

## How to Run on adaptec1

```bash
./global_placer --aux /path/to/ispd2005/adaptec1/adaptec1.aux
```

This repository also includes a helper script:

```bash
./scripts/run_adaptec1.sh testbench/ispd2005/adaptec1/adaptec1.aux
```

## Expected Summary Output

The output format is:

```text
========== PlacementDB Summary ==========
Number of cells:
  total: ...
  movable: ...
  terminal/fixed: ...

Number of nets:
  total: ...

Number of pins:
  total: ...
  average pins per net: ...
  max pins per net: ...

Rows:
  total: ...

Placement region:
  min x: ...
  max x: ...
  min y: ...
  max y: ...

Sample cells:
  [0] name=..., width=..., height=..., x=..., y=..., fixed=...

Sample nets:
  [0] name=..., degree=...
========================================
```

## Notes for Future Development

- Optimizers should consume `PlacementDB::cells()`, `PlacementDB::nets()`, and `PlacementDB::pins()` rather than Limbo data structures.
- HPWL evaluators can traverse each `Net::pin_ids`, then use `Pin::cell_id` and `Cell` coordinates.
- Density and legalization modules can read `Cell` sizes, fixed flags, and `Row` geometry.
- Any future algorithm should be added outside this PlacementDB/parser scope and should preserve the ID stability guarantees described above.
