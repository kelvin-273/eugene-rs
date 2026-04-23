# Eugene-rs

`eugene-rs` is a Rust workspace for modelling and solving plant breeding programs.
It focuses on constructing crossing schedules for single-chromosome diploid
instances, with solvers for objectives such as:

- minimum generations
- minimum crossings
- minimum resources
- distribute-array instances

The workspace contains a core Rust library, `eugene_core`, and a PyO3-based Python
extension, `eugene_pywrapper`, that exposes selected solvers to Python.

## Workspace layout

- `eugene_core`: core plant/genotype data structures, solver implementations,
  instance generators, and visualization/resource helpers
- `eugene_pywrapper`: Python bindings for selected `eugene_core` solver APIs
- `.github/workflows/ci.yml`: CI checks for formatting, linting, and tests

## Quick start

Build the whole workspace:

```bash
cargo build
```

Run the example binary for a distribute-instance minimum-crossings solver:

```bash
cargo run -p eugene_core --bin example-to-profile
```

Run the same checks as CI:

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace
```

## Using the Rust API

The main entry points live under:

- `eugene_core::plants` for genotype and gamete representations
- `eugene_core::solvers` for solver implementations
- `eugene_core::extra` for instance generation, visualization, and resource helpers
- `eugene_core::solution` for crossing-schedule types

A minimal example using the segment-based minimum-generations solver:

```rust
use eugene_core::plants::bit_array::SingleChromGenotype;
use eugene_core::solvers::base_min_generations_segment::breeding_program;
use rand::rng;

let mut rng = rng();
let n_loci = 10;
let n_pop = 6;
let pop_0 = SingleChromGenotype::init_pop_random(&mut rng, n_loci, n_pop);

let result = breeding_program(n_loci, &pop_0);
assert!(result.is_some());
```

Notable solver families in `eugene_core` include:

- `base_min_generations_enumerator`
- `base_min_generations_enumerator_dominance`
- `base_min_generations_segment`
- `base_min_crossings_astar`
- `base_min_crossings_distribute_astar`
- `base_min_resources_sampling`
- `base_min_resources_greedy_dom`
- `greedy_multichrom_ti`

## Using the Python wrapper

The Python extension builds a native module named `eugene_pywrapper`.

If you want to install it into an active virtual environment during development,
`maturin` is the simplest route:

```bash
python -m pip install maturin
maturin develop -m eugene_pywrapper/Cargo.toml
```

You can also build the shared library directly with Cargo:

```bash
cargo build -p eugene_pywrapper
```

The Python module exposes submodules including:

- `min_gen`: minimum-generations solvers (`naive1`, `naive2`, `segment`)
- `min_cross`: minimum-crossings solvers (`astar`, `distribute_astar`)
- `min_res`: minimum-resource heuristics (`sampling`, `greedy_dom`)
- `utils`: helper classes such as recombination-rate wrappers

A small Python example:

```python
import eugene_pywrapper as eugene

n_loci = 4
pop_0 = [
    [[False, False, True, False], [False, False, False, False]],
    [[True, False, False, False], [False, False, False, False]],
]

answer = eugene.min_gen.segment.mingen_answer(n_loci, pop_0, None)
print(answer)
```

For richer return objects, some APIs return a serialized base solution tuple, while
resource-oriented APIs return `PyCrossingSchedule` instances.

## Where to look next

If you are new to the codebase, these files are good starting points:

- `eugene_core/src/lib.rs`: crate-level overview and top-level module exports
- `eugene_core/src/solvers.rs`: solver families available in the Rust crate
- `eugene_pywrapper/src/lib.rs`: Python module layout and exposed submodules
- `eugene_pywrapper/src/solution.rs`: Python-facing wrapper types

## Development notes

- The workspace currently uses Rust 2021 edition.
- CI runs on every push and pull request.
- The repository includes both unit tests and doctests; `cargo test --workspace`
  exercises both.
