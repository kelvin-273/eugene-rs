# Eugene-rs

Eugene-rs is a Rust workspace for modeling and solving plant breeding programs. The
workspace is split into a core Rust library (`eugene_core`) and a Python extension
module (`eugene_pywrapper`) that exposes selected solvers via PyO3.

## Workspace layout

- `eugene_core`: Core domain types (genotypes, gametes, distribute arrays), solvers,
  and helper utilities for generating/visualizing instances.
- `eugene_pywrapper`: PyO3 wrapper that exposes solver APIs to Python as a native
  extension module.

## Building the Rust crates

From the workspace root:

```bash
cargo build
```

To run the binary that solves an example Distribute instance with Distribute A*:

```bash
cargo run -p eugene_core --bin example-to-profile
```

## Using the Rust API

The core crate exposes plant representations in `plants` and solvers in `solvers`.
A simple workflow is to generate a random initial population and run a solver:

```rust
use eugene_core::plants::bit_array::SingleChromGenotype;
use eugene_core::solvers::base_min_generations_segment::breeding_program;
use rand::prelude::*;

let mut rng = thread_rng();
let n_loci = 10;
let n_pop = 6;
let pop_0 = SingleChromGenotype::init_pop_random(&mut rng, n_loci, n_pop);

let result = breeding_program(n_loci, &pop_0);
```

## Python wrapper

The Python wrapper builds a native module named `eugene_pywrapper` and exposes
solver groups under `min_gen` (minimum generations) and `min_cross` (minimum
crossings). To build the extension, use a PyO3-compatible build tool (for example,
`maturin`) or `cargo build -p eugene_pywrapper` if you are embedding the resulting
shared library yourself.

## Additional resources

- The `eugene_core` crate is organized by module:
  - `plants` for genotype/gamete data structures.
  - `solvers` for algorithm implementations.
  - `extra` for instance generation and visualization helpers.

If you are new to the codebase, start with `eugene_core/src/lib.rs` for the module
entry points and `eugene_pywrapper/src/lib.rs` for the Python-facing API surface.
