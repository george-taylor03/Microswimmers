# MicroSwimmers.jl — Code Archive

This repository contains the Julia packages and analysis scripts accompanying the paper. It provides:

- **MicroSwimmers.jl** — the core simulation package for microswimmer hydrodynamics (boundary element / regularised Stokeslet method)
- **MicroSwimmersPlots.jl** — visualisation utilities built on Makie
- **BeatPlanes** — organism definitions, flagellar beat models, and simulation runners specific to the paper
- **scripts/** — figure-generation scripts and analysis routines
- **results/** — pre-computed simulation data (`.jld2`) used by the figure scripts

---

## Citation

If you use this code in your research please cite the accompanying paper, titled *"Simulation-driven discovery of morphology-function relationships in microswimmers"*. If this work is not yet published, a preprint is available on [bioRxiv](https://doi.org/10.64898/2026.06.01.729272).

---

## Requirements

- **Julia 1.12** or later. Download from [https://julialang.org/downloads/](https://julialang.org/downloads/) and follow the platform installer instructions. On Windows or macOS the installer adds `julia` to your PATH automatically; on Linux add the `bin/` directory to your PATH manually.

---

## Installation

1. Download and extract this archive.
2. Open a terminal in the root directory of the archive.
3. Start Julia:
   ```
   julia
   ```
4. Activate the environment and install all dependencies:
   ```julia
   using Pkg
   Pkg.activate(".")
   Pkg.instantiate()
   ```
   This reads `Manifest.toml` and installs the exact package versions used to produce the paper results. It may take several minutes on first run.

---

## Reproducing paper figures

Each figure script loads the pre-computed simulation data from `results/` and produces the corresponding figure (finishing touches were added in Inkscape). Run from the archive root:

| Script | Output |
|---|---|
| `scripts/fig1.jl` | Figure 1 — package overview (organisms, problems, flow fields) |
| `scripts/fig2.jl` | Figure 2 — Chlamydomonas beat-plane analysis |
| `scripts/fig3.jl` | Figure 3 - Filter feeding |
| `scripts/fig4.jl` | Figure 4 - Morphological transition |

```julia
include("scripts/fig1.jl")
```
The SI videos (`SI/SI_1.mp4`, `SI/SI_2.mp4`) can be regenerated with:
```julia
include("SI/SI_1.jl")
include("SI/SI_2.jl")
```
You can move the camera around (hold left-click to rotate, right-click to translate).

---

## Getting started

`scripts/getting_started.jl` is a self-contained introduction to MicroSwimmers.jl. It walks through:

- Constructing flagellum models (`PlanarFlagellum`) and discretising them
- Solving the instantaneous swimming problem (`SwimmingProblem`)
- Computing velocity fields (`PlanarVelocityField`, `TimeAveragedPlanarVelocityField`)
- Running time-dependent trajectory simulations (`SwimmingTrajectoryProblem`)
- Building a cell body and attaching flagella to create a `Flagellate`

This tutorial should ideally be run line-by-line. You can do this in Visual Studio Code on Windows by pressing shift+enter on a line to run it. Copying and pasting into the REPL will work if all else fails.

The script uses `GLMakie` for interactive figures; a display (not a headless server) is required.


---

## Repository structure

```
.
├── MicroSwimmers/          # Core package: BEM solver, flagellum/cell-body models, trajectory integration
├── MicroSwimmersPlots/     # Visualisation: viz(), stream(), animate(), mesh generation
├── BeatPlanes/             # Paper-specific: organism constructors (chlamy, filter feeder, excavate→chlamy)
│   └── src/
│       ├── chlamy.jl           # Chlamydomonas model and ChlamySimulation structure
│       ├── filter_feeder.jl    # Filter feeder model etc.
│       ├── excavate_2_chlamy.jl
│       └── simulations.jl      # parameter sweep infrastructure
├── scripts/
│   ├── getting_started.jl  # Introductory tutorial
│   ├── chlamy_analysis.jl  # Analysis routines loaded by figure scripts
│   ├── filter_feeder_analysis.jl
│   ├── excavate2chlamy_analysis.jl
│   ├── convergence_checks.jl
│   └── fig1.jl – fig4.jl   # generate the paper figures
├── results/                # Pre-computed .jld2 simulation outputs
├── SI/                     # Supplementary video files and generation script
├── Project.toml            # Top-level environment (all three packages + script dependencies)
└── Manifest.toml           # Pinned dependency versions
```

### Core concepts in MicroSwimmers.jl

The package implements a regularised boundary element method for Stokes flow around microscale swimmers.

- **Models** (`PlanarFlagellum`, `StandingWaveFlagellum`, `EllipsoidBody`, …) — parametric descriptions of geometry and kinematics.
- **Discretised structures** (`Flagellum`, `CellBody`, `Flagellate`) — force-point / quadrature-point representations built from a model.
- **Problems** (`SwimmingProblem`, `ResistanceProblem`, `SwimmingTrajectoryProblem`, `ParticleTrajectoryProblem`) — set up and solve the linear system for forces, torques, and rigid-body velocities, or integrate the swimmer trajectory in time.
- **Fluid** (`PlanarVelocityField`, `TimeAveragedPlanarVelocityField`, `FluidVelocity`) — evaluate and average the induced flow field on a grid.

---

## Re-running simulations

The figure scripts use the pre-saved `results/*.jld2` files and do not re-run simulations. To reproduce the simulations from scratch, use the constructors in `BeatPlanes`:

```julia
using BeatPlanes
sim = ChlamySimulation("run1", beat_plane_tilt=(0.2, 0.0))
run!(sim)
chlamy = generate_chlamy(sim)
# ... see scripts/chlamy_analysis.jl for analysis functions

Numerical convergence of the BEM solver can be checked with:
```julia
include("scripts/convergence_checks.jl")
```

---

## Contact

James Cass — j.cass@exeter.ac.uk

