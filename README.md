# Microswimmers

Systematic parameter analysis of **Jakoba** microswimmers using the
[MicroSwimmers.jl](https://github.com/micromotility-lab/MicroSwimmers.jl).

Our analysis focuses on the Jakoba microorganism and studies its change in swimming performance and feeding efficiency.

---

## Features

- Parameter sweeps over flagellar movement paramteres.
- Swimming velocity analysis.
- Feeding feeding efficiency analysis.
- Power consumption and efficiency measurements.
- Heatmap visualisations using GLMakie.
- Helical trajectory analysis.

---

## Example Outputs

Some analyses performed by this repository include:

- Swimming velocity heatmaps across parameter spaces.
- Feeding efficiency plots.
- Helical trajectory reconstruction.

## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/Microswimmers.git
cd Microswimmers
```

Open Julia and install the required packages.

```julia
using Pkg

Pkg.add([
    "GLMakie",
    "FastGaussQuadrature",
    "Statistics"
])

Pkg.add(url="https://github.com/micromotility-lab/MicroSwimmer.jl")
Pkg.add(url="https://github.com/micromotility-lab/MicroSwimmersPlots.jl")
```

---

## Dependencies

| Package | Purpose |
|---------|---------|
| [MicroSwimmers.jl](https://github.com/micromotility-lab/MicroSwimmer.jl) | Microswimmer simulation framework |
| [MicroSwimmersPlots.jl](https://github.com/micromotility-lab/MicroSwimmersPlots.jl) | Visualisation utilities |
| `GLMakie` | Interactive plotting |
| `FastGaussQuadrature` | Numerical integration |
| `Statistics` | Statistical analysis |

---
