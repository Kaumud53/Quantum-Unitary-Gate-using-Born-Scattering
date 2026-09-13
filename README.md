# Born–Lippmann–Schwinger Framework for Quantum Gate Synthesis via Scattering Operators

**Authors:** Kumar Gautam, Harish Parthasarathy, Kaumud Sharma, Namisha Gupta, Divya Punia, Ajay K. Sharma

## Overview

This repository contains the manuscript, numerical code, and figures for a paper that develops a rigorous scattering-theoretic framework for **quantum gate synthesis**: instead of implementing a gate through time-dependent control fields or a sequence of elementary gates on a fixed platform, the gate is encoded **statically and permanently** in the spatial geometry of a scattering potential.

A particle prepared in a fixed-energy plane-wave state scatters off a static potential V(r); the energy-shell scattering matrix S(E), restricted to a discrete set of incident directions on the momentum sphere S², is interpreted directly as a finite-dimensional unitary quantum gate. Different logical basis states correspond to different incident directions rather than different internal quantum numbers.

## Key Contributions

- A complete, first-principles derivation of the **Born series expansion** of the scattered state, its **exact resolvent (operator) resummation**, and the **Møller wave operators** Ω±(E)
- Proof that the resulting **energy-shell scattering matrix** S(E) = Ω₊(E)*Ω₋(E) is unitary on L²(S²) — the key property motivating its interpretation as a quantum gate
- Explicit **plane-wave matrix elements** ⟨ψᵢ₂|S(E)|ψᵢ₁⟩ to first order in the coupling, derived both in the far-field (radiation-zone) limit and without that approximation
- Formulation of the **inverse scattering / gate-design problem**: choosing a potential V(r), incident directions, and energy so that the restricted scattering matrix approximates a target gate
- A **multi-channel optimization framework** (adapting a perturbed-Hamiltonian, Lagrange-multiplier method from prior work on driven harmonic oscillators) that reduces gate design to an exactly solvable, finite-dimensional Schrödinger evolution on the energy-degenerate truncated subspace
- A **numerical realization of the Hadamard gate** using a two-direction (single-qubit) encoding and a three-channel Gaussian potential (two isotropic + one dipole shape), demonstrating that the noise-to-signal ratio can be driven to machine-precision zero given sufficient design time — with **exact** (not merely perturbative) unitarity at every design time, since the direction-truncated Hilbert space is finite-dimensional

## Repository Contents

```
.
├── paper/              # Manuscript source (LaTeX, IEEE format) and compiled PDF
├── simulations/        # Numerical code for the multi-channel gate-design optimization and NSR(T) computation
├── figures/            # Generated figures (scattering-based gate schematic, mathematical pipeline diagram, NSR(T) vs T)
└── README.md
```

*(Adjust the folder names above to match the actual repository layout.)*

## Numerical Example: Single-Qubit Hadamard Gate

The paper's numerical demonstration (Section IV) uses:

| Quantity | Value |
|---|---|
| Incident directions | n̂₁ = (0,0,1), n̂₂ = (1,0,0) |
| Wavenumber | k = 1 |
| Channel potentials | 2 isotropic Gaussians (a₁ = 1, a₂ = 0.5) + 1 dipole (p-wave) Gaussian |
| Target gate | Hadamard, U_d = (1/√2)[[1,1],[1,−1]] |
| Scattering-strength budget | E₀ = 3.0 |
| Time discretization | K = 6 bins, optimized with L-BFGS + random restarts |

The two isotropic channels independently control the identity and σₓ components of the effective Hamiltonian; the dipole channel independently controls the σ_y component. Since [σₓ, σ_y] = 2iσ_z closes the Lie algebra to the full u(2), every single-qubit gate is in principle reachable given sufficient design time — confirmed numerically by NSR(T) decreasing monotonically to exactly zero for design time T ≳ 2.2.

### Reproducing the Results

```bash
# example — update to match actual script names
pip install numpy scipy matplotlib
python simulations/optimize_hadamard_gate.py
python simulations/plot_nsr_vs_T.py
```

This regenerates the NSR(T)-vs-design-time figure (Fig. 3) reported in the paper.

## Scope and Open Questions

The paper is explicit about what remains open for future work:

- The numerical demonstration uses the exact, non-perturbative, finite-time evolution operator U_N(T) rather than the strictly first-order stationary matrix from the fixed-energy S(E); making the adiabatic-switching limit U_N(T) → U_ℓℓ′ as T → ∞ quantitative is left for future work.
- The example is limited to N = 2 directions (a single qubit); scaling the multi-channel design procedure to larger direction sets (e.g. N = 8 for a three-qubit register) and to genuinely non-separable gates (e.g. controlled-unitary gates) remains to be tested.
- The box-normalized direction basis used for numerical simulation regularizes a formally divergent term in the plane-wave matrix elements; a fully rigorous continuum-to-discrete treatment as the number and density of encoded directions grows is not yet carried out.

## Citation

If you use this work, please cite:

```bibtex
@article{gautam_et_al_born_scattering_gates,
  title   = {Born--Lippmann--Schwinger Framework for Quantum Gate Synthesis via Scattering Operators},
  author  = {Gautam, Kumar and Parthasarathy, Harish and Sharma, Kaumud and Gupta, Namisha and Punia, Divya and Sharma, Ajay K.},
}
```

*(Update with the final venue, volume, page numbers, and DOI once published.)*

## Author Contributions

K.G. and H.P. wrote the main manuscript; K.S. and N.G. developed the mathematical framework; D.P. and A.S. provided overall guidance and supervision.

## Funding & Data Availability

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors. No new data were created or analyzed in this study.
