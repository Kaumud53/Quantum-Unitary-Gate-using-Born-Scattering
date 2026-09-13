# A Rigorous Quantum Communication Framework for Optical Fibre Channels

**Integrated Control, Computational Hardness, and Noise-Assisted Security**

**Authors:** K. Gautam, Kaumud Sharma

## Overview

This repository contains the manuscript, simulation code, and figures for a paper that develops a unified, five-layer mathematical framework for quantum communication over optical fibres — connecting the physics of electromagnetic propagation directly to information-theoretic, computational, and cryptographic limits.

The five interlocking layers are:

1. **Quantum field quantisation** over modal decompositions of the fibre (retarded-potential formulation, atom–field interaction Hamiltonian)
2. **Classical–Quantum (Cq) channel theory** with Holevo capacity optimisation under realistic fibre parameters (BPSK coherent-state encoding)
3. **GKSL noise modelling** with exact Kraus-operator solutions and the Quantum Data Processing Inequality
4. **Integrated quantum control** applied *during* channel evolution, proven to strictly increase capacity beyond the Data Processing Inequality ceiling
5. **Multi-parameter quantum channel estimation** via the Quantum Fisher Information Matrix (QFIM) and Quantum Cramér–Rao Bound (QCRB)

## Key Contributions

- A complete derivation path from Maxwell's equations (retarded potentials, modal decomposition, atom–field Hamiltonian) to the Holevo channel capacity, with all fibre parameters (core radius, attenuation, dispersion, temperature) made explicit
- A proof (Theorem IV.1, IV.2) that a control Hamiltonian acting *concurrently* with channel dissipation can strictly increase Cq channel capacity — a regime the Quantum Data Processing Inequality does not constrain, since that inequality only bounds post-channel CPTP processing
- Introduction of the **Semigroup Julia Inversion Problem (SJIP)** and a rigorous, unconditional polynomial-time reduction from Subset-Sum proving SJIP is **NP-hard** (Theorem IV.3)
- A carefully flagged **conjecture** (not a theorem) that optimal coherent-state ML decoding is at least as hard as SJIP, together with an explicit account of exactly which step in that correspondence remains unproven
- Closed-form QFIM/QCRB derivations for a combined loss-plus-dephasing bosonic channel, with an adaptive Bayesian estimation protocol (Algorithm 1) shown to be QCRB-saturating
- A **dual-layer security bound**: an unconditional, exponentially-decaying information-theoretic term (from channel non-injectivity) plus a computational term that is negligible *conditional on* the SJIP-decoding conjecture
- Four explicit, experimentally testable predictions relating capacity to core radius, propagation length, and temperature

The paper is explicit throughout about which results are unconditional theorems versus the single conjectural step (Conjecture IV.5), and deliberately avoids overclaiming the cryptographic guarantee.

## Repository Contents

```
.
├── paper/              # Manuscript source (LaTeX) and compiled PDF
├── simulations/        # Python implementation of Algorithm 1 (adaptive Bayesian QFIM estimator)
├── figures/            # Generated figures (framework diagram, posterior convergence, QCRB comparison)
└── README.md
```

*(Adjust the folder names above to match the actual repository layout.)*

## Simulation Details

Numerical validation of Algorithm 1 (adaptive Bayesian channel parameter estimation) was implemented in **Python** using a sequential Monte Carlo / particle-filter approximation.

| Quantity | Value |
|---|---|
| Loss rate, γ_loss | 0.02 ns⁻¹ |
| Dephasing rate, γ_φ | 0.005 ns⁻¹ |
| Propagation length, L | 10 km |
| Particles, N_p | 300 |
| True parameters, (n̄, φ) | (5.0, 0.30 rad) |
| Prior mean, (n̄₀, φ₀) | (4.0, 0.50 rad) (deliberately misspecified) |
| Measurement rounds, M | 50 |

The simulation confirms consistent parameter recovery from a misspecified prior, and shows the posterior variance dropping below the frequentist QCRB threshold — consistent with Bayesian efficiency, since the prior supplies information beyond the channel measurements themselves (see Remark VI.1 in the paper).

### Reproducing the Results

```bash
# example — update to match actual script names
pip install numpy scipy matplotlib
python simulations/run_algorithm1.py
```

This regenerates the posterior mean convergence plots and the posterior-variance-vs-QCRB comparison reported in the paper.

## Scope and Honesty Notes

- The explicit worked examples use **BPSK** coherent-state encoding; the general Cq-channel, GKSL, control, and QFIM machinery is stated for arbitrary alphabets, but extending the closed-form capacity/QFIM expressions to M-ary constellations is left for future work.
- The SJIP **NP-hardness proof (Theorem IV.3) is unconditional**. The link between SJIP and optimal quantum decoding (**Conjecture IV.5**) is explicitly *not* claimed as proven — the paper states precisely which algebraic correspondence would need to be established to upgrade it to a theorem.
- The security guarantee (Corollary V.9) therefore has one unconditional term (information-theoretic, from channel non-injectivity) and one conditional term (computational, contingent on Conjecture IV.5 and NP ⊄ BQP).

## Citation

If you use this work, please cite:

```bibtex
@article{gautam_sharma_cq_fibre_framework,
  title   = {A Rigorous Quantum Communication Framework for Optical Fibre Channels:
             Integrated Control, Computational Hardness, and Noise-Assisted Security},
  author  = {Gautam, K. and Sharma, Kaumud},
}
```

*(Update with the final venue, volume, page numbers, and DOI once published.)*

## Acknowledgements

The authors thank Prof. K. R. Parthasarathy for invaluable guidance on quantum stochastic differential equations.
