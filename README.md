# Born-Lippmann-Schwinger Framework for Quantum Gate Synthesis via Scattering Operators

**Authors:** Kumar Gautam, Harish Parthasarathy, Kaumud Sharma, Namisha Gupta, Divya Punia, Ajay K. Sharma
**Affiliations:** Department of Quantum Computing, QRACE, New Delhi · ECE Division, Netaji Subhas University of Technology (NSUT), Dwarka · Department of Computer Science and Engineering, NIT Delhi · Department of Electronics and Communication Engineering, NIT Delhi

> Status: Under review

## Abstract

This paper investigates how fixed-energy Born scattering theory can be used to realize quantum gates. The energy-shell scattering matrix S(E) is interpreted as a finite-dimensional unitary quantum gate by restricting it to a discrete computational basis of plane-wave momentum states, connecting quantum computation to scattering theory. We establish the unitarity of S(E) on the fixed-energy Hilbert space L²(S²), starting from the Lippmann-Schwinger formalism, developing the full Born-series expansion, deriving its exact resolvent representation, and constructing the Møller wave operators. We further formulate the **inverse scattering problem** of designing a potential that realizes a prescribed quantum gate and develop a multi-channel optimization framework for potential design. A numerical realization of the Hadamard gate using a Gaussian scattering potential demonstrates the feasibility of the framework.

## Core Idea

Instead of encoding a quantum gate as a time-dependent control sequence applied to internal degrees of freedom (spin, energy levels, or a circuit of elementary gates), this framework encodes the gate **permanently and statically** in the spatial geometry of a scattering potential *V(r)*. A particle prepared in a fixed-energy plane-wave state, incident along direction *n̂ℓ*, scatters off *V(r)* and emerges as a superposition over a discrete set of outgoing directions *{n̂₁, …, n̂_N}* on the momentum sphere S². The matrix of amplitudes

```
U_ℓℓ' = ⟨ψ_i,ℓ' | S(E) | ψ_i,ℓ⟩
```

is the realized N-dimensional unitary gate — fixed once and for all by the shape of *V(r)*, not switched on/off in time.

## Key Contributions

1. **Full Born-series expansion** of the scattered state, resummed into the exact resolvent (operator) form |ψ_f⟩ = [I − ε(H₁ − E)⁻¹V]|ψ_i⟩ (Sec. II-A–C).
2. **Møller wave operators** Ω±(E) and the energy-shell **scattering matrix** S(E) = Ω₊(E)*Ω₋(E), proven unitary on L²(S²) (Sec. II-D–E).
3. **Plane-wave matrix elements** ⟨ψ_i2|S(E)|ψ_i1⟩ derived to first order in the coupling ε, both in the far-field (radiation-zone) limit and without that approximation (Sec. II-F).
4. **Multi-channel potential-design framework** (Sec. III): the scattering potential is written as a linear combination of M fixed spatial shape functions with time-dependent design amplitudes εₐ(t); projecting onto a truncated N-direction basis gives N×N Hermitian channel matrices V₁,…,V_M and an exactly-solvable, exactly-unitary finite-dimensional time-dependent Hamiltonian H_N(t).
5. **Gate-error energy / scattering-strength budget optimization** (Lagrange-multiplier stationarity condition, Eq. 52), adapting the perturbed-Hamiltonian method of Gautam et al. (2015) from a driven harmonic oscillator to the scattering setting.
6. **Numerical realization of a single-qubit Hadamard gate** (Sec. IV): a two-direction encoding (N=2) with a three-channel Gaussian potential (two isotropic s-wave shapes + one dipole p-wave shape) reaches the target gate to machine precision (NSR → 0) for sufficient design time T, and is *exactly* unitary at every T (not just perturbatively), since the truncated Hilbert space is finite-dimensional.

## Repository Contents

```
.
├── paper/                     # Manuscript source (LaTeX) and compiled PDF
├── src/
│   ├── channel_matrices.py     # Gaussian shape functions and their Fourier transforms (Eq. 41)
│   ├── hadamard_synthesis.py   # Multi-channel Lagrange/L-BFGS gate-synthesis optimizer (Sec. III-IV)
│   └── plot_nsr_vs_T.py        # Reproduces Fig. 3 (NSR(T) convergence curve)
├── figures/
│   └── fig3_nsr_vs_T.png
└── README.md
```

*(Adjust the tree above to match your actual repo layout before pushing.)*

## Numerical Example: Hadamard Gate

| Parameter | Value |
|---|---|
| Wavenumber k | 1 |
| Incident directions | n̂₁ = (0,0,1), n̂₂ = (1,0,0) |
| Channel shapes | 2 isotropic Gaussians (a₁=1, a₂=0.5) + 1 dipole Gaussian |
| Channel matrices | V₁ₐ = 15.75 I₂ + 5.79 σₓ, V₁ᵦ = 1.969 I₂ + 1.533 σₓ, V₂ = −5.79 σ_y |
| Target gate | Hadamard, U_d = (1/√2)[[1,1],[1,−1]] |
| Time bins K | 6 |
| Scattering-strength budget E₀ | 3.0 |
| Penalty weight Γ | 80 |
| Optimizer | L-BFGS with multiple random restarts |

**Result:** NSR(T) ≈ 2 for short design times (T ≲ 0.3, gate ≈ identity), decreasing smoothly and monotonically to **machine-precision zero for T ≳ 2.2**. Because the two-direction truncation is exactly (not perturbatively) unitary, NSR(T) reaches exactly zero rather than saturating at a residual floor, unlike the infinite-dimensional driven-oscillator case of Gautam (2015).

## Code

The script below reproduces the Sec. IV numerical example: it builds the calibrated channel Hamiltonian H₂(t) = E·I₂ + c₀(t)I₂ + c_x(t)σₓ + c_y(t)σ_y, discretizes the design window [0,T] into K=6 piecewise-constant bins, forms the exact time-ordered propagator U₂(T) as a product of matrix exponentials (exactly unitary at every step — no Dyson-series truncation), and minimizes NSR(T) + Γ(E_diss − E₀)² over the bin amplitudes via L-BFGS with random restarts, sweeping T to reproduce Fig. 3.

```python
"""
hadamard_synthesis.py

Numerical realization of a single-qubit Hadamard gate via multi-channel
Born-scattering potential design (Sec. III-IV of the paper).

Reproduces Fig. 3: NSR(T) vs design time T.
"""

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

# ----------------------------------------------------------------------
# Pauli matrices
# ----------------------------------------------------------------------
I2 = np.eye(2, dtype=complex)
sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma_z = np.array([[1, 0], [0, -1]], dtype=complex)

# ----------------------------------------------------------------------
# Target gate: single-qubit Hadamard (Eq. 58)
# ----------------------------------------------------------------------
U_d = (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)

# ----------------------------------------------------------------------
# Physical setup (Sec. IV-A)
# ----------------------------------------------------------------------
k = 1.0
E = k**2 / 2.0  # collision energy (hbar = m = 1 units), E = k^2/2 (non-relativistic)

# Channel matrices calibrated in the {|0>,|1>} basis (Eqs. 54-57)
# H2(t) = E*I2 + c0(t)*I2 + cx(t)*sigma_x + cy(t)*sigma_y
# We optimize directly over the calibrated triad (c0, cx, cy) per Sec. IV-A.


def time_ordered_propagator(c0, cx, cy, T, K):
    """
    Build the exact time-ordered product U2(T) = prod_{j=K}^{1} exp(-i*(T/K)*H_j)
    (Eq. 59). Each bin's Hamiltonian is exactly Hermitian, so each factor is
    exactly unitary -- no Dyson/Born truncation is used here.

    c0, cx, cy : arrays of length K (piecewise-constant bin amplitudes)
    T          : total design time
    K          : number of time bins
    """
    dt = T / K
    U = np.eye(2, dtype=complex)
    for j in range(K):
        Hj = E * I2 + c0[j] * I2 + cx[j] * sigma_x + cy[j] * sigma_y
        Uj = expm(-1j * dt * Hj)
        U = Uj @ U  # left-multiply: last bin applied last
    return U


def gate_error_energy(U2_T, T):
    """
    E(T) = || U_d' - U2(T) ||_F^2,  U_d' = exp(i*E*T) * U_d   (Eq. 47-48)
    """
    U_d_prime = np.exp(1j * E * T) * U_d
    diff = U_d_prime - U2_T
    return np.real(np.trace(diff.conj().T @ diff))


def dissipation_energy(c0, cx, cy, T, K):
    """
    E_diss = sum_a integral_0^T eps_a(t)^2 dt, approximated on the piecewise-
    constant discretization as sum_j (c_j^2) * dt   (Eq. 49).
    """
    dt = T / K
    return dt * np.sum(c0**2 + cx**2 + cy**2)


def objective(params, T, K, E0, Gamma):
    """
    Minimize NSR(T) + Gamma*(E_diss - E0)^2   (Sec. IV-B)
    NSR(T) = E(T) / ||U_d||_F^2   (Eq. 50); ||U_d||_F^2 = 2 for a unitary 2x2 gate.
    """
    c0, cx, cy = np.split(params, 3)
    U2_T = time_ordered_propagator(c0, cx, cy, T, K)
    E_T = gate_error_energy(U2_T, T)
    NSR_T = E_T / 2.0  # ||U_d||_F^2 = Tr[U_d^dagger U_d] = 2
    E_diss = dissipation_energy(c0, cx, cy, T, K)
    return NSR_T + Gamma * (E_diss - E0) ** 2


def optimize_gate(T, K=6, E0=3.0, Gamma=80.0, n_restarts=10, seed=0):
    """
    Minimize the objective over the 3K real design parameters using L-BFGS
    with multiple random restarts, keeping the best-found solution.
    """
    rng = np.random.default_rng(seed)
    best_val = np.inf
    best_params = None
    for _ in range(n_restarts):
        x0 = rng.normal(scale=1.0, size=3 * K)
        res = minimize(objective, x0, args=(T, K, E0, Gamma), method="L-BFGS-B")
        if res.fun < best_val:
            best_val = res.fun
            best_params = res.x
    c0, cx, cy = np.split(best_params, 3)
    U2_T = time_ordered_propagator(c0, cx, cy, T, K)
    NSR_T = gate_error_energy(U2_T, T) / 2.0
    return NSR_T, best_params


def sweep_design_time(T_values, K=6, E0=3.0, Gamma=80.0, n_restarts=10, seed=0):
    """
    Sweep design time T and record the best-achieved NSR(T), reproducing Fig. 3.
    """
    nsr_values = []
    for T in T_values:
        nsr_T, _ = optimize_gate(T, K=K, E0=E0, Gamma=Gamma,
                                  n_restarts=n_restarts, seed=seed)
        nsr_values.append(nsr_T)
        print(f"T = {T:5.2f}   NSR(T) = {nsr_T:.6e}")
    return np.array(nsr_values)


if __name__ == "__main__":
    T_values = np.linspace(0.05, 8.0, 40)
    nsr_curve = sweep_design_time(T_values)

    # Optional plotting (matplotlib)
    try:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 4))
        plt.plot(T_values, nsr_curve, "o-", lw=1.5, ms=4)
        plt.xlabel("design time T")
        plt.ylabel("NSR(T)")
        plt.title("Noise-to-signal ratio vs design time (Hadamard gate synthesis)")
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig("fig3_nsr_vs_T.png", dpi=150)
        print("Saved figures/fig3_nsr_vs_T.png")
    except ImportError:
        pass
```

Run it directly:

```bash
pip install numpy scipy matplotlib
python src/hadamard_synthesis.py
```

This sweeps design time T from ~0 to 8, optimizing the 18 real design parameters (K=6 bins × 3 channels) at each T via L-BFGS with random restarts, and reports NSR(T) — expect it to fall from ≈2 (near-identity gate) toward ≈0 (exact Hadamard) around T ≈ 2–2.5, matching Fig. 3 of the paper.

> **Note:** the channel matrices V₁ₐ, V₁ᵦ, V₂ (Eqs. 54–56) in the paper are derived from the first-order Born matrix element (Eq. 41) applied to specific Gaussian potential shapes; the script above works directly with the already-calibrated (c₀, c_x, c_y) triad as described in Sec. IV-A, which is the quantity actually optimized in the paper's numerical example. A `channel_matrices.py` module deriving V₁ₐ, V₁ᵦ, V₂ from the raw Gaussian shape functions v₁ₐ(r), v₁ᵦ(r), v₂(r) via Eq. (41) can be added if you want the full potential-to-gate pipeline rather than starting from the calibrated triad.

## Open Questions / Future Work (Sec. V)

- Making the adiabatic-switching argument (Sec. III-B) quantitative: establishing the rate at which U_N(T) → U_ℓℓ' as T → ∞, and relating it to the higher-order Born corrections needed for the strictly stationary S(E) to be unitary beyond O(ε²).
- Scaling the multi-channel design procedure beyond N=2 directions (e.g., N=8 for a three-qubit register) and to non-separable gates (e.g., a controlled-unitary gate), testing the scalability of the Lagrange-multiplier stationarity condition.
- A fully rigorous treatment of the continuum-to-discrete limit as the number and density of encoded directions grows (the box-normalized direction basis used here regularizes the formally divergent δ(n̂₂ − n̂₁) term only for numerical purposes).

## Citation

```bibtex
@article{gautam_born_scattering_gates,
  title   = {Born-Lippmann-Schwinger Framework for Quantum Gate Synthesis via Scattering Operators},
  author  = {Gautam, Kumar and Parthasarathy, Harish and Sharma, Kaumud and Gupta, Namisha and Punia, Divya and Sharma, Ajay K.},
  journal = {Under review},
  year    = {2026}
}
```

## Related Work

- K. Gautam, "Study and implementation of unitary gates in quantum computation using Schrödinger dynamics," *Quantum Information Processing*, 2015. (The perturbed-Hamiltonian, Lagrange-multiplier optimization method adapted here from a driven harmonic oscillator.)
- K. Gautam et al., "Electrodynamics-based quantum gate optimization with Born scattering," *Scientific Reports*, vol. 14, p. 27838, 2024.

## Author Contributions

K.G. and H.P. wrote the main manuscript; K.S. and N.G. developed the mathematical framework; D.P. and A.K.S. provided overall guidance and supervision.

## Funding

This research received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors.
