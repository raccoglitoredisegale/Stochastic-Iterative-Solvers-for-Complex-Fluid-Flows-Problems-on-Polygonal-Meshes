# Stochastic Iterative Solvers for Complex Fluid Flows Problems on Polygonal Meshes

**Authors:** Christian Capelli and Matteo Morini  
**Programme:** MSc in Mathematical Engineering — Numerical Analysis for Partial Differential Equations  
**Academic year:** 2025/2026  
**Tutors:** Paola Francesca Antonietti, Gabriele Ciaramella, and Alessandra Cancrini

## Overview

This project investigates iterative solvers for the **pseudo-stress formulation of the unsteady Stokes problem**, which models incompressible viscous flow. The equations are discretized using a **polytopal discontinuous Galerkin method (PolyDG)** in space and **Crank–Nicolson** time integration.

The central objective is to assess whether Kaczmarz-type inner iterations can provide effective preconditioning for **flexible AB-GMRES**. The report compares classical, randomized, and greedy row-selection strategies, identifies their limitations on PolyDG systems, and introduces an adaptive solver to reduce unproductive inner work.

## Mathematical model and discretization

The primary unknown is the pseudo-stress tensor

$$
\sigma = \mu\nabla u - pI,
$$

where $u$ is velocity, $p$ is pressure, and $\mu$ is viscosity. PolyDG approximates this tensor with discontinuous piecewise polynomials on polygonal elements, using interface terms and penalties to couple neighboring elements.

The spatial discretization produces

$$
M\dot{\sigma}_h + A\sigma_h = f.
$$

Crank–Nicolson then requires solving the following system at each time step:

$$
\left(M+\frac{\Delta t}{2}A\right)\sigma_h^{n+1}
=\left(M-\frac{\Delta t}{2}A\right)\sigma_h^n
+\frac{\Delta t}{2}\left(f^{n+1}+f^n\right).
$$

The expected energy-norm discretization error is $O(\Delta t^2)+O(h^q)$, with mesh size $h$ and polynomial degree $q$. Here, $q$ distinguishes the approximation degree from pressure $p$. At fixed spatial resolution, reducing the time step eventually reaches a spatial-error plateau. Very small time steps can also worsen conditioning as the system increasingly reflects the singular mass operator.

## Solver approach

Flexible AB-GMRES uses an inner Kaczmarz solver to generate a preconditioned search direction at each outer iteration. The inner solve stops when its relative residual reaches a threshold $\eta$, or when it reaches an iteration cap. Flexibility allows the preconditioning action and amount of inner work to vary between outer iterations.

| Method | Row-selection strategy |
| --- | --- |
| Classical Kaczmarz (K) | Visits rows cyclically. |
| Randomized Kaczmarz (RK) | Samples rows with probabilities proportional to their squared norms. |
| Greedy Kaczmarz (GK) | Selects the row with the largest absolute residual component. |
| Greedy Randomized Kaczmarz (GRK) | Randomly selects among rows with sufficiently large residual contributions. |
| NE-SOR baseline | Uses relaxed cyclic row sweeps within AB-GMRES. |

A parameter-tuning procedure selects the relaxation parameter $\omega$ and the maximum number of inner iterations. Smaller $\eta$ values demand more accurate, but more expensive, inner solves.

## Numerical experiments and findings

The implementation uses **MATLAB** and the **lymph** library. Experiments cover two settings:

- **Synthetic linear systems:** underdetermined and overdetermined matrices of sizes $500\times5000$ and $5000\times500$, with 20% density and condition numbers approximately $10$ and $10^6$.
- **PolyDG flow systems:** a manufactured pseudo-stress solution on the unit square, with mixed boundary conditions, polygonal meshes, and varying time steps. Computed solutions are compared with the exact solution and MATLAB's backslash solver.

On the underdetermined ill-conditioned benchmark, GK gives the lowest reported CPU times, while GRK can require fewer outer iterations at tighter inner tolerances. These advantages are not uniform across every benchmark: NE-SOR remains competitive and has lower reported CPU times in the overdetermined RANDL6 tests. RK generally requires substantially more inner work.

On PolyDG systems, the inner residual often stagnates before reaching its target, causing repeated use of the full iteration budget. The report attributes this behavior to row-wise updates inadequately exploiting the matrix's block structure. Too little inner work can also produce ineffective search directions, so simply lowering the cap is insufficient.

## Adaptive flexible AB-GMRES

The proposed adaptive strategy adjusts the inner iteration budget using residual behavior and increases the cap when additional work is needed. It also uses restarts and checks the true outer residual to continue toward the requested tolerance.

For the reported five-step comparison with $\Delta t=0.05$ and outer tolerance $10^{-3}$, the adaptive solver starts with an inner cap of 2,048 and a growth factor of 2; the standard solver uses a fixed cap of 60,000.

| Metric | Adaptive solver | Standard solver |
| --- | ---: | ---: |
| Total outer iterations | 3,600 | 413 |
| Total inner iterations | 12,224,948 | 24,515,059 |
| Total CPU time | 719.93 s | 1,014.75 s |
| Reported maximum error | 0.07581016 | 0.07614793 |

The adaptive method reduces CPU time by approximately **29%** and inner iterations by approximately **50%**, while maintaining a similar maximum error. It achieves this with more outer iterations. These figures describe the specific experiment in Table 8 of the report, rather than a general performance guarantee.

## Conclusions and future work

Kaczmarz-type inner preconditioners can be effective on standard linear systems, particularly with greedy row selection, but their performance depends strongly on parameter choices and matrix structure. Their direct row-wise application to PolyDG systems remains computationally expensive and is not yet scalable.

Adaptivity partially mitigates stagnation. The main proposed extension is **Block Kaczmarz**, which would update groups of coupled rows to exploit the PolyDG matrix structure. This extension is proposed as future work and was not implemented in the report.

## Report guide

- **Section 2:** AB-GMRES, Kaczmarz variants, parameter tuning, and benchmark tests.
- **Section 3:** Pseudo-stress formulation and PolyDG–Crank–Nicolson discretization.
- **Section 4:** Flow simulations, stagnation analysis, and adaptive solver results.
- **Section 5:** Conclusions and proposed extensions.
