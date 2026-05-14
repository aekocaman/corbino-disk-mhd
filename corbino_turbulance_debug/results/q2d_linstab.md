# Q2D Linear Stability of the Laminar Corbino Flow

Generalized eigenvalue solution of the Q2D Sommeria-Moreau linearised system around the laminar Corbino base $\langle U_{\theta,0}\rangle_z(r) = K/r$. Discretization: Chebyshev collocation on $N+1 = 33$ nodes; clamped streamfunction BCs ($\psi = \psi' = 0$) imposed by null-space projection of the 4-row constraint matrix.

- Hartmann number:       $Ha = 94.88$
- Hartmann time:         $\tau_H = 706.40$ ms
- Analytical prediction: $\sigma\tau_A \approx -1/Ha = -0.01054$

## Computed `sigma_max * tau_A`

| $Re_a/Ha^2$ | $U_\infty$ (cm/s) | $\sigma_{\max}\tau_A$ (max over $m$) | $|\sigma_{\max}\tau_A - (-1/Ha)|$ |
|------------:|------------------:|-----------------------------------:|----------------------------------:|
| 0.50 | 2.099 | $-0.01070$ | $1.65e-04$ |
| 0.75 | 3.148 | $-0.01074$ | $2.04e-04$ |
| 1.00 | 4.197 | $-0.01024$ | $3.01e-04$ |
| 1.19 | 4.995 | $-0.00978$ | $7.63e-04$ |
| 1.50 | 6.296 | $-0.00879$ | $1.75e-03$ |
| 2.00 | 8.395 | $-0.00723$ | $3.31e-03$ |
| 3.00 | 12.592 | $-0.00339$ | $7.15e-03$ |
| 4.00 | 16.790 | $+0.00096$ | $1.15e-02$ |
| 5.00 | 20.987 | $+0.00604$ | $1.66e-02$ |
| 7.00 | 29.382 | $+0.01570$ | $2.62e-02$ |
| 10.00 | 41.974 | $+0.02953$ | $4.01e-02$ |

## Verdict

At low $Re_a/Ha^2$ the dominant decay rate matches the Hartmann-friction scale $-1/Ha = -0.01054$, as expected from the Hermitian part of the linearised operator. As $Re_a/Ha^2$ increases, the imaginary advection term $-imK/r^2$ does **not** commute with $L_m$, and the non-normal spectral shift moves $\Re(\sigma_{\max})$ monotonically toward zero. By extrapolation of the trend across the table above, the Q2D marginal stability condition $\sigma_{\max} = 0$ is approached at a moderate $Re_a/Ha^2$ above the experimentally-observed transition value of $\sim 1$. Q2D therefore does carry non-trivial linear-stability information, but the threshold it reports does not coincide with the PK14 transition; the discrepancy is consistent with the Pothérat-Klein interpretation that the actual transition is a 3D phenomenon.

The full 3D linear-stability problem in cylindrical coordinates -- which couples the Hartmann-resolved $z$ structure to the in-plane modal structure -- is computationally non-trivial and is deferred to the cluster runs of the production package.
