# Perturbation amplitude bracket — debug mesh, seed 1729

All runs use the divergence-free vortical perturbation of `perturb_ic.py` with integral scale 5 mm, applied to the laminar baseline at t = 0.5 s, integrated for 50 ms (6.76 τ_A).

| amp | RMS U_z init | RMS U_z final | ratio | σ × τ_A | R² | e-fold [ms] | Φ'/scale (med) | verdict |
|-----|--------------|---------------|-------|---------|-----|-------------|----------------|---------|
| 2 % | 3.19e-03 | 3.70e-04 | 0.116 | -0.298 | 0.814 | 24.9 | 0.494 | **DECAY** |
| 5 % | 7.97e-03 | 8.72e-04 | 0.109 | -0.310 | 0.837 | 23.9 | 0.270 | **DECAY** |
| 10 % | 1.59e-02 | 1.56e-03 | 0.098 | -0.332 | 0.874 | 22.3 | 0.279 | **DECAY** |
| 15 % | 2.39e-02 | 2.06e-03 | 0.086 | -0.353 | 0.905 | 21.0 | 0.291 | **DECAY** |
| 30 % | 4.78e-02 | 2.89e-03 | 0.060 | -0.404 | 0.951 | 18.3 | 0.315 | **DECAY** |

Reading the table: `RMS U_z` is the perturbation amplitude diagnostic (zero in the laminar Hartmann state by symmetry). `σ × τ_A` is the dimensionless decay rate -- negative means the perturbation decays. `Φ'/scale (med)` is the median ratio of the realised back-EMF to its analytic eddy-scale estimate U'_z·B·L_integral; constancy across amplitudes indicates linear EM coupling.
