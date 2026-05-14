# Mesh-sensitivity comparison: debug vs Hartmann-resolved

Three-panel mesh-sensitivity comparison.

**(a)** $\sigma\tau_A$ vs probe-RMS amplitude. The five filled circles (viridis-coloured by $R^2$) are the debug-mesh same-seed bracket, $\alpha\in[2,30]\%$. The grey shaded band is the linear extrapolation of the debug trend with a $\pm\sigma_{\rm seed}\tau_A = \pm 0.10$ uncertainty estimate from the original independent-realisation bracket. Solid red star: the single Hartmann-resolved data point at its effective amplitude ($\sim 85\%$ at probes); $\sigma\tau_A = -0.227$, $R^2 = 0.98$. Open red star: the bias-corrected position the same data point would occupy at its requested 5% amplitude if `perturb_ic.py` used volume-weighted RMS rescaling. The dotted arrow links the two positions, indicating the magnitude of the amplitude bias caused by the cell-weighting bug. **(b)** Single-exponential fit quality $R^2$ vs amplitude. The Hartmann-resolved value of 0.978 is above the debug-mesh trend at all amplitudes, consistent with the cleaner exponential expected at large effective amplitude. **(c)** Back-EMF efficiency $\langle\Phi'\rangle_{\rm RMS}/(U'_z B L_{\rm int})$ vs amplitude. The debug-mesh values (blue, three points) are amplitude-independent at $\sim 0.27$, the signature of linear electromagnetic coupling. The Hartmann-resolved value is $\sim 14$ (off-scale; box at lower right), reflecting the same cell-weighting bias that inflated the panel-(a) amplitude: a $5\times$ over-estimate in the perturbation magnitude $U'_z$ used as the denominator of the efficiency.

## Data table

| mesh | amplitude | $\sigma\tau_A$ | $R^{2}$ | e-fold (ms) |
|------|-----------|-----------------|---------|-------------|
| debug | 2\% | $-0.298$ | 0.81 | 24.9 |
| debug | 5\% | $-0.310$ | 0.84 | 23.9 |
| debug | 10\% | $-0.332$ | 0.87 | 22.3 |
| debug | 15\% | $-0.353$ | 0.91 | 21.0 |
| debug | 30\% | $-0.404$ | 0.95 | 18.3 |
| Hartmann-resolved | 85\% effective (5\% requested) | $-0.227$ | 0.98 | 32.6 |
