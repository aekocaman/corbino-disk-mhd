# corbino_turbulance_debug_hartmann — Hartmann-resolved laptop sibling case

A sibling of `corbino_turbulance_debug/`: same topology,  same
solver, same protocol  but only the z-grading is finer so that the Hartmann
boundary layer δ_H ≈ 53 μm is resolved with ≥ 5 cells.  Intended for an  laminar baseline + 
perturbation bracket.

> ** This is not the DNS mesh.**  r and θ resolutions are unchanged
> from the under-resolved debug case.  This mesh adds Hartmann-layer
> resolution only.

## Sizing

  - Cells per block: 200 × 96 × 70  (vs debug 200 × 96 × 30, vs production 600 × 384 × 70).
  - Total: 4 × 200 × 96 × 70 = **5.376 × 10⁶ cells** (~2.3× larger than debug).
  - Working set: ~8 GB at ~1.5 KB/cell 
  - Decomposition: 8 Scotch subdomains → ~0.67 M cells / rank.

## z-grading (Hartmann-resolved)

  - Wall layer: `Δz_min = 7 μm`, `q = 1.2`, 11 cells/side.
  - Side span = 7 μm × (1.2¹¹ − 1) / 0.2 ≈ 225 μm.
  - Cells inside δ_H = 53 μm: **8** (cells 1–8 cumulative widths 7…52 μm).
    "≥ 5 cells in δ_H" criterion satisfied.
  - Bulk: 48 uniform cells × 199 μm.
  - Total: 11 + 48 + 11 = 70.

## Compute budget (vs under-resolved debug)

| metric                       | debug (1 cell in δ_H) | hartmann (8 cells) |
|------------------------------|-----------------------|--------------------|
| total cells                  | 2.30 M                | 5.38 M             |
| Δt cap (CFL@dz_min, U=0.1)   | 5e-4 s                | 3.5e-5 s           |
| 0 → 0.5 s laminar baseline   | ~50 min               | ~7-9 h             |
| 50 ms perturbation restart   | ~3 min                | ~50 min            |
| 50 snapshot disk             | ~4.5 GB binary        | ~10.5 GB binary    |


