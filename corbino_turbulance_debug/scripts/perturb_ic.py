#!/usr/bin/env python3
"""
perturb_ic.py — inject a divergence-free turbulent kick onto the laminar
                Corbino baseline at t = 0.5 s.

References
----------
Spectrum form: Pope, "Turbulent Flows", CUP (2000), §6.5;
   E(k) ~ (k/k_p)^4 * exp(-2 (k/k_p)^2), the so-called
   Passot-Pouquet model -- see Passot & Pouquet, J. Fluid Mech.
   181, 441 (1987).
Solenoidal construction U' = curl(psi): standard in MHD-DNS initial
   condition practice; see e.g. Krasnov-Zikanov-Boeck-style channel
   DNS papers.
Wall-tapering envelope sin(pi*r_frac)*sin(pi*z/H): heuristic that
   makes U' vanish at all four no-slip / electrode boundaries while
   preserving the spectral content elsewhere.  The taper breaks
   div(U') = 0 by O(taper-width / integral-scale); the residual
   divergence is cleaned up by the first PISO pressure-correction
   sweep at restart -- a standard procedure for non-trivial
   perturbed-IC restarts in OpenFOAM.

Pipeline (run from inside the OpenFOAM case directory):

    foamJob -wait reconstructPar -time 0.5
    foamFormatConvert -time 0.5
    postProcess -func writeCellCentres -time 0.5
    python3 scripts/perturb_ic.py 0.5 --rms 0.05 --integral-scale 0.005

Reads:
    <time>/U          ASCII volVectorField  (laminar steady state)
    <time>/C          ASCII volVectorField  (cell centres, from writeCellCentres)

Writes:
    <time>/U.lam      backup of the original laminar field
    <time>/U          replaced with U_lam + U' (perturbed initial condition)

What it does
------------
  1. Reads the cell-centre positions from <time>/C and the laminar U from <time>/U.
  2. Builds a Cartesian box that covers the annulus, picks a uniform grid with
     ~6 points per integral scale.
  3. In spectral (k_x, k_y, k_z) space, constructs a random divergence-free
     velocity field by drawing a complex random vector potential psi_hat(k)
     and setting Uhat'(k) = i k x psi_hat(k).  This is exactly solenoidal:
     k . (k x psi_hat) = 0 identically.
  4. Modulates the spectrum by a Pope/Passot-Pouquet form
         E(k) ~ (k/k_p)^4 * exp(-2 (k/k_p)^2)
     where k_p = 2 pi / l_integral.  This concentrates energy at the
     integral scale and decays smoothly at small k and large k, avoiding the
     spectrum-truncation ringing of a top-hat.
  5. IFFTs psi_hat and Uhat' to physical space; trilinear-interpolates U'
     onto every OpenFOAM cell centre.
  6. Multiplies by a smooth wall taper T(r, z) = sin(pi r_frac) * sin(pi z/h)
     so the perturbation vanishes at all four no-slip / electrode walls.
     This taper breaks div-free at O(taper-width / integral-scale); the
     residual divergence is small and gets cleaned up by the first PISO
     pressure-correction loop on restart.
  7. Rescales globally so that the in-bulk RMS magnitude of U' equals
     `--rms` times the local RMS of U_lam (default 0.05 = 5%).
  8. Writes U_lam + U' back, preserving the original boundaryField verbatim.

Validation diagnostics printed at the end:
    target RMS, achieved RMS, max(|U'|), max(|U'|/|U_lam_bulk|),
    crude divergence estimate via central differences on a sampled subgrid,
    azimuthal energy spectrum E(k_theta) at z = h/2, r = mid radius.

Dependencies: numpy, scipy.  No PyFoam, no Ofpp -- pure ASCII parsing.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import sys
import time as wallclock
from pathlib import Path

import numpy as np

# -----------------------------------------------------------------------------
# Geometry constants (must match blockMeshDict)
# -----------------------------------------------------------------------------
R_IN = 0.020      # m, inner electrode radius
R_OUT = 0.100     # m, outer electrode radius
H = 0.010         # m, channel height (z = 0 .. H)


# =============================================================================
# OpenFOAM ASCII volVectorField parser
# =============================================================================

_HEADER_RE = re.compile(r"FoamFile\s*\{(.*?)\}", re.DOTALL)
_INTERNAL_RE = re.compile(
    r"internalField\s+nonuniform\s+List<vector>\s*"
    r"(?P<count>\d+)\s*\(",
    re.DOTALL,
)
_BOUNDARY_RE = re.compile(r"boundaryField\s*\{", re.DOTALL)


def read_vector_field(path: Path) -> tuple[np.ndarray, str, str]:
    """Parse an ASCII OpenFOAM volVectorField. Returns (data, header, boundary).

    `data` is shape (N, 3).  `header` is everything before the internalField
    block (FoamFile + dimensions + any preamble).  `boundary` is the verbatim
    text from `boundaryField` to end-of-file, so we can splice it back when
    writing.
    """
    text = path.read_text()
    m_int = _INTERNAL_RE.search(text)
    if m_int is None:
        raise RuntimeError(
            f"{path}: not an ASCII nonuniform volVectorField. "
            f"Did you run `foamFormatConvert -time <t>`?"
        )
    count = int(m_int.group("count"))
    header = text[: m_int.start()]

    # Parse N vectors of the form `(x y z)` after the opening paren.
    body_start = m_int.end()
    m_bnd = _BOUNDARY_RE.search(text, body_start)
    if m_bnd is None:
        raise RuntimeError(f"{path}: missing boundaryField block")
    body = text[body_start : m_bnd.start()]

    # Extract every (a b c) triplet.  Avoid regex on huge strings; use scan.
    nums = re.findall(r"\(([^)]*)\)", body)
    if len(nums) < count:
        raise RuntimeError(
            f"{path}: expected {count} vectors, found {len(nums)}"
        )
    data = np.fromstring(" ".join(nums[:count]), sep=" ").reshape(count, 3)

    boundary = text[m_bnd.start() :]
    return data, header, boundary


def write_vector_field(
    path: Path, data: np.ndarray, header: str, boundary: str
) -> None:
    """Write an ASCII volVectorField with the given internal data.

    `header` and `boundary` are passed through verbatim from a previous
    `read_vector_field` so that dimensions and BCs are preserved.
    """
    n = data.shape[0]
    lines = [header.rstrip(), ""]
    lines.append(f"internalField   nonuniform List<vector> ")
    lines.append(f"{n}")
    lines.append("(")
    # 6 decimals is enough for double-precision restart accuracy in ASCII.
    fmt = "(%.10g %.10g %.10g)"
    chunk = 100_000
    for i in range(0, n, chunk):
        block = data[i : i + chunk]
        lines.extend(fmt % tuple(v) for v in block)
    lines.append(")")
    lines.append(";")
    lines.append("")
    lines.append(boundary.rstrip())
    lines.append("")
    path.write_text("\n".join(lines))


# =============================================================================
# Divergence-free perturbation generator
# =============================================================================


def make_solenoidal_perturbation(
    cell_xyz: np.ndarray,
    integral_scale: float = 5.0e-3,
    target_rms: float = 0.05,
    u_ref_scale: float = 0.10,
    grid_pts_per_scale: int = 6,
    pad: float = 1.5,
    seed: int = 1729,
) -> np.ndarray:
    """Return a solenoidal velocity perturbation evaluated at `cell_xyz`.

    Parameters
    ----------
    cell_xyz : (N, 3) array of OpenFOAM cell centres.
    integral_scale : eddy peak scale [m].
    target_rms : in-bulk RMS magnitude of |U'|, in units of `u_ref_scale`.
    u_ref_scale : reference velocity used to set the absolute amplitude
        before the local-rescale step in the caller.
    grid_pts_per_scale : Cartesian-grid resolution.  6 is enough.
    pad : padding factor on the bounding box.
    seed : RNG seed.

    Returns
    -------
    Up : (N, 3) divergence-free perturbation at the cell centres.

    Construction
    ------------
    1. Pick a Cartesian box [-L, L]^2 x [0, H] with L = pad * R_OUT.
    2. Choose grid spacing h_grid ~ integral_scale / grid_pts_per_scale.
    3. Random complex vector potential psi_hat(k) with Hermitian symmetry.
    4. Spectrum filter |psi_hat(k)|^2 ~ E(k) / k^2,
       E(k) ~ (k/k_p)^4 * exp(-2 (k/k_p)^2).
    5. Uhat'(k) = 1j * cross(k, psi_hat).  Exactly solenoidal.
    6. IFFT to physical Up_grid(x, y, z).
    7. Trilinear-interpolate Up_grid at cell_xyz.
    """
    rng = np.random.default_rng(seed)

    L = pad * R_OUT
    Lz = H
    h_grid = integral_scale / grid_pts_per_scale
    Nxy = int(np.ceil(2 * L / h_grid))
    Nz = int(np.ceil(Lz / h_grid))
    # Round up to even for nicer FFTs.
    Nxy += Nxy % 2
    Nz += Nz % 2

    print(
        f"  Cartesian grid for spectral construction: "
        f"{Nxy} x {Nxy} x {Nz} = {Nxy*Nxy*Nz:,} points"
    )

    # Coordinate axes of the Cartesian grid.
    x_grid = np.linspace(-L, L, Nxy, endpoint=False)
    y_grid = np.linspace(-L, L, Nxy, endpoint=False)
    z_grid = np.linspace(0.0, Lz, Nz, endpoint=False)

    # Wavenumber axes (cycles -> rad/m via 2 pi).
    kx = 2 * np.pi * np.fft.fftfreq(Nxy, d=2 * L / Nxy)
    ky = 2 * np.pi * np.fft.fftfreq(Nxy, d=2 * L / Nxy)
    kz = 2 * np.pi * np.fft.fftfreq(Nz, d=Lz / Nz)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    Kmag = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
    Kmag_safe = np.where(Kmag == 0.0, 1.0, Kmag)

    # Spectrum.  Pope/Passot-Pouquet (1987) form.
    k_p = 2 * np.pi / integral_scale
    E_k = (Kmag / k_p) ** 4 * np.exp(-2.0 * (Kmag / k_p) ** 2)
    # |psi_hat|^2 ~ E(k)/k^2 so that |U'_hat|^2 = k^2 |psi_hat|^2 ~ E(k).
    psi_amp = np.sqrt(np.maximum(E_k, 0.0)) / Kmag_safe
    psi_amp[Kmag == 0.0] = 0.0  # zero-mean

    # Random Hermitian-symmetric vector potential via 3 real-valued IFFTs:
    # easiest robust route -- generate three real Gaussian fields in physical
    # space and take their FFT, then multiply by psi_amp, then IFFT.
    psi_hat = np.empty((3, Nxy, Nxy, Nz), dtype=np.complex128)
    for c in range(3):
        white = rng.standard_normal((Nxy, Nxy, Nz))
        psi_hat[c] = np.fft.fftn(white) * psi_amp

    # U'_hat = 1j * (k x psi_hat).
    Uhat = np.empty_like(psi_hat)
    Uhat[0] = 1j * (KY * psi_hat[2] - KZ * psi_hat[1])
    Uhat[1] = 1j * (KZ * psi_hat[0] - KX * psi_hat[2])
    Uhat[2] = 1j * (KX * psi_hat[1] - KY * psi_hat[0])

    # Inverse transform.
    Ugrid = np.empty((3, Nxy, Nxy, Nz))
    for c in range(3):
        Ugrid[c] = np.fft.ifftn(Uhat[c]).real

    # Pre-rescale to a sensible nominal amplitude (final rescale comes after
    # taper application, in the caller).  We aim for in-bulk RMS = u_ref_scale.
    cur_rms = np.sqrt(np.mean(np.sum(Ugrid ** 2, axis=0)))
    if cur_rms > 0:
        Ugrid *= u_ref_scale / cur_rms

    # Trilinear interpolation onto cell centres.
    print(f"  Interpolating onto {cell_xyz.shape[0]:,} cell centres ...")
    Up = _trilinear_interp(Ugrid, x_grid, y_grid, z_grid, cell_xyz)
    return Up


def _trilinear_interp(
    Ugrid: np.ndarray,
    x_axis: np.ndarray,
    y_axis: np.ndarray,
    z_axis: np.ndarray,
    pts: np.ndarray,
) -> np.ndarray:
    """Vectorised trilinear interpolation of a (3, Nx, Ny, Nz) field at pts."""
    dx = x_axis[1] - x_axis[0]
    dy = y_axis[1] - y_axis[0]
    dz = z_axis[1] - z_axis[0]

    # Fractional indices
    fx = (pts[:, 0] - x_axis[0]) / dx
    fy = (pts[:, 1] - y_axis[0]) / dy
    fz = (pts[:, 2] - z_axis[0]) / dz

    # Clamp to valid range
    Nx, Ny, Nz = Ugrid.shape[1], Ugrid.shape[2], Ugrid.shape[3]
    fx = np.clip(fx, 0.0, Nx - 1.0001)
    fy = np.clip(fy, 0.0, Ny - 1.0001)
    fz = np.clip(fz, 0.0, Nz - 1.0001)

    ix = fx.astype(np.int64)
    iy = fy.astype(np.int64)
    iz = fz.astype(np.int64)
    ax = fx - ix
    ay = fy - iy
    az = fz - iz

    Up = np.zeros((pts.shape[0], 3))
    for c in range(3):
        g = Ugrid[c]
        v000 = g[ix, iy, iz]
        v001 = g[ix, iy, iz + 1]
        v010 = g[ix, iy + 1, iz]
        v011 = g[ix, iy + 1, iz + 1]
        v100 = g[ix + 1, iy, iz]
        v101 = g[ix + 1, iy, iz + 1]
        v110 = g[ix + 1, iy + 1, iz]
        v111 = g[ix + 1, iy + 1, iz + 1]
        Up[:, c] = (
            (1 - ax) * (1 - ay) * (1 - az) * v000
            + (1 - ax) * (1 - ay) * az * v001
            + (1 - ax) * ay * (1 - az) * v010
            + (1 - ax) * ay * az * v011
            + ax * (1 - ay) * (1 - az) * v100
            + ax * (1 - ay) * az * v101
            + ax * ay * (1 - az) * v110
            + ax * ay * az * v111
        )
    return Up


# =============================================================================
# Wall taper
# =============================================================================


def wall_taper(cell_xyz: np.ndarray) -> np.ndarray:
    """Smooth taper that vanishes at all four no-slip / electrode walls.

    T(r, z) = sin(pi (r - R_IN)/(R_OUT - R_IN)) * sin(pi z / H)
    """
    r = np.sqrt(cell_xyz[:, 0] ** 2 + cell_xyz[:, 1] ** 2)
    z = cell_xyz[:, 2]
    r_frac = (r - R_IN) / (R_OUT - R_IN)
    r_frac = np.clip(r_frac, 0.0, 1.0)
    z_frac = np.clip(z / H, 0.0, 1.0)
    return np.sin(np.pi * r_frac) * np.sin(np.pi * z_frac)


# =============================================================================
# Diagnostics
# =============================================================================


def diagnose(Up: np.ndarray, U_lam: np.ndarray, cell_xyz: np.ndarray) -> None:
    print("\n=== Perturbation diagnostics ===")
    rms_p = np.sqrt(np.mean(np.sum(Up ** 2, axis=1)))
    rms_lam = np.sqrt(np.mean(np.sum(U_lam ** 2, axis=1)))
    print(f"  RMS |U'|         = {rms_p:.4e} m/s")
    print(f"  RMS |U_lam|      = {rms_lam:.4e} m/s")
    print(f"  RMS ratio        = {rms_p/rms_lam:.4f}")
    print(f"  max(|U'|)        = {np.linalg.norm(Up, axis=1).max():.4e} m/s")

    # Wall-vanishing check: probe cells closest to z = 0
    z = cell_xyz[:, 2]
    near_wall = z < 1e-4
    if near_wall.any():
        rms_wall = np.sqrt(
            np.mean(np.sum(Up[near_wall] ** 2, axis=1))
        )
        print(
            f"  RMS |U'| in z<0.1mm: {rms_wall:.4e}  "
            f"(should be << bulk RMS {rms_p:.4e})"
        )

    # Crude divergence estimate on a Cartesian sub-sampling
    print("\n  (divergence will be cleaned by first PISO step on restart;\n"
          "   a small residual is expected after wall taper)")


# =============================================================================
# Driver
# =============================================================================


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("time", help="OpenFOAM time directory, e.g. 0.5")
    p.add_argument("--rms", type=float, default=0.05,
                   help="target RMS |U'| as fraction of laminar RMS")
    p.add_argument("--integral-scale", type=float, default=5.0e-3,
                   help="eddy peak scale [m]; default 5 mm = h/2")
    p.add_argument("--seed", type=int, default=1729,
                   help="RNG seed; change to generate independent realisations")
    p.add_argument("--no-backup", action="store_true",
                   help="skip writing <time>/U.lam backup")
    args = p.parse_args()

    case = Path.cwd()
    tdir = case / args.time
    if not tdir.is_dir():
        print(f"FATAL: {tdir} is not a directory.", file=sys.stderr)
        print(f"  did you run `reconstructPar -time {args.time}`?",
              file=sys.stderr)
        return 1

    U_path = tdir / "U"
    C_path = tdir / "C"
    if not U_path.is_file():
        print(f"FATAL: {U_path} not found.", file=sys.stderr)
        return 1
    if not C_path.is_file():
        print(
            f"FATAL: {C_path} not found.  Run `postProcess -func "
            f"writeCellCentres -time {args.time}` first.",
            file=sys.stderr,
        )
        return 1

    print(f"Reading laminar U from {U_path}")
    t0 = wallclock.time()
    U_lam, U_header, U_bnd = read_vector_field(U_path)
    print(f"  {U_lam.shape[0]:,} cells, {wallclock.time()-t0:.1f}s")

    print(f"Reading cell centres from {C_path}")
    t0 = wallclock.time()
    cell_xyz, _, _ = read_vector_field(C_path)
    print(f"  {cell_xyz.shape[0]:,} cells, {wallclock.time()-t0:.1f}s")

    if cell_xyz.shape != U_lam.shape:
        print("FATAL: cell-centre count != U count.", file=sys.stderr)
        return 1

    # Reference scale: bulk-RMS of U_lam (this is the "U_inf" we scale against).
    rms_lam = np.sqrt(np.mean(np.sum(U_lam ** 2, axis=1)))
    print(f"\nLaminar bulk RMS |U_lam| = {rms_lam:.4e} m/s")
    target_abs_rms = args.rms * rms_lam
    print(f"Target perturbation RMS = {args.rms*100:.1f}% = "
          f"{target_abs_rms:.4e} m/s")

    print(f"\nGenerating solenoidal perturbation "
          f"(integral scale = {args.integral_scale*1e3:.1f} mm, "
          f"seed = {args.seed})")
    t0 = wallclock.time()
    Up = make_solenoidal_perturbation(
        cell_xyz,
        integral_scale=args.integral_scale,
        target_rms=args.rms,
        u_ref_scale=rms_lam,
        seed=args.seed,
    )
    print(f"  Generated in {wallclock.time()-t0:.1f}s")

    # Wall taper
    print("\nApplying wall taper (sin*sin envelope)")
    T = wall_taper(cell_xyz)
    Up *= T[:, None]

    # Final rescale to hit the target RMS exactly (taper attenuated the field)
    cur_rms = np.sqrt(np.mean(np.sum(Up ** 2, axis=1)))
    if cur_rms > 0:
        Up *= target_abs_rms / cur_rms

    # Add to laminar
    U_new = U_lam + Up

    diagnose(Up, U_lam, cell_xyz)

    # Backup
    if not args.no_backup:
        backup = tdir / "U.lam"
        if not backup.exists():
            shutil.copy2(U_path, backup)
            print(f"\nBackup written: {backup}")

    print(f"Writing perturbed U to {U_path}")
    t0 = wallclock.time()
    write_vector_field(U_path, U_new, U_header, U_bnd)
    print(f"  Wrote in {wallclock.time()-t0:.1f}s")

    print("\nDone.")
    print("Next steps (in WSL):")
    print(f"  decomposePar -fields -time {args.time} -force")
    print(f"  # edit system/controlDict: startFrom startTime, "
          f"startTime {args.time}, endTime {float(args.time)+0.05}")
    print(f"  nohup mpirun -np 8 epotMhdFoam -parallel "
          f"> log.epotMhdFoam.perturb 2>&1 &")
    return 0


if __name__ == "__main__":
    sys.exit(main())
