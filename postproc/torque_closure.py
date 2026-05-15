#!/usr/bin/env python3
"""
Steady-state torque-closure check for the Corbino disk.

In steady state for the laminar Corbino flow:

    T_EM  +  T_visc  =  0

where

    T_EM   = ∫_V  r · (J × B)_θ  dV     (volume integral, body force)
    T_visc = ∫_∂V r · τ_wθ        dA     (surface integral, on insulatingPlates)

Closure ratio   |T_visc / T_EM|  should be 1 ± a few percent at near-steady
state.  Mismatch sources:

  - residual transient (we're at 99.3% of asymptote, so expect ~1% offset)
  - non-axisymmetric numerical noise (should be << 1% on this 360° mesh)
  - wall-gradient discretisation error (we use one-sided FD into the BL)

Run with the case's foamReader so paths resolve correctly. 

Usage (from the case directory):

    python3 /path/to/torque_closure.py /case/dir 
"""

import sys
import numpy as np
import pyvista as pv
from pathlib import Path


RHO  = 6440.0     # kg/m^3
NU   = 3.73e-7    # m^2/s
MU   = RHO * NU   # dynamic viscosity, Pa·s


# Helpers

def cell_volume_field(mesh: pv.UnstructuredGrid) -> np.ndarray:
    """Per-cell volume from PyVista's compute_cell_sizes filter."""
    sized = mesh.compute_cell_sizes(volume=True, area=False, length=False)
    return np.asarray(sized.cell_data["Volume"])


def cell_centers(mesh: pv.UnstructuredGrid) -> np.ndarray:
    return np.asarray(mesh.cell_centers().points)


def face_areas_and_centers(patch: pv.PolyData):
    """Per-face area and centroid for a polydata patch."""
    sized = patch.compute_cell_sizes(area=True, volume=False, length=False)
    areas = np.asarray(sized.cell_data["Area"])
    centers = np.asarray(patch.cell_centers().points)
    return areas, centers


def az_unit_vector_at(points: np.ndarray) -> np.ndarray:
    """θ̂ unit vector at each (x,y,z); θ̂ = (-sin φ, cos φ, 0)."""
    x, y = points[:, 0], points[:, 1]
    r_xy = np.hypot(x, y)
    r_xy[r_xy == 0] = 1.0  # axis singularity guard (no points should be there)
    return np.column_stack((-y / r_xy, x / r_xy, np.zeros_like(x)))


# ----------------------------------------------------------------------- #
# Main
# ----------------------------------------------------------------------- #
def main(case_dir: Path):
    foam_file = case_dir / "case.foam"
    if not foam_file.exists():
        foam_file.touch()

    print(f"[load] {foam_file}")
    reader = pv.OpenFOAMReader(str(foam_file))
    times = reader.time_values
    print(f"[load] {len(times)} time steps available; using last = {times[-1]}")
    reader.set_active_time_value(times[-1])

    # Make sure we get the patches as well as the internal mesh
    reader.disable_all_patch_arrays()
    reader.enable_patch_array("internalMesh")
    reader.enable_patch_array("patch/innerElectrode")
    reader.enable_patch_array("patch/outerElectrode")
    reader.enable_patch_array("patch/insulatingPlates")

    multi = reader.read()
    internal = multi["internalMesh"]
    print(f"[load] internalMesh: {internal.n_cells} cells, "
          f"fields = {list(internal.cell_data.keys())}")

   
    def find_patch(mb, name):
        try:
            keys = list(mb.keys())
        except Exception:
            return None
        for k in keys:
            if k.endswith(name) or k == name:
                return mb[k]
        for k in keys:
            sub = mb[k]
            if hasattr(sub, "keys"):
                hit = find_patch(sub, name)
                if hit is not None:
                    return hit
        return None

    plates = find_patch(multi, "insulatingPlates")
    if plates is None:
        print("[fatal] could not locate insulatingPlates patch.")
        print("[fatal] Top-level MultiBlock keys:")
        for k in multi.keys():
            print(f"        {k!r}")
            sub = multi[k]
            if hasattr(sub, "keys"):
                for kk in sub.keys():
                    print(f"          - {kk!r}")
        sys.exit(1)


    # Volume integral T_EM = ∫ r · (J × B)_θ dV

    J = np.asarray(internal.cell_data["J"])
    B = np.asarray(internal.cell_data["B"])
    JxB = np.cross(J, B)                              # body force / σ-units
    centers = cell_centers(internal)
    vols = cell_volume_field(internal)

    that = az_unit_vector_at(centers)
    JxB_theta = np.einsum("ij,ij->i", JxB, that)      # azimuthal component
    r_xy = np.hypot(centers[:, 0], centers[:, 1])
    T_EM = float(np.sum(r_xy * JxB_theta * vols))

    print()
    print(f"[EM] (J × B)_θ stats:  min={JxB_theta.min():.3e}  "
          f"max={JxB_theta.max():.3e}  mean={JxB_theta.mean():.3e}")
    print(f"[EM] T_EM = ∫ r·(JxB)_θ dV = {T_EM:+.4e} N·m")

    # ------------------------------------------------------------------- #
    # Wall shear torque on insulatingPlates
    #
    # For each wall face:
    #   - find the adjacent cell centre (nearest neighbour via KDTree)
    #   - U_θ_cell = U_cell · θ̂
    #   - dn = |z_cell − z_wall|
    #   - τ_wθ_on_fluid = -μ · U_θ_cell / dn
    #
    # Sign derivation: traction t_θ = σ_zθ · n_z with σ_zθ = μ ∂u_θ/∂z.
    #   - top wall (z=h, n_z = +1):  ∂u_θ/∂z|_w ≈ (0 − U_θ_cell)/(h − z_cell)
    #     → t_θ = +μ · (−U_θ_cell)/dn = −μ·U_θ_cell/dn
    #   - bottom wall (z=0, n_z = −1):  ∂u_θ/∂z|_w ≈ (U_θ_cell − 0)/z_cell
    #     → t_θ = −μ · (+U_θ_cell)/dn = −μ·U_θ_cell/dn
    # Both walls collapse to the same closed-form expression.
    #
    # First-cell height is 6.6 μm and δ_H ≈ 53 μm, so the FD samples
    # well inside the linear sub-region of the Hartmann BL.
    # ------------------------------------------------------------------- #
    print()
    print("[visc] one-sided FD into the Hartmann BL ...")
    print(f"[visc] insulatingPlates: {plates.n_cells} faces")

    areas, fcenters = face_areas_and_centers(plates)
    that_face = az_unit_vector_at(fcenters)
    r_face = np.hypot(fcenters[:, 0], fcenters[:, 1])

    # KDTree over internal cell centres → nearest cell for each wall face
    from scipy.spatial import cKDTree
    U_cells = np.asarray(internal.cell_data["U"])
    tree = cKDTree(centers)
    _, idxs = tree.query(fcenters, k=1)
    U_near = U_cells[idxs]
    z_near = centers[idxs, 2]
    z_face = fcenters[:, 2]
    dn = np.abs(z_near - z_face)

    print(f"[visc] dn stats:  min={dn.min()*1e6:.2f} μm  "
          f"max={dn.max()*1e6:.2f} μm  mean={dn.mean()*1e6:.2f} μm")

    # azimuthal component of U at the near-wall cell centre
    U_theta_near = np.einsum("ij,ij->i", U_near, that_face)

    # τ_wθ on fluid (positive when fluid swirls −θ̂ — wall drags it +θ̂)
    tau_wall_theta = -MU * U_theta_near / dn

    T_visc = float(np.sum(r_face * tau_wall_theta * areas))

    print(f"[visc] τ_wθ stats:  min={tau_wall_theta.min():.3e}  "
          f"max={tau_wall_theta.max():.3e}  mean={tau_wall_theta.mean():.3e}")
    print(f"[visc] T_visc = ∫ r·τ_wθ dA = {T_visc:+.4e} N·m")


    #  Closure

    total = T_EM + T_visc
    rel_err = abs(total) / max(abs(T_EM), abs(T_visc))

    print()
    print("=" * 60)
    print(f"  T_EM         = {T_EM:+.4e} N·m   (volumetric driving)")
    print(f"  T_visc       = {T_visc:+.4e} N·m   (wall drag)")
    print(f"  T_EM + T_visc= {total:+.4e} N·m")
    print(f"  relative err = {100 * rel_err:.2f} %")
    print("=" * 60)

    # Diagnostic
 
    R_MID = 0.5 * (0.02 + 0.10)   # halfway between r_in and r_out
    inner_side = r_xy < R_MID
    T_EM_inner = float(np.sum(r_xy[inner_side] * JxB_theta[inner_side]
                              * vols[inner_side]))
    T_EM_outer = float(np.sum(r_xy[~inner_side] * JxB_theta[~inner_side]
                              * vols[~inner_side]))
    print(f"[diag] T_EM r<{R_MID*100:.0f} cm  = {T_EM_inner:+.4e}")
    print(f"[diag] T_EM r>{R_MID*100:.0f} cm  = {T_EM_outer:+.4e}")
    print(f"[diag] (split is volumetric, not by patch)")


if __name__ == "__main__":
    case = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()
    main(case.resolve())
