"""
Charge carriers other than gas ions: electrosprayed droplets ("wet emitter" stage).

corona.py assumes the carriers are molecular ions: they drift at mu*E ~ 100 m/s, the
air moves at ~2 m/s, so the air is ignored and one unit-voltage solution gives the
whole I-V curve by similarity.  A charged droplet is 10x slower, which is the point
(thrust per watt is 1 / carrier speed) but breaks both shortcuts:

    div(grad(phi))               = -rho / eps0
    div(rho * (Z * E + u))       = 0            Z = droplet mobility, u = air velocity

This module adds, on top of CoronaProblem and on the same mesh,

  * u: potential flow through the strip row (same stiffness matrix, impermeable
    electrodes), scaled so that the superficial duct velocity is u0;
  * an outflow boundary: charge that gets past the middle of the collector is taken
    out of the cell and reported as `leak` instead of being wrapped round into the next
    (periodic) stage.  The periodic cell of corona.py accelerates such charge a second
    time, which is right for a stack whose potential keeps falling but not for the
    alternating-polarity stack of this design, where it meets a reversed field;
  * the emitter law for a wetted fibre: emission starts when the electric stress on
    the liquid exceeds its capillary pressure, E_on = sqrt(2*gamma / (eps0*r)), and
    the surface field is taken to stay there while it sprays (the same closure as
    Kaptzov's for a corona, with a different threshold).

In unit-voltage variables the carrier velocity is Z*V*(E' + w*u_hat) with
w = u0 / (Z*V), so every (V, u0) pair needs its own solve; stage_curve_carrier() finds,
for each voltage, the current at which the emitter field sits at its threshold.

Droplet mobility follows from Stokes-Cunningham drag and the Rayleigh charge limit:
droplet_mobility(), rayleigh_charge(), mobility_band().
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from corona import EPS0, CoronaProblem, Mesh, StageCurve, _finish_gmsh

E_CHARGE = 1.602176634e-19
ETA_AIR = 1.81e-5       # Pa s
MFP_AIR = 67.3e-9       # m, mean free path at 20 C, 1 atm


# ------------------------------------------------------------------ droplet physics
def cunningham(d_p: float) -> float:
    """Slip correction for a sphere of diameter d_p in air (Allen & Raabe coefficients)."""
    kn = 2.0 * MFP_AIR / d_p
    return 1.0 + kn * (1.142 + 0.558 * np.exp(-0.999 / kn))


def droplet_mobility(d_p: float, q: float) -> float:
    """Electrical mobility [m^2/(V s)] of a droplet of diameter d_p carrying charge q."""
    return q * cunningham(d_p) / (3.0 * np.pi * ETA_AIR * d_p)


def rayleigh_charge(d_p: float, gamma: float) -> float:
    """Largest charge a droplet of diameter d_p and surface tension gamma can hold."""
    return np.pi * np.sqrt(8.0 * EPS0 * gamma * d_p**3)


def spray_onset_field(r: float, gamma: float) -> float:
    """Field at which a liquid film on a fibre of radius r starts to spray.

    Electric stress eps0*E^2/2 against the capillary pressure gamma/r of the wetted
    cylinder.  For gamma ~ 0.04 N/m this equals Peek's corona field for a thin wire,
    which is the classical reason why only liquids below ~0.05 N/m can be
    electrosprayed in air without a discharge.
    """
    return np.sqrt(2.0 * gamma / (EPS0 * r))


def mobility_band(gamma: float, d_range=(50e-9, 300e-9), frac=(0.5, 1.0)) -> tuple[float, float]:
    """Lowest and highest mobility of droplets in a size range charged to a fraction of the Rayleigh limit."""
    z = [droplet_mobility(d, f * rayleigh_charge(d, gamma)) for d in np.geomspace(*d_range, 12) for f in frac]
    return float(min(z)), float(max(z))


# ---------------------------------------------------------------------------- mesh
def mesh_slab(d: float, h: float, n: int = 60) -> Mesh:
    """Rectangle 0<x<d, 0<y<h: emitting plane at x=0, collecting plane at x=d (1-D test case)."""
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        occ = gmsh.model.occ
        box = occ.addRectangle(0, 0, 0, d, h)
        occ.synchronize()
        groups = {"emitter": [], "collector": []}
        for dim, tag in gmsh.model.getBoundary([(2, box)], oriented=False):
            lo, hi = gmsh.model.getParametrizationBounds(1, tag)
            q = np.asarray(gmsh.model.getValue(1, tag, list(np.linspace(lo[0], hi[0], 5)))).reshape(-1, 3)
            if np.allclose(q[:, 0], 0.0, atol=1e-12):
                groups["emitter"].append(tag)
            elif np.allclose(q[:, 0], d, atol=1e-12):
                groups["collector"].append(tag)
        gmsh.option.setNumber("Mesh.MeshSizeMin", d / n)
        gmsh.option.setNumber("Mesh.MeshSizeMax", d / n)
        gmsh.option.setNumber("Mesh.Algorithm", 5)
        gmsh.model.mesh.generate(2)
        return _finish_gmsh(gmsh, groups, None, 1, {"kind": "slab", "d": d, "pitch": h, "h": h})
    finally:
        gmsh.finalize()


# -------------------------------------------------------------------------- solver
class CarrierProblem(CoronaProblem):
    """CoronaProblem with carriers that are also convected by the air.

    self.w = u0 / (Z * V) [1/m] sets the strength of convection relative to drift in
    the unit-voltage problem; w = 0 and outflow=False reproduce CoronaProblem exactly.
    Currents returned by diagnostics() are normalised as in corona.py:  I = Z * V^2 * i.
    """

    def __init__(self, mesh: Mesh, outflow: bool = True):
        super().__init__(mesh)
        self.w = 0.0
        periodic = bool(np.any(mesh.joff != 0))
        self.eU = self._unit_flow_flux(periodic)
        # wrap[e] = +1 / -1 on edges that lead from the periodic boundary into the next cell
        # (sign: direction ei -> ej or ej -> ei), 0 elsewhere.  Carriers leaving through them
        # are discarded (outflow) and counted as leak.
        self.wrap = self._wrap_edges() if (outflow and periodic) else np.zeros(len(self.ei))

    def _wrap_edges(self) -> np.ndarray:
        m = self.mesh
        t, rid, n = m.t, m.rid, self.n
        on_bnd = np.zeros(n, bool)
        on_bnd[rid[m.joff != 0]] = True                      # reduced ids of the periodic boundary
        left_orig = on_bnd[rid] & (m.joff == 0) & (np.abs(m.p[:, 0] - m.info["x0"]) < 1e-9)
        ea = t[:, [0, 1, 2]].T.ravel()
        eb = t[:, [1, 2, 0]].T.ravel()
        hit = left_orig[ea] ^ left_orig[eb]                  # one end on the left boundary, the other inside
        hit &= ~(on_bnd[rid[ea]] & on_bnd[rid[eb]])
        ra, rb = rid[ea[hit]], rid[eb[hit]]
        bnd_first = left_orig[ea[hit]]                       # True: boundary node is `a`
        lo, hi = np.minimum(ra, rb), np.maximum(ra, rb)
        b = np.where(bnd_first, ra, rb)                      # reduced id of the boundary end
        sign = np.where(b == lo, 1.0, -1.0)
        idx = np.searchsorted(self.ei.astype(np.int64) * n + self.ej, lo.astype(np.int64) * n + hi)
        wrap = np.zeros(len(self.ei))
        wrap[idx] = sign
        return wrap

    # --- air flow
    def _unit_flow_flux(self, periodic: bool) -> np.ndarray:
        """Volume flux of air along every mesh edge [m^2/s per m depth] for 1 m/s superficial velocity.

        Periodic cells: potential flow, velocity potential quasi-periodic, all solid
        surfaces impermeable (natural boundary condition), so the flux is discretely
        divergence free at every node.  Non-periodic meshes (the slab test): uniform flow
        along +x.
        """
        m = self.mesh
        if not periodic:
            pot = -m.p[:, 0]                       # u = -grad(pot) = (1, 0)
            pr = np.zeros(self.n)
            pr[m.rid] = pot
            return -(self.eK * (pr[self.ei] - pr[self.ej]))
        # K pot = -gjump with one dof pinned (pure Neumann problem)
        free = np.arange(1, self.n)
        pot = np.zeros(self.n)
        pot[free] = spla.spsolve(self.K[free][:, free].tocsc(), -self.gjump[free])
        flux = -(self.eK * (pot[self.ei] - pot[self.ej]) + self.eKj)
        # superficial velocity = (1 / (S * H)) * integral(u_x dA)
        ph = pot[m.rid] + m.joff
        p, t = m.p, m.t
        x, y = p[t, 0], p[t, 1]
        b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], axis=1)
        ux_area = -0.5 * (b * ph[t]).sum(1)        # u_x * area of each triangle
        S = m.info["stage_pitch"]
        H = 0.5 * m.info["pitch"]
        return flux / (ux_area.sum() / (S * H))

    def total_flux(self, phi):
        return self.edge_flux(phi) + self.w * self.eU

    # --- overrides: same structure as CoronaProblem, with drift + convection
    def emitter_state(self, phi, rho):
        tot = self.total_flux(phi)
        out = (np.bincount(self.ei, weights=np.maximum(tot, 0.0), minlength=self.n)
               + np.bincount(self.ej, weights=np.maximum(-tot, 0.0), minlength=self.n))[self.em].sum()
        flux = self.edge_flux(phi)
        e = (self.node_outflux(flux)[self.em] - (self.M * rho / EPS0)[self.em]).sum() / self.len_em
        return out, e

    def solve_rho(self, phi, rho_old, rho_e):
        flux = self.edge_flux(phi)
        tot = flux + self.w * self.eU
        fp = np.maximum(tot, 0.0)
        fm = np.maximum(-tot, 0.0)
        diag = np.bincount(self.ei, weights=fp, minlength=self.n) + np.bincount(self.ej, weights=fm, minlength=self.n)
        fb = self.M * rho_old / EPS0 - self.node_outflux(flux)   # electric flux into electrode surfaces
        diag[self.absorb] += np.maximum(fb[self.absorb], 0.0)
        diag += 1e-12 * (np.abs(tot).max() + 1e-300)
        rows = np.concatenate([self.ej, self.ei])
        cols = np.concatenate([self.ei, self.ej])
        vals = np.concatenate([-np.where(self.wrap > 0, 0.0, fp), -np.where(self.wrap < 0, 0.0, fm)])
        keep = ~self.is_em[rows]
        d = diag.copy()
        d[self.em] = 1.0
        A = sp.coo_matrix((vals[keep], (rows[keep], cols[keep])), shape=(self.n, self.n)).tocsr() + sp.diags(d)
        rhs = np.zeros(self.n)
        rhs[self.em] = rho_e
        return np.maximum(spla.spsolve(A.tocsc(), rhs), 0.0)

    def diagnostics(self, phi, rho):
        s = self.mesh.symmetry
        flux = self.edge_flux(phi)
        tot = flux + self.w * self.eU
        up = np.where(tot > 0, rho[self.ei], rho[self.ej])        # upwind density on each edge
        lost = (self.wrap * tot) > 0                               # leaves through the outflow boundary
        i_leak = np.abs(tot * up)[lost].sum()
        
        i_em = self.node_outflux(tot * up)[self.em].sum()
        fb = self.M * rho / EPS0 - self.node_outflux(flux)
        i_co = (np.maximum(fb[self.co], 0.0) * rho[self.co]).sum()
        i_gd = (np.maximum(fb[self.gd], 0.0) * rho[self.gd]).sum() if len(self.gd) else 0.0
        fx = (np.where(lost, 0.0, flux * up) * self.edx).sum()     # integral(rho E_x dA) inside the cell
        _, e_em = self.emitter_state(phi, rho)
        with np.errstate(divide="ignore", invalid="ignore"):
            e_co = np.where(self.len_co_nodes[self.co] > 0, fb[self.co] / self.len_co_nodes[self.co], 0.0)
        # i_gd also carries the leak so that CoronaProblem.solve_unit's charge balance still closes
        return dict(i_em=s * i_em, i_co=s * i_co, i_gd=s * (i_gd + i_leak), i_back=s * i_gd, i_leak=s * i_leak,
                    fx=s * fx, e_em=float(e_em), e_co_max=float(e_co.max()), e_co=e_co)


def stage_curve_carrier(prob: CarrierProblem, E_on: float, Z: float, u0: float, V=None, n_pts: int = 14,
                        v_max_over_onset: float = 3.0, verbose: bool = False) -> StageCurve:
    """I-V-thrust curve of one stage whose carriers have mobility Z in air moving at u0.

    For each voltage the unit-voltage problem is solved at w = u0/(Z V) for the current
    that leaves the emitter surface field at E_on.  Returns a StageCurve whose meta
    holds the leak fraction (charge blown through the collector row) at every point.
    """
    prob.w = 0.0
    _, _, d0, _ = prob.solve_unit(0.0)
    f0 = d0["e_em"]
    v_on = E_on / f0
    if V is None:
        V = v_on * np.geomspace(1.0, v_max_over_onset, n_pts + 1)[1:]
    V = np.sort(np.asarray(V, float))
    info = prob.mesh.info
    i_ref = 9.0 / 8.0 * EPS0 * info["pitch"] / info["d"] ** 3
    Vs, Is, Fs, EC, LK = [v_on], [0.0], [0.0], [d0["e_co_max"] * v_on], [0.0]
    rho, i_prev = None, 0.0
    for v in V:
        if v <= v_on * 1.0005:
            continue
        prob.w = u0 / (Z * v)
        target = E_on / (v * f0)                 # wanted f / f0  (< 1)
        lo, f_lo = 0.0, 1.0                      # f(i) decreases monotonically with i
        hi, f_hi = None, None
        i = max(1.5 * i_prev, 0.004 * i_ref)
        best = None
        for _ in range(40):
            phi, rho_new, dg, it = prob.solve_unit(i, rho0=rho)
            f = dg["e_em"] / f0
            if it < 0 or f <= 0.0:               # beyond the space-charge limit for this w
                hi, f_hi = i, 0.0
            else:
                rho = rho_new
                if abs(f - target) < 2e-4:
                    best = (i, dg)
                    break
                if f > target:
                    lo, f_lo = i, f
                else:
                    hi, f_hi = i, f
            if hi is None:
                slope = (f_lo - 1.0) / lo if lo > 0 else None
                i = lo + (target - f_lo) / slope if slope and slope < 0 else 2.0 * i
                i = min(i, 3.0 * max(lo, 0.004 * i_ref))
            else:
                t = (f_lo - target) / (f_lo - f_hi) if f_lo > f_hi else 0.5
                i = lo + (hi - lo) * float(np.clip(t, 0.15, 0.85))
            if hi is not None and (hi - lo) < 1e-5 * hi:
                best = (i, dg) if it >= 0 else None
                break
        if best is None:
            if verbose:
                print(f"  V={v/1e3:6.2f} kV: no solution (space-charge limit)")
            break
        i, dg = best
        i_prev = i
        Vs.append(v)
        Is.append(Z * dg["i_em"] * v * v)
        Fs.append(dg["fx"] * v * v)
        EC.append(dg["e_co_max"] * v)
        LK.append(dg["i_leak"] / max(dg["i_em"], 1e-300))
        if verbose:
            print(f"  V={v/1e3:6.2f} kV  w={prob.w:7.2f}/m  I={Is[-1]*1e6:8.3f} uA/m  F={Fs[-1]*1e3:7.3f} mN/m  "
                  f"F/P={Fs[-1]/max(Is[-1]*v,1e-300)*1e3:7.1f} N/kW  leak={LK[-1]:.3f}  back={dg['i_back']/max(dg['i_em'],1e-300):.3f}  "
                  f"Ecoll={EC[-1]/1e6:5.2f} MV/m")
    return StageCurve(V=np.array(Vs), I=np.array(Is), F=np.array(Fs), E_coll=np.array(EC), V_onset=v_on,
                      meta=dict(info, E_c=E_on, mu=Z, u0=u0, nodes=prob.n, symmetry=prob.mesh.symmetry, leak=LK))


# ---------------------------------------------------------- analytic reference (1-D)
def slab_analytic(d: float, E0: float, Z: float, u: float, J: np.ndarray):
    """Planar unipolar flow with drift and convection: (V, pressure rise) for current density J.

    (Z E + u) dE/dx = J / eps0 with E(0) = E0  ->  Z E^2/2 + u E = J x / eps0 + c.
    """
    J = np.atleast_1d(np.asarray(J, float))
    x = np.linspace(0.0, d, 4001)
    c = 0.5 * Z * E0**2 + u * E0
    E = (-u + np.sqrt(u * u + 2.0 * Z * (J[:, None] * x[None, :] / EPS0 + c))) / Z
    V = np.trapezoid(E, x, axis=1)
    dp = 0.5 * EPS0 * (E[:, -1] ** 2 - E0**2)
    return V, dp
