"""
2-D drift-region model of a DC corona discharge (the "plasma" in an
electroaerodynamic thruster).

Physics
-------
Outside the thin ionisation sheath that hugs the emitter wire, the gap is filled
with a unipolar cloud of ions drifting along the electric field:

    div(grad(phi))    = -rho / eps0              (Poisson)
    div(rho * mu * E) = 0,   E = -grad(phi)      (steady ion continuity, drift only)

Ion drift (~100 m/s) is far faster than the air (~3 m/s), so convection and
diffusion of ions are neglected.  The sheath is represented by Kaptzov's
hypothesis: once the corona is lit, the field on the emitter surface stays
pinned at Peek's onset field E_c, which fixes how much charge is injected.

The Coulomb body force on the air is f = rho * E.  Its integral is the thrust
reacted on the electrodes (before aerodynamic losses, handled in duct.py).

Numerics
--------
* P1 finite elements on an unstructured gmsh triangulation for phi.
* Node-centred, first-order upwind finite volumes for rho, using the *same*
  discrete fluxes as the FEM stiffness matrix (Phi_ij = -K_ij (phi_i - phi_j)),
  so the discrete Gauss law and charge conservation hold to round-off.
* Similarity: (phi, rho) -> (a phi, a rho) maps solutions to solutions, so the
  problem is solved at unit voltage for a sweep of ion current and every solve
  yields one point of the I-V curve:  V = E_c / f(i).

All quantities are per unit electrode length (2-D planar).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

EPS0 = 8.8541878128e-12
MU_ION = 2.0e-4  # m^2/(V s), positive ions in ambient air (literature 1.8-2.2e-4)


def peek_field(r: float, delta: float = 1.0, m: float = 1.0) -> float:
    """Peek's corona-onset surface field for a round wire of radius r [m] in air.

    E_c = 3.1 MV/m * delta * m * (1 + 0.0308 / sqrt(delta * r))
    delta = relative air density, m = surface roughness factor (1 = polished).
    """
    return 3.1e6 * delta * m * (1.0 + 0.0308 / np.sqrt(delta * r))


# --------------------------------------------------------------------------- mesh
@dataclass
class Mesh:
    """Triangulation plus boundary tags.  `rid` maps nodes to reduced (periodic) dofs."""

    p: np.ndarray  # (N, 2) node coordinates
    t: np.ndarray  # (M, 3) triangles
    emitter: np.ndarray  # nodes on the emitter wire (phi = V, injects ions)
    guard: np.ndarray  # nodes on non-emitting conductors at emitter potential
    collector: np.ndarray  # nodes on conductors at collector potential (phi = 0)
    rid: np.ndarray  # (N,) reduced dof id for each node
    n_red: int
    joff: np.ndarray  # (N,) potential offset of each node relative to its dof (quasi-periodicity)
    symmetry: int = 1  # how many copies of the mesh make one physical unit
    info: dict = field(default_factory=dict)


def _curve_samples(gmsh, tag: int, n: int = 7) -> np.ndarray:
    """Points along a model curve (robust alternative to padded OCC bounding boxes)."""
    lo, hi = gmsh.model.getParametrizationBounds(1, tag)
    u = np.linspace(lo[0], hi[0], n)
    return np.asarray(gmsh.model.getValue(1, tag, list(u))).reshape(-1, 3)[:, :2]


def _size_fields(gmsh, specs, h_max: float, grow: float = 0.25):
    """Graded mesh size: h = h0 on each listed curve set, growing linearly with distance."""
    f = gmsh.model.mesh.field
    ids = []
    for curves, h0 in specs:
        if not curves:
            continue
        d = f.add("Distance")
        f.setNumbers(d, "CurvesList", curves)
        f.setNumber(d, "Sampling", 200)
        m = f.add("MathEval")
        f.setString(m, "F", f"Min({h0} + {grow}*F{d}, {h_max})")
        ids.append(m)
    mn = f.add("Min")
    f.setNumbers(mn, "FieldsList", ids)
    f.setAsBackgroundMesh(mn)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.Algorithm", 5)  # Delaunay -> K_ij <= 0 almost everywhere


def _finish_gmsh(gmsh, groups: dict[str, list[int]], pair, symmetry: int, info: dict) -> Mesh:
    """Extract nodes/triangles/boundary node sets from the current gmsh model.

    pair = None, or (key_left(p), key_right(p), jump): boundary nodes on the 'left' and
    'right' curves are identified after sorting by the given keys, and the right copy
    carries the potential offset `jump` relative to the left one.
    """
    tags, coords, _ = gmsh.model.mesh.getNodes()
    order = np.argsort(tags)
    tags = np.asarray(tags, dtype=np.int64)[order]
    p = np.asarray(coords).reshape(-1, 3)[order][:, :2].copy()
    lookup = np.full(int(tags.max()) + 1, -1, dtype=np.int64)
    lookup[tags] = np.arange(len(tags))
    _, _, enodes = gmsh.model.mesh.getElements(2)
    t = lookup[np.asarray(enodes[0], dtype=np.int64)].reshape(-1, 3)

    def nodes_of(curves):
        out = []
        for c in curves:
            nt, _, _ = gmsh.model.mesh.getNodes(1, c, includeBoundary=True)
            out.append(lookup[np.asarray(nt, dtype=np.int64)])
        return np.unique(np.concatenate(out)) if out else np.array([], dtype=np.int64)

    n = len(p)
    rep = np.arange(n)
    joff = np.zeros(n)
    if pair is not None:
        key_l, key_r, jump = pair
        left = nodes_of(groups["left"])
        right = nodes_of(groups["right"])
        left = left[np.argsort(key_l(p[left]))]
        right = right[np.argsort(key_r(p[right]))]
        if len(left) != len(right) or not np.allclose(key_l(p[left]), key_r(p[right]), atol=1e-8):
            raise RuntimeError("periodic boundaries do not match")
        rep[right] = left
        joff[right] = jump
    uniq, rid = np.unique(rep, return_inverse=True)
    a = p[t[:, 1]] - p[t[:, 0]]
    b = p[t[:, 2]] - p[t[:, 0]]
    flip = (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0]) < 0
    t[flip] = t[flip][:, [0, 2, 1]]
    return Mesh(p=p, t=t, emitter=nodes_of(groups["emitter"]), guard=nodes_of(groups.get("guard", [])),
                collector=nodes_of(groups["collector"] + groups.get("far", [])), rid=rid,
                n_red=len(uniq), joff=joff, symmetry=symmetry, info=info)


def mesh_coaxial(r0: float, R: float, n_wire: int = 16, h_max_frac: float = 0.03) -> Mesh:
    """Quarter cross-section of a wire (radius r0) inside a coaxial cylinder (radius R)."""
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        occ = gmsh.model.occ
        big = occ.addDisk(0, 0, 0, R, R)
        small = occ.addDisk(0, 0, 0, r0, r0)
        box = occ.addRectangle(0, 0, 0, 2 * R, 2 * R)
        ann, _ = occ.cut([(2, big)], [(2, small)])
        q, _ = occ.intersect(ann, [(2, box)])
        occ.synchronize()
        groups = {"emitter": [], "collector": []}
        for dim, tag in gmsh.model.getBoundary(q, oriented=False):
            rr = np.linalg.norm(_curve_samples(gmsh, tag), axis=1)
            if np.allclose(rr, r0, rtol=1e-3):
                groups["emitter"].append(tag)
            elif np.allclose(rr, R, rtol=1e-3):
                groups["collector"].append(tag)
        h_w = 0.5 * np.pi * r0 / n_wire  # quarter circle -> n_wire segments
        _size_fields(gmsh, [(groups["emitter"], h_w)], h_max_frac * R, grow=0.22)
        gmsh.model.mesh.generate(2)
        return _finish_gmsh(gmsh, groups, None, 4, {"r0": r0, "R": R, "kind": "coaxial"})
    finally:
        gmsh.finalize()


def _nose_blade(occ, x_nose: float, x_end: float, yc: float, t_c: float):
    """Round-nosed strip from x_nose to x_end (square cut at x_end)."""
    rc = 0.5 * t_c
    parts = [(2, occ.addDisk(x_nose + rc, yc, 0, rc, rc))]
    if x_end > x_nose + rc + 1e-9:
        parts.append((2, occ.addRectangle(x_nose + rc, yc - rc, 0, x_end - x_nose - rc, t_c)))
    return parts


def mesh_cell(
    d: float,
    pitch: float,
    r_e: float,
    t_c: float,
    c_c: float,
    c_d: float = 0.0,
    mode: str = "cascade",
    upstream: float | None = None,
    downstream: float | None = None,
    n_wire: int = 14,
    h_max: float | None = None,
    h_coll: float | None = None,
    split: float = 0.5,
) -> Mesh:
    """Half unit cell of one thruster stage.

    x is the flow/thrust axis, y the transverse axis.  The emitter wire (radius r_e) is
    at the origin.  The collector is a round-nosed strip of thickness t_c whose nose is
    at x = d and whose conductive face extends c_c downstream.  Rows repeat in y every
    `pitch`; by symmetry only 0 <= y <= pitch/2 is meshed.

    mode="isolated": one bare wire row in front of one collector row (c_d unused).  A
        grounded plane `upstream` ahead of the wires stands in for distant surroundings.

    mode="cascade": an infinite stack of identical "integrated strips".  Each strip is
        [conductive collector, length c_c] + [bare insulating tail, length c_d] + [the
        emitter wire of the next stage bonded along the trailing edge], with the wire
        wired to the strip's own collector.  The potential drops by V per stage, so the
        cell (cut through the conductive parts) is quasi-periodic: phi(x+S) = phi(x) - V
        with S = c_c + c_d + d.  The wire's rear half is buried in the insulating tail
        and its collector upstream is at the same potential, so no ions can flow back
        towards the previous stage.  The tail is treated as charged-up insulation
        (zero normal field, no ion flux); its permittivity (foam, eps_r ~ 1.03) is
        neglected.  `split` places the periodic cut at that fraction of the collector's
        conductive chord (0.5 = mid-chord).
    """
    import gmsh

    H = 0.5 * pitch
    rc = 0.5 * t_c
    h_max = h_max or pitch / 40.0
    h_coll = h_coll or min(t_c / 5.0, h_max)
    c = c_c + c_d
    if mode == "cascade":
        if c_d < 2 * rc or c_c < 3 * rc:
            raise ValueError("need c_d >= t_c and c_c >= 1.5 t_c")
        S = c + d
        x0 = -c + split * c_c
        x1 = x0 + S
    elif mode == "isolated":
        x0 = -(upstream if upstream is not None else 4.0 * d)
        x1 = d + c_c + (downstream if downstream is not None else 2.0 * d)
        S = None
    else:
        raise ValueError(mode)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        occ = gmsh.model.occ
        box = occ.addRectangle(x0, 0, 0, x1 - x0, H)
        tools = [(2, occ.addDisk(0, 0, 0, r_e, r_e))]
        if mode == "cascade":
            pad = 0.25 * c_c
            tools.append((2, occ.addRectangle(d + rc, -rc, 0, x1 + pad - d - rc, t_c)))  # next collector
            tools.append((2, occ.addDisk(d + rc, 0, 0, rc, rc)))
            tools.append((2, occ.addRectangle(x0 - pad, -rc, 0, -c_d - (x0 - pad), t_c)))  # own collector
            tools.append((2, occ.addRectangle(-c_d, -rc, 0, c_d - rc, t_c)))  # insulating tail
            tools.append((2, occ.addDisk(-rc, 0, 0, rc, rc)))
        else:
            tools += _nose_blade(occ, d, d + c_c - rc, 0.0, t_c)
            if c_c > t_c * 1.001:
                tools.append((2, occ.addDisk(d + c_c - rc, 0, 0, rc, rc)))
        dom, _ = occ.cut([(2, box)], tools)
        occ.synchronize()
        groups = {"emitter": [], "guard": [], "collector": [], "far": [], "left": [], "right": [],
                  "insul": []}
        x_ce = (x1 if mode == "cascade" else d + c_c - rc)
        for dim, tag in gmsh.model.getBoundary(dom, oriented=False):
            q = _curve_samples(gmsh, tag)
            qi = q[1:-1]  # interior samples (end points may sit on a junction)
            dist_c = np.hypot(qi[:, 0] - np.clip(qi[:, 0], d + rc, x_ce), qi[:, 1])
            if np.allclose(q[:, 0], x0, atol=1e-9):
                groups["left"].append(tag)
            elif np.allclose(q[:, 0], x1, atol=1e-9):
                groups["right"].append(tag)
            elif np.allclose(np.hypot(q[:, 0], q[:, 1]), r_e, rtol=1e-3):
                groups["emitter"].append(tag)
            elif np.allclose(dist_c, rc, rtol=1e-3) and np.all(qi[:, 0] > 0.5 * d):
                groups["collector"].append(tag)
            elif mode == "cascade" and np.allclose(qi[:, 1], rc, atol=1e-9) and np.all(qi[:, 0] < -c_d):
                groups["guard"].append(tag)
            elif mode == "cascade" and np.all(qi[:, 0] > -c_d - 1e-9) and np.all(qi[:, 0] < r_e) \
                    and np.all(qi[:, 1] > 1e-9):
                groups["insul"].append(tag)
        if not groups["emitter"] or not groups["collector"]:
            raise RuntimeError(f"boundary tagging failed: {groups}")
        pair = None
        if mode == "cascade":
            if not groups["guard"] or len(groups["left"]) != 1 or len(groups["right"]) != 1:
                raise RuntimeError(f"boundary tagging failed: {groups}")
            aff = [1, 0, 0, S, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
            pair = (lambda q: q[:, 1], lambda q: q[:, 1], -1.0)
            gmsh.model.mesh.setPeriodic(1, groups["right"], groups["left"], aff)
        else:
            groups["far"] = groups["left"]
        h_w = np.pi * r_e / n_wire  # half circle -> n_wire segments
        _size_fields(gmsh, [(groups["emitter"], h_w), (groups["collector"], h_coll),
                            (groups["guard"], h_coll), (groups["insul"], 2 * h_coll)], h_max)
        gmsh.model.mesh.generate(2)
        info = dict(kind=mode, d=d, pitch=pitch, r_e=r_e, t_c=t_c, c_c=c_c, c_d=c_d,
                    x0=x0, x1=x1, stage_pitch=S, split=split)
        return _finish_gmsh(gmsh, groups, pair, 2, info)
    finally:
        gmsh.finalize()


# ------------------------------------------------------------------------- solver
class CoronaProblem:
    """Discrete operators for one mesh; solves the unit-voltage problem for a given current."""

    def __init__(self, mesh: Mesh):
        self.mesh = mesh
        p, t, rid, n = mesh.p, mesh.t, mesh.rid, mesh.n_red
        x = p[t, 0]
        y = p[t, 1]
        b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], axis=1)
        c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1)
        area = 0.5 * (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])
        if np.any(area <= 0):
            raise RuntimeError("inverted element")
        ke = (b[:, :, None] * b[:, None, :] + c[:, :, None] * c[:, None, :]) / (4.0 * area[:, None, None])
        rt = rid[t]
        ii = np.repeat(rt, 3, axis=1).ravel()
        jj = np.tile(rt, (1, 3)).ravel()
        self.K = sp.coo_matrix((ke.ravel(), (ii, jj)), shape=(n, n)).tocsr()
        self.M = np.bincount(rt.ravel(), weights=np.repeat(area / 3.0, 3), minlength=n)
        # quasi-periodic jump: K (phi + joff) -> extra load vector
        jo = mesh.joff[t]  # (M, 3)
        self.gjump = np.bincount(rt.ravel(), weights=np.einsum("eab,eb->ea", ke, jo).ravel(), minlength=n)

        # unique edges: summed K_ij, jump term and geometric offset (x_j - x_i, y_j - y_i)
        la = np.array([0, 1, 2])
        lb = np.array([1, 2, 0])
        ea = t[:, la].T.ravel()
        eb = t[:, lb].T.ravel()
        kv = ke[:, la, lb].T.ravel()
        ra, rb = rid[ea], rid[eb]
        dx = p[eb, 0] - p[ea, 0]
        dy = p[eb, 1] - p[ea, 1]
        dj = mesh.joff[ea] - mesh.joff[eb]  # (phi_a - phi_b) = phi[ra] - phi[rb] + dj
        swap = ra > rb
        ra, rb = np.where(swap, rb, ra), np.where(swap, ra, rb)
        dx = np.where(swap, -dx, dx)
        dy = np.where(swap, -dy, dy)
        dj = np.where(swap, -dj, dj)
        key = ra.astype(np.int64) * n + rb
        uk, inv = np.unique(key, return_inverse=True)
        self.ei = (uk // n).astype(np.int64)
        self.ej = (uk % n).astype(np.int64)
        self.eK = np.bincount(inv, weights=kv, minlength=len(uk))
        self.eKj = np.bincount(inv, weights=kv * dj, minlength=len(uk))
        cnt = np.bincount(inv, minlength=len(uk))
        self.edx = np.bincount(inv, weights=dx, minlength=len(uk)) / cnt
        self.edy = np.bincount(inv, weights=dy, minlength=len(uk)) / cnt

        self.em = np.unique(rid[mesh.emitter])
        self.gd = np.unique(rid[mesh.guard]) if len(mesh.guard) else np.array([], dtype=np.int64)
        self.co = np.unique(rid[mesh.collector])
        self.is_em = np.zeros(n, bool)
        self.is_em[self.em] = True
        self.dirichlet = np.zeros(n, bool)
        self.dirichlet[self.em] = True
        self.dirichlet[self.gd] = True
        self.dirichlet[self.co] = True
        self.hi = np.concatenate([self.em, self.gd])  # nodes held at emitter potential
        self.absorb = np.concatenate([self.co, self.gd])  # conductor surfaces that take up ions
        self.free = np.flatnonzero(~self.dirichlet)
        self.dir_idx = np.flatnonzero(self.dirichlet)
        self.Kff = self.K[self.free][:, self.free].tocsc()
        self.Kfd = self.K[self.free][:, self.dir_idx].tocsr()
        self.lu = spla.splu(self.Kff)
        self.len_em = self._boundary_length(mesh.emitter)
        self.len_co_nodes = self._boundary_node_lengths(mesh.collector)
        self.n = n

    # --- geometry helpers
    def _boundary_edges(self, nodes):
        m = self.mesh
        mark = np.zeros(len(m.p), bool)
        mark[nodes] = True
        t = m.t
        out = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            sel = mark[t[:, a]] & mark[t[:, b]]
            out.append(np.stack([t[sel, a], t[sel, b]], axis=1))
        e = np.concatenate(out)
        key = np.sort(e, axis=1)
        _, idx, cnt = np.unique(key, axis=0, return_index=True, return_counts=True)
        return key[idx[cnt == 1]]  # true boundary edges appear in exactly one triangle

    def _boundary_length(self, nodes):
        e = self._boundary_edges(nodes)
        return float(np.linalg.norm(self.mesh.p[e[:, 0]] - self.mesh.p[e[:, 1]], axis=1).sum())

    def _boundary_node_lengths(self, nodes):
        e = self._boundary_edges(nodes)
        L = np.linalg.norm(self.mesh.p[e[:, 0]] - self.mesh.p[e[:, 1]], axis=1)
        out = np.zeros(self.mesh.n_red)
        np.add.at(out, self.mesh.rid[e[:, 0]], 0.5 * L)
        np.add.at(out, self.mesh.rid[e[:, 1]], 0.5 * L)
        return out

    # --- field solves
    def solve_phi(self, rho, v_em=1.0):
        phi = np.zeros(self.n)
        phi[self.hi] = v_em
        rhs = (self.M * rho / EPS0 - v_em * self.gjump)[self.free] - self.Kfd @ phi[self.dir_idx]
        phi[self.free] = self.lu.solve(rhs)
        return phi

    def edge_flux(self, phi, v_em=1.0):
        """Electric flux Phi_ij leaving node i's dual cell towards node j (per unit depth)."""
        return -(self.eK * (phi[self.ei] - phi[self.ej]) + v_em * self.eKj)

    def node_outflux(self, flux):
        out = np.bincount(self.ei, weights=flux, minlength=self.n)
        out -= np.bincount(self.ej, weights=flux, minlength=self.n)
        return out

    def emitter_state(self, phi, rho):
        """(total outgoing electric flux, mean surface field) of the emitter, unit voltage."""
        flux = self.edge_flux(phi)
        fp = np.maximum(flux, 0.0)
        fm = np.maximum(-flux, 0.0)
        out = (np.bincount(self.ei, weights=fp, minlength=self.n)
               + np.bincount(self.ej, weights=fm, minlength=self.n))[self.em].sum()
        e = (self.node_outflux(flux)[self.em] - (self.M * rho / EPS0)[self.em]).sum() / self.len_em
        return out, e

    def solve_rho(self, phi, rho_old, rho_e):
        """Upwind transport solve for rho with density rho_e injected at the emitter."""
        flux = self.edge_flux(phi)
        fp = np.maximum(flux, 0.0)  # i -> j
        fm = np.maximum(-flux, 0.0)  # j -> i
        diag = np.bincount(self.ei, weights=fp, minlength=self.n) + np.bincount(self.ej, weights=fm, minlength=self.n)
        fb = self.M * rho_old / EPS0 - self.node_outflux(flux)  # flux out through electrode surfaces
        diag[self.absorb] += np.maximum(fb[self.absorb], 0.0)
        diag += 1e-12 * (np.abs(flux).max() + 1e-300)
        rows = np.concatenate([self.ej, self.ei])
        cols = np.concatenate([self.ei, self.ej])
        vals = np.concatenate([-fp, -fm])
        keep = ~self.is_em[rows]
        d = diag.copy()
        d[self.em] = 1.0
        A = sp.coo_matrix((vals[keep], (rows[keep], cols[keep])), shape=(self.n, self.n)).tocsr() + sp.diags(d)
        rhs = np.zeros(self.n)
        rhs[self.em] = rho_e
        return np.maximum(spla.spsolve(A.tocsc(), rhs), 0.0)

    def diagnostics(self, phi, rho):
        """Unit-voltage, unit-mobility integrals for the *full* physical unit (all symmetry copies)."""
        s = self.mesh.symmetry
        flux = self.edge_flux(phi)
        j = np.where(flux > 0, flux * rho[self.ei], flux * rho[self.ej])  # charge flux i->j (/mu)
        i_em = self.node_outflux(j)[self.em].sum()
        fb = self.M * rho / EPS0 - self.node_outflux(flux)
        i_co = (np.maximum(fb[self.co], 0.0) * rho[self.co]).sum()
        # ions that drift back onto conductors at emitter potential: no supply current, no net work
        i_gd = (np.maximum(fb[self.gd], 0.0) * rho[self.gd]).sum() if len(self.gd) else 0.0
        fx = (j * self.edx).sum()  # = integral(rho E_x dA) exactly, by telescoping
        _, e_em = self.emitter_state(phi, rho)
        with np.errstate(divide="ignore", invalid="ignore"):
            e_co = np.where(self.len_co_nodes[self.co] > 0, fb[self.co] / self.len_co_nodes[self.co], 0.0)
        return dict(i_em=s * i_em, i_co=s * i_co, i_gd=s * i_gd, fx=s * fx, e_em=float(e_em),
                    e_co_max=float(e_co.max()), e_co=e_co)

    def solve_unit(self, i_target, rho0=None, tol=1e-5, max_iter=600, relax=0.5):
        """Self-consistent (phi, rho) at 1 V on the emitter carrying total ion current i_target.

        i_target is the current of the full physical unit divided by the mobility (C/m^3 at
        1 V).  Driving the problem with the current rather than the injected density keeps
        the Picard iteration well conditioned all the way to the space-charge limit.
        Charge is injected at uniform density over the exposed emitter surface, and
        Kaptzov's condition is applied to the mean surface field.
        """
        s = self.mesh.symmetry
        rho = np.zeros(self.n) if rho0 is None else rho0.copy()
        phi = self.solve_phi(rho)
        if i_target <= 0:
            return phi, rho, self.diagnostics(phi, rho), 0
        last, last_err, dg = None, np.inf, None
        for it in range(1, max_iter + 1):
            out, _ = self.emitter_state(phi, rho)
            if out <= 0:
                return phi, rho, dict(self.diagnostics(phi, rho), e_em=-1.0), -it
            rho_new = self.solve_rho(phi, rho, i_target / (s * out))
            rho = (1 - relax) * rho + relax * rho_new
            phi = self.solve_phi(rho)
            dg = self.diagnostics(phi, rho)
            cur = np.array([dg["i_co"], dg["fx"], dg["e_em"]])
            if last is not None:
                err = np.max(np.abs(cur - last) / (np.abs(cur) + 1e-300))
                if err < tol * relax and abs(dg["i_em"] - dg["i_co"] - dg["i_gd"]) < 1e-2 * abs(dg["i_em"]):
                    return phi, rho, dg, it
                if err > last_err * 1.02 and it > 5:
                    relax = max(0.08, 0.8 * relax)
                last_err = err
            last = cur
        return phi, rho, dg, -max_iter

    def nodal_fields(self, phi, rho):
        """phi, rho, E on the *original* mesh nodes (for plotting), at unit voltage."""
        m = self.mesh
        ph = phi[m.rid] + m.joff
        p, t = m.p, m.t
        x, y = p[t, 0], p[t, 1]
        b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], axis=1)
        c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], axis=1)
        area = 0.5 * (b[:, 0] * c[:, 1] - b[:, 1] * c[:, 0])
        ex = -(b * ph[t]).sum(1) / (2 * area)
        ey = -(c * ph[t]).sum(1) / (2 * area)
        w = np.bincount(t.ravel(), weights=np.repeat(area, 3), minlength=len(p))
        exn = np.bincount(t.ravel(), weights=np.repeat(ex * area, 3), minlength=len(p)) / w
        eyn = np.bincount(t.ravel(), weights=np.repeat(ey * area, 3), minlength=len(p)) / w
        return ph, rho[m.rid], exn, eyn


@dataclass
class StageCurve:
    """I-V-thrust characteristic of one electrode pair (per metre of electrode span)."""

    V: np.ndarray  # applied voltage [V]
    I: np.ndarray  # current per unit span [A/m]
    F: np.ndarray  # EHD body force (ideal thrust) per unit span [N/m]
    E_coll: np.ndarray  # peak field on the collector surface [V/m]
    V_onset: float
    meta: dict

    def at(self, V):
        """Interpolate (I, F, E_coll) at voltage V (zero below onset)."""
        V = np.asarray(V, float)
        I = np.interp(V, self.V, self.I, left=0.0)
        F = np.interp(V, self.V, self.F, left=0.0)
        Ec = np.interp(V, self.V, self.E_coll)
        return I, F, Ec

    def to_dict(self):
        meta = {k: (v.item() if isinstance(v, np.generic) else v) for k, v in self.meta.items()}
        return dict(V=self.V.tolist(), I=self.I.tolist(), F=self.F.tolist(),
                    E_coll=self.E_coll.tolist(), V_onset=float(self.V_onset), meta=meta)

    @staticmethod
    def from_dict(d):
        return StageCurve(V=np.array(d["V"]), I=np.array(d["I"]), F=np.array(d["F"]),
                          E_coll=np.array(d["E_coll"]), V_onset=d["V_onset"], meta=d["meta"])


def stage_curve(prob: CoronaProblem, E_c: float, mu: float = MU_ION, n_pts: int = 14,
                v_max_over_onset: float = 7.0, verbose: bool = False) -> StageCurve:
    """Sweep the ion current and convert each unit-voltage solve to a point on the I-V curve.

    At 1 V with normalised current i the mean emitter field is f(i) [1/m].  The real
    discharge pins that field at E_c, so V = E_c / f, I = mu * i * V^2, F = fx * V^2.
    """
    phi, rho, d0, _ = prob.solve_unit(0.0)
    f0 = d0["e_em"]  # mean emitter field per volt (Laplace)
    v_on = E_c / f0
    V, I, F, EC = [v_on], [0.0], [0.0], [d0["e_co_max"] * v_on]
    info = prob.mesh.info
    if "d" in info:  # planar space-charge limit sets the current scale
        i_ref = 9.0 / 8.0 * EPS0 * info["pitch"] / info["d"] ** 3
    else:
        i_ref = 2 * np.pi * EPS0 / info["R"] ** 2
    i = 0.004 * i_ref
    f_target = np.geomspace(1.0, 1.0 / v_max_over_onset, n_pts + 1)[1:]  # wanted f/f0 values
    rho, pts, k, fails = None, [(0.0, 1.0)], 0, 0
    while k < len(f_target) and fails < 10:
        phi, rho_new, dg, it = prob.solve_unit(i, rho0=rho)
        f = dg["e_em"] / f0
        if it < 0 or f <= 0.5 * f_target[-1]:
            fails += 1  # too close to (or beyond) the space-charge limit: step back
            if verbose:
                print(f"  (no convergence at i/iref={i/i_ref:.4f}, it={it}, f/f0={f:.3f})")
            i = pts[-1][0] + 0.4 * (i - pts[-1][0])
            continue
        rho = rho_new
        pts.append((i, f))
        v = v_on / f
        V.append(v)
        I.append(mu * dg["i_co"] * v * v)
        F.append(dg["fx"] * v * v)
        EC.append(dg["e_co_max"] * v)
        if verbose:
            print(f"  i/iref={i/i_ref:7.4f} it={it:4d} V={v/1e3:6.2f} kV I={I[-1]*1e3:7.4f} mA/m "
                  f"F={F[-1]*1e3:7.3f} mN/m  F/(I d/mu)={F[-1]*mu/max(I[-1],1e-300)/info.get('d',1.0):.3f}"
                  f" Ecoll={EC[-1]/1e6:5.2f} MV/m leak={dg['i_gd']/max(dg['i_em'],1e-300):.4f}")
        while k < len(f_target) and f <= f_target[k]:
            k += 1
        if k >= len(f_target):
            break
        (i1, f1), (i2, f2) = pts[-2], pts[-1]
        slope = (i2 - i1) / (f2 - f1) if f2 != f1 else -i_ref
        i = float(np.clip(i2 + slope * (f_target[k] - f2), 1.05 * i2, 3.0 * i2))
    order = np.argsort(V)
    return StageCurve(V=np.array(V)[order], I=np.array(I)[order], F=np.array(F)[order],
                      E_coll=np.array(EC)[order], V_onset=v_on,
                      meta=dict(info, E_c=E_c, mu=mu, nodes=prob.n, symmetry=prob.mesh.symmetry))


# ------------------------------------------------------------- analytic reference
def coaxial_analytic(r0: float, R: float, E_c: float, I: np.ndarray, mu: float = MU_ION) -> np.ndarray:
    """Townsend's exact solution: voltage needed to drive current I [A/m] in a coaxial corona."""
    from scipy.integrate import quad

    out = []
    for i in np.atleast_1d(I):
        a = i / (2 * np.pi * EPS0 * mu)
        out.append(quad(lambda r: np.sqrt(a * (1 - (r0 / r) ** 2) + (E_c * r0 / r) ** 2), r0, R,
                        limit=400, points=np.geomspace(r0, R, 12)[1:-1])[0])
    return np.array(out)
