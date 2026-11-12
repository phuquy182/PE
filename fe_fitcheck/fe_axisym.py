# -*- coding: utf-8 -*-
"""
fe_axisym.py — Axisymmetric linear-elastic FE model of a shaft – inner-ring
interference fit, used to check the fit branch of bearing_fit (Paper 2, task 3).

Model
-----
* Coordinates (r, z); half model z >= 0 (mid-plane symmetry).
* Two bodies (solid shaft, inner ring), conforming meshes on the fit surface.
* Interference: u_r(ring) - u_r(shaft) = delta_eff/2 at every radial DOF of the
  fit surface, imposed with Lagrange multipliers; u_z of the two bodies is
  independent (frictionless). Closed contact is assumed over the full width and
  checked afterwards (all multipliers must be compressive).
* BCs: u_r = 0 on r = 0 (shaft), u_z = 0 on z = 0 (both bodies).
* Quadratic triangles (P2), r-weighted weak form without the factor 2*pi.
  Hence sum(nodal contact forces) = p_mean * (d/2) * (B/2).
* Hoop stress sigma_theta = lambda*tr(eps) + 2*mu*eps_theta, L2-projected onto a
  continuous P2 field per body before point evaluation.

Units: mm, N, MPa.
"""
import math
from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import triangle
from skfem import (MeshTri, Basis, ElementTriP2, ElementVector,
                   BilinearForm, LinearForm)

# ----------------------------------------------------------------------------- geometry


@dataclass
class Section:
    """Axial half-section (z >= 0) of the shaft – ring assembly."""
    a: float                 # bore radius d/2
    half_B: float            # half ring width B/2
    L_shaft: float           # shaft half length (>= half_B)
    ring_outline: list       # closed polyline of the ring, counter-clockwise, list of (r, z)
    r_groove: float          # groove-bottom radius (None for the plain ring)
    r_outer_mid: float       # outer radius of the ring at z = 0
    focus: list = field(default_factory=list)   # points that need the finest mesh


def _arc_points(rc, R, r_bottom, z_mouth, ds):
    """Points on the groove arc (centre (rc, 0), radius R) from the bottom (r_bottom, 0)
    up to the mouth at z = z_mouth; spacing about ds."""
    phi_end = math.asin(z_mouth / R)           # angle from the -r direction
    n = max(8, int(math.ceil(R * phi_end / ds)))
    phis = np.linspace(0.0, phi_end, n + 1)
    return [(rc - R * math.cos(p), R * math.sin(p)) for p in phis]


def plain_ring_section(d, B, D_out, L_shaft):
    """Case A: plain ring with outer diameter D_out."""
    a, hB, ro = 0.5 * d, 0.5 * B, 0.5 * D_out
    ring = [(a, 0.0), (ro, 0.0), (ro, hB), (a, hB)]
    return Section(a, hB, L_shaft, ring, None, ro, focus=[(ro, 0.0)])


def grooved_ring_section(d, B, D_IR, Dw, fi=0.52, shoulder=0.2, L_shaft=None, ds_arc=0.04):
    """
    Cases B – D: inner ring with the raceway groove.
      bore radius d/2, groove bottom at D_IR/2 on the mid-plane,
      groove arc radius fi*Dw, shoulder radius D_IR/2 + shoulder*Dw, no chamfers.
    """
    a, hB = 0.5 * d, 0.5 * B
    R = fi * Dw
    rb = 0.5 * D_IR
    h = shoulder * Dw
    rs = rb + h
    z_m = math.sqrt(R * R - (R - h) ** 2)       # half width of the groove mouth
    if z_m >= hB:
        raise ValueError("Groove mouth wider than the ring.")
    arc = _arc_points(rb + R, R, rb, z_m, ds_arc)
    arc[-1] = (rs, z_m)
    ring = [(a, 0.0)] + arc + [(rs, hB), (a, hB)]
    return Section(a, hB, L_shaft if L_shaft is not None else hB, ring, rb, rb,
                   focus=[(rb, 0.0)])


def ring_area(sec: Section) -> float:
    """Area of the ring half-section (shoelace)."""
    p = np.array(sec.ring_outline)
    x, y = p[:, 0], p[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


# ----------------------------------------------------------------------------- mesh


def _size_function(sec: Section, h_glob, h_fine, h_edge, grade=0.25):
    """Target element size: h_fine at the focus points, h_edge at the bore and the
    end of the fit surface, growing linearly with distance, capped at h_glob."""
    a, hB = sec.a, sec.half_B
    fx = np.array(sec.focus)

    def h(x, y):
        d_f = np.min(np.hypot(x[:, None] - fx[None, :, 0], y[:, None] - fx[None, :, 1]), axis=1)
        hh = np.minimum(h_glob, h_fine + grade * d_f)
        # fit surface (bore line r = a, 0 <= z <= hB) and the corner (a, hB)
        d_bore = np.hypot(x - a, np.clip(y - hB, 0.0, None))
        hh = np.minimum(hh, h_edge + grade * d_bore)
        d_corner = np.hypot(x - a, y - hB)
        hh = np.minimum(hh, 0.5 * h_edge + grade * d_corner)
        return hh
    return h


def build_mesh(sec: Section, h_glob, h_fine, h_edge, max_iter=12):
    """
    Conforming triangulation of shaft (region 1) and ring (region 2).
    Returns (points (N,2), triangles (M,3), region (M,)).
    """
    a, hB, Ls = sec.a, sec.half_B, sec.L_shaft
    hfun = _size_function(sec, h_glob, h_fine, h_edge)

    # densify straight boundary lines so that boundary resolution follows h
    def line(p, q):
        p, q = np.array(p, float), np.array(q, float)
        L = np.linalg.norm(q - p)
        n = max(1, int(math.ceil(L / max(1e-9, 0.7 * float(np.min(hfun(
            np.linspace(p[0], q[0], 50), np.linspace(p[1], q[1], 50))))))))
        return [tuple(p + (q - p) * t) for t in np.linspace(0, 1, n + 1)]

    verts, segs = [], []
    index = {}

    def vid(pt):
        key = (round(pt[0], 9), round(pt[1], 9))
        if key not in index:
            index[key] = len(verts)
            verts.append(key)
        return index[key]

    def add_poly(pts):
        for p, q in zip(pts[:-1], pts[1:]):
            i, j = vid(p), vid(q)
            if i != j:
                segs.append((i, j))

    # shaft outline (counter-clockwise): (0,0) -> (a,0) -> (a,hB) -> (a,Ls) -> (0,Ls) -> (0,0)
    add_poly(line((0, 0), (a, 0)))
    add_poly(line((a, 0), (a, hB)))       # fit surface (shared)
    if Ls > hB + 1e-9:
        add_poly(line((a, hB), (a, Ls)))
    add_poly(line((a, Ls), (0, Ls)))
    add_poly(line((0, Ls), (0, 0)))
    # ring outline without the fit surface (already added)
    ro = sec.ring_outline
    ring_pts = []
    for p, q in zip(ro[:-1], ro[1:]):
        if abs(p[0] - a) < 1e-12 and abs(q[0] - a) < 1e-12:
            continue
        seg = line(p, q) if math.hypot(q[0] - p[0], q[1] - p[1]) > 0.2 else [p, q]
        ring_pts.append(seg)
    for s in ring_pts:
        add_poly(s)
    segs = list({tuple(sorted(s)) for s in segs})

    pslg = dict(vertices=np.array(verts), segments=np.array(segs),
                regions=np.array([[0.5 * a, 0.5 * min(hB, Ls), 1, 0],
                                  [a + 1e-3 * (sec.r_outer_mid - a) + 0.5 * (sec.r_outer_mid - a),
                                   0.5 * hB, 2, 0]]))
    t = triangle.triangulate(pslg, 'pq30Aa%.6f' % (0.433 * h_glob ** 2))
    for _ in range(max_iter):
        P, T = t['vertices'], t['triangles']
        c = P[T].mean(axis=1)
        target = 0.433 * hfun(c[:, 0], c[:, 1]) ** 2
        e1, e2 = P[T[:, 1]] - P[T[:, 0]], P[T[:, 2]] - P[T[:, 0]]
        area = 0.5 * np.abs(e1[:, 0] * e2[:, 1] - e1[:, 1] * e2[:, 0])
        if np.all(area <= 1.05 * target):
            break
        t['triangle_max_area'] = target
        t = triangle.triangulate(t, 'rpq30Aa')
    reg = t['triangle_attributes'][:, 0].astype(int)
    return t['vertices'], t['triangles'], reg


def _submesh(P, T, mask):
    tri = T[mask]
    used = np.unique(tri)
    remap = -np.ones(len(P), int)
    remap[used] = np.arange(len(used))
    return MeshTri(P[used].T.copy(), remap[tri].T.copy())


# ----------------------------------------------------------------------------- solver


@dataclass
class FEResult:
    sec: Section
    delta: float              # effective diametral interference, mm
    n_dof: int
    p_mean: float             # mean fit pressure from the multipliers, MPa
    p_nodes: np.ndarray       # (z, nodal pressure) along the fit surface
    all_compressive: bool
    min_p_node: float
    sig_bore_mid: float       # sigma_theta at (a, 0), ring side
    sig_bore_avg: float       # width average of sigma_theta along the bore, ring side
    u_r_groove: float         # u_r at the groove bottom (or outer surface) at z = 0, mm
    ring: dict                # mesh, bases and fields of the ring (for plotting / probing)
    shaft: dict

    def sig_theta_at(self, r, z=0.0):
        return probe(self.ring, 'sig_t', r, z)

    def u_r_at(self, r, z=0.0):
        return probe(self.ring, 'u_r', r, z)


def probe(body, name, r, z):
    r = np.atleast_1d(np.asarray(r, float))
    z = np.broadcast_to(np.asarray(z, float), r.shape)
    pts = np.vstack([r, np.maximum(z, 1e-9)])
    # keep probe points strictly inside the ring section
    pts[0] = np.clip(pts[0], body['rmin'] + 1e-9, None)
    M = body['basis_s'].probes(pts)
    v = M @ body[name]
    return v if v.size > 1 else float(v[0])


def _elastic_forms(E, nu):
    lam = E * nu / ((1 + nu) * (1 - 2 * nu))
    mu = E / (2 * (1 + nu))

    def strains(u, r):
        return (u.grad[0][0], u.grad[1][1], u[0] / r, u.grad[0][1] + u.grad[1][0])

    @BilinearForm
    def stiffness(u, v, w):
        r = w.x[0]
        eu, ev = strains(u, r), strains(v, r)
        tru, trv = eu[0] + eu[1] + eu[2], ev[0] + ev[1] + ev[2]
        return (lam * tru * trv + 2 * mu * (eu[0] * ev[0] + eu[1] * ev[1] + eu[2] * ev[2])
                + mu * eu[3] * ev[3]) * r

    def sig_theta_qp(uh, r):
        err, ezz, ett = uh.grad[0][0], uh.grad[1][1], uh[0] / r
        return lam * (err + ezz + ett) + 2 * mu * ett

    return stiffness, sig_theta_qp


def solve(sec: Section, delta, E=205878.0, nu=0.30, h_glob=1.0, h_fine=0.05, h_edge=0.25,
          intorder=4):
    """Solve one load case. delta = effective diametral interference, mm."""
    P, T, reg = build_mesh(sec, h_glob, h_fine, h_edge)
    stiffness, sig_theta_qp = _elastic_forms(E, nu)
    e = ElementVector(ElementTriP2())
    bodies = {}
    for name, rid in (('shaft', 1), ('ring', 2)):
        m = _submesh(P, T, reg == rid)
        bv = Basis(m, e, intorder=intorder)
        bs = Basis(m, ElementTriP2(), intorder=intorder)
        bodies[name] = dict(mesh=m, basis=bv, basis_s=bs, K=stiffness.assemble(bv),
                            rmin=m.p[0].min())

    sh, rg = bodies['shaft'], bodies['ring']
    n_s, n_r = sh['basis'].N, rg['basis'].N
    K = sp.block_diag([sh['K'], rg['K']]).tocsr()
    N = n_s + n_r
    a, hB = sec.a, sec.half_B
    tol = 1e-7

    # radial DOFs on the fit surface, matched by location
    def fit_dofs(body):
        b = body['basis']
        loc = b.doflocs
        # component 0 (u_r) dofs: nodal_dofs[0], facet_dofs[0]
        cand = np.concatenate([b.nodal_dofs[0], b.facet_dofs[0]])
        x = loc[:, cand]
        sel = (np.abs(x[0] - a) < tol) & (x[1] <= hB + tol)
        return cand[sel], x[:, sel]

    ds, xs = fit_dofs(sh)
    dr, xr = fit_dofs(rg)
    order_s, order_r = np.lexsort((xs[0], xs[1])), np.lexsort((xr[0], xr[1]))
    ds, xs, dr, xr = ds[order_s], xs[:, order_s], dr[order_r], xr[:, order_r]
    if len(ds) != len(dr) or np.max(np.abs(xs - xr)) > 1e-8:
        raise RuntimeError("Fit-surface DOFs of shaft and ring do not match.")
    m_c = len(ds)
    rows = np.repeat(np.arange(m_c), 2)
    cols = np.ravel(np.column_stack([dr + n_s, ds]))
    vals = np.tile([1.0, -1.0], m_c)
    C = sp.csr_matrix((vals, (rows, cols)), shape=(m_c, N))
    g = np.full(m_c, 0.5 * delta)

    # Dirichlet DOFs (homogeneous)
    fixed = []
    for off, body, axis_too in ((0, sh, True), (n_s, rg, False)):
        b = body['basis']
        loc = b.doflocs
        for comp in (0, 1):
            cand = np.concatenate([b.nodal_dofs[comp], b.facet_dofs[comp]])
            x = loc[:, cand]
            if comp == 1:
                fixed.append(cand[np.abs(x[1]) < tol] + off)
            elif axis_too:
                fixed.append(cand[np.abs(x[0]) < tol] + off)
    fixed = np.unique(np.concatenate(fixed))
    free = np.setdiff1d(np.arange(N), fixed)

    A = sp.bmat([[K[free][:, free], C[:, free].T], [C[:, free], None]]).tocsc()
    rhs = np.concatenate([np.zeros(len(free)), g])
    sol = spla.spsolve(A, rhs)
    u = np.zeros(N)
    u[free] = sol[:len(free)]
    lam_c = sol[len(free):]
    res_c = np.max(np.abs(C @ u - g))
    if res_c > 1e-6 * delta:
        raise RuntimeError("Constraint not satisfied: residual %.3e mm." % res_c)

    # contact force on the ring dof = -lambda; positive = compressive
    f_c = -lam_c
    p_mean = f_c.sum() / (a * hB)
    # nodal pressures: force / integral of the basis function (r-weighted) over the
    # fit surface, from the same P2 trace (vertex: L/6, midside: 2L/3 on each edge)
    z_c = xr[1]
    zs = np.unique(np.round(z_c, 12))
    w = np.zeros(m_c)
    # fit-surface edges: consecutive vertex pairs; vertices are the even-indexed sorted points
    zv = z_c[0::2]
    for k in range(len(zv) - 1):
        L = zv[k + 1] - zv[k]
        w[2 * k] += a * L / 6
        w[2 * k + 2] += a * L / 6
        w[2 * k + 1] += a * 2 * L / 3
    p_nodes = np.column_stack([z_c, f_c / w])

    u_s, u_r = u[:n_s], u[n_s:]
    for body, ub in ((sh, u_s), (rg, u_r)):
        bv, bs = body['basis'], body['basis_s']
        uh = bv.interpolate(ub)

        @LinearForm
        def rhs_sig(v, w_):
            return sig_theta_qp(w_['uh'], w_.x[0]) * v

        @BilinearForm
        def mass(uu, vv, w_):
            return uu * vv

        M = mass.assemble(bs)
        bvec = rhs_sig.assemble(bs, uh=uh)
        body['sig_t'] = spla.spsolve(M.tocsc(), bvec)
        # u_r as a scalar P2 field (same DOF layout as the scalar basis)
        body['u_r'] = ub[bv.nodal_dofs[0]].tolist()
        ur_s = np.zeros(bs.N)
        ur_s[bs.nodal_dofs[0]] = ub[bv.nodal_dofs[0]]
        ur_s[bs.facet_dofs[0]] = ub[bv.facet_dofs[0]]
        body['u_r'] = ur_s
        body['u'] = ub

    sig_mid = probe(rg, 'sig_t', a, 0.0)
    zq = np.linspace(0.0, hB, 401)
    s_line = probe(rg, 'sig_t', np.full_like(zq, a), zq)
    sig_avg = np.trapezoid(s_line, zq) / hB
    ur_g = probe(rg, 'u_r', sec.r_outer_mid - 1e-9, 0.0)

    return FEResult(sec=sec, delta=delta, n_dof=N, p_mean=p_mean, p_nodes=p_nodes,
                    all_compressive=bool(np.all(f_c > 0)), min_p_node=float(np.min(f_c / w)),
                    sig_bore_mid=sig_mid, sig_bore_avg=sig_avg, u_r_groove=ur_g,
                    ring=rg, shaft=sh)


# ----------------------------------------------------------------------------- Lamé


def lame_pressure(delta, d, D_out, E):
    """Plane-stress Lamé fit pressure, same material, solid shaft (= fitmodel.fit_pressure)."""
    return (E * delta / d) * (D_out ** 2 - d ** 2) / (2.0 * D_out ** 2)


def lame_hoop(p, d, D_out, r):
    a2, b2 = (0.5 * d) ** 2, (0.5 * D_out) ** 2
    return p * a2 / (b2 - a2) * (1.0 + b2 / np.asarray(r) ** 2)


def lame_ur(p, d, D_out, r, E, nu):
    a2, b2 = (0.5 * d) ** 2, (0.5 * D_out) ** 2
    r = np.asarray(r)
    return r * p * a2 / (E * (b2 - a2)) * ((1 - nu) + (1 + nu) * b2 / r ** 2)
