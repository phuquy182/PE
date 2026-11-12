# -*- coding: utf-8 -*-
"""
model.py — Runs the thirteen steps of Fig. 3.1 (section 3.7.1).
The operating clearance P_d,op is kept constant (section 3.1.3); it is not updated
with the interference.

Top-level module: takes a bearing, a fit class and a load, and returns all
intermediate quantities together with the fit life factor f_fit.
"""
import math
from dataclasses import dataclass

from . import iso286, hertz, loaddist, fitmodel, life
from .geometry import Bearing


@dataclass(frozen=True)
class Params:
    """Material constants and model parameters, as in Table 4.1."""
    E: float = 205878.0        # Young's modulus of bearing steel, MPa
    nu: float = 0.30           # Poisson's ratio
    roughness: float = 4.0     # interference loss from surface roughness, both surfaces, µm
    c_exp: float = 9.0         # shear-stress – life exponent
    e_weib: float = 1.11       # Weibull slope
    # Lundberg – Palmgren constant A in the raceway capacity Q_c.
    # A appears in BOTH Q_ci and Q_ce and the model only uses the ratio
    # X = (Q_ce/Q_ee)^3 / (Q_ci/Q_ei)^3, so A cancels and affects no result.
    # The value below only gives Q_c the right magnitude when printed.
    # This is the classical Lundberg – Palmgren Q_c, NOT equations (19)/(20) of
    # ISO/TS 16281, which derive Q_ci, Q_ce from the catalogue C_r.
    A_LP: float = 98.1         # N/mm^1.8
    sigma_res: float = 0.0     # residual stress, MPa (negative when compressive)

    def __post_init__(self):
        if not (self.E > 0 and 0 <= self.nu < 0.5):
            raise ValueError("Need E > 0 and 0 <= nu < 0.5.")
        if self.roughness < 0 or self.c_exp <= 0 or self.e_weib <= 0:
            raise ValueError("Need roughness >= 0, c > 0 and e > 0.")


@dataclass
class Result:
    """Full result of one model evaluation."""
    # --- input ---
    fit: str              # fit class, or "δ = … µm" when the interference is given directly
    Fr: float             # radial load, N
    delta_min: float      # smallest geometric interference, µm
    delta_max: float      # largest geometric interference, µm
    delta_mean: float     # mean geometric interference, µm
    delta_eff: float      # effective interference, µm
    # --- fit branch ---
    D_eff: float          # effective outer diameter of the inner ring, mm
    p_f: float            # fit pressure, MPa
    sigma_theta: float    # hoop stress at depth z0 — USED IN THE FATIGUE CORRECTION, MPa
    sigma_bore: float     # MAXIMUM hoop stress at the bore — equation (3.11), Table 4.3, MPa
    # --- load and contact branch ---
    delta_r: float        # radial displacement of the inner ring, mm
    Qmax: float           # maximum ball load, N
    a: float              # semi-major axis of the contact ellipse, mm
    b: float              # semi-minor axis of the contact ellipse, mm
    b_over_a: float       # ellipse ratio
    p_H: float            # maximum Hertz pressure, MPa
    tau_H: float          # maximum subsurface shear stress, MPa
    z0: float             # depth of the critical point, mm
    # --- combined ---
    tau_h: float          # corrected shear stress, MPa
    LR: float             # inner-race life ratio
    X: float              # outer-to-inner race life ratio
    f_fit: float          # fit life factor
    L10: float            # basic rating life, million revolutions
    L_fit: float          # life including the fit, million revolutions


def _check_load(Fr: float):
    if not (isinstance(Fr, (int, float)) and math.isfinite(Fr) and Fr > 0):
        raise ValueError(f"Radial load Fr must be positive, got {Fr!r}.")


def _branch_load(brg: Bearing, Fr: float, prm: Params):
    """Load and rolling-contact branch — steps 5 to 9."""
    _check_load(Fr)
    Ep = hertz.equivalent_modulus(prm.E, prm.nu)
    cg_i = hertz.contact_geometry(brg.Dw, brg.D_IR, brg.fi, inner=True)
    cg_e = hertz.contact_geometry(brg.Dw, brg.dm + brg.Dw, brg.fe, inner=False)
    Ki = hertz.contact_stiffness(cg_i, Ep)
    Ke = hertz.contact_stiffness(cg_e, Ep)
    Kn = loaddist.total_stiffness(Ki, Ke)

    # step 7: solve the equilibrium for the ball load distribution
    delta_r = loaddist.solve_delta_r(Fr, brg.Z, brg.Pd, Kn)
    loads = loaddist.ball_loads(delta_r, brg.Z, brg.Pd, Kn)
    Qmax = max(loads)

    # steps 8 and 9: contact ellipse, Hertz pressure, maximum shear stress and its depth
    pH, tau, z0, a, b = hertz.subsurface_shear(Qmax, cg_i, Ep)

    # for step 13: life ratio of the two races
    Qei, Qee = loaddist.equivalent_ball_loads(loads, brg.Z)
    Qc_i = life.Qc_raceway(prm.A_LP, brg.fi, brg.gamma, brg.Dw, brg.Z, inner=True)
    Qc_e = life.Qc_raceway(prm.A_LP, brg.fe, brg.gamma, brg.Dw, brg.Z, inner=False)
    X = life.race_life_ratio_X(Qc_i, Qc_e, Qei, Qee)

    u = z0 / b                                     # dimensionless depth, equation (3.25)

    return dict(delta_r=delta_r, loads=loads, Qmax=Qmax, p_H=pH, tau_H=tau,
                z0=z0, u=u, a=a, b=b, X=X, Kn=Kn, Ki=Ki, Ke=Ke, cg_i=cg_i, cg_e=cg_e)


def evaluate(brg: Bearing, fit: str = None, Fr: float = None, prm: Params = None,
             delta_geo_um: float = None) -> Result:
    """
    Run the thirteen steps for one fit class and one load.

    Give the interference in exactly ONE of two ways:
      fit          : ISO 286 fit class, e.g. 'm6'; evaluated at the mean interference.
      delta_geo_um : geometric interference in µm, given directly — for parameter
                     sweeps and the verification cases; then leave fit = None.
    """
    prm = prm or Params()
    _check_load(Fr)

    # --- steps 1 and 2: geometric and effective interference ---
    if (fit is None) == (delta_geo_um is None):
        raise ValueError("Give exactly one of: fit (fit class) or delta_geo_um (interference, µm).")
    if delta_geo_um is None:
        itf = iso286.interference(brg.d, fit)
        dmin, dmax, dmean = itf['dmin'], itf['dmax'], itf['dmean']
    else:
        if not math.isfinite(delta_geo_um):
            raise ValueError(f"Invalid interference: {delta_geo_um!r}.")
        dmin = dmax = dmean = float(delta_geo_um)
        fit = f'δ = {dmean:g} µm'
    delta_eff = fitmodel.effective_interference(dmean, prm.roughness)

    # --- steps 3 and 4: effective diameter and fit pressure ---
    D_eff = brg.D_eff
    p_f = fitmodel.fit_pressure(delta_eff, brg.d, D_eff, prm.E)

    # --- steps 5 to 9: load and contact branch ---
    L = _branch_load(brg, Fr, prm)

    # --- step 10: hoop stress at the location of the maximum shear stress ---
    # Must follow step 9 because it needs the depth z0 from the load branch.
    k3 = fitmodel.k3_parameter(prm.E, prm.nu, brg.D_IR, brg.Dw, L['p_H'])
    sigma_theta = fitmodel.hoop_stress_at_depth(p_f, brg.d, D_eff, brg.D_IR, L['u'], k3)
    # Maximum stress at the bore, equation (3.11), Table 4.3 — not used for fatigue.
    sigma_bore = fitmodel.hoop_stress_at_bore(p_f, brg.d, D_eff)

    # --- steps 11 to 13: superposition, life ratio and f_fit ---
    tau_h = life.corrected_shear(L['tau_H'], sigma_theta, prm.sigma_res)
    LR = life.life_ratio(L['tau_H'], tau_h, prm.c_exp)
    f_fit = life.f_fit_factor(LR, L['X'], prm.e_weib)

    L10 = life.L10_iso281(brg.C, Fr)

    return Result(fit=fit, Fr=Fr,
                  delta_min=dmin, delta_max=dmax, delta_mean=dmean,
                  delta_eff=delta_eff * 1000.0,
                  D_eff=D_eff, p_f=p_f, sigma_theta=sigma_theta, sigma_bore=sigma_bore,
                  delta_r=L['delta_r'], Qmax=L['Qmax'],
                  a=L['a'], b=L['b'], b_over_a=L['b'] / L['a'],
                  p_H=L['p_H'], tau_H=L['tau_H'], z0=L['z0'],
                  tau_h=tau_h, LR=LR, X=L['X'], f_fit=f_fit,
                  L10=L10, L_fit=f_fit * L10)


def life_hours(brg: Bearing, Fr: float, delta_geo_um: float, n_rpm: float,
               prm: Params = None) -> float:
    """Life including the fit, L_fit,h = f_fit · L10h, hours — full thirteen-step chain."""
    r = evaluate(brg, None, Fr, prm, delta_geo_um=delta_geo_um)
    return life.million_rev_to_hours(r.L_fit, n_rpm)


def closed_form_constants(brg: Bearing, Fr: float, prm: Params = None) -> dict:
    """
    Constants of the closed-form life function (3.37) at one load: K_δ (mm⁻¹) and q = 2X^(−e).

    With F_r and P_d,op fixed, τ_H, D_z and X do not depend on the interference and
    σ_θ is linear in δ_eff, so σ_θ/(2|τ_H|) = K_δ·δ_eff holds exactly (for σ_res = 0).
    """
    prm = prm or Params()
    if prm.sigma_res != 0.0:
        raise ValueError("The closed form (3.37) requires sigma_res = 0.")
    L = _branch_load(brg, Fr, prm)
    k3 = fitmodel.k3_parameter(prm.E, prm.nu, brg.D_IR, brg.Dw, L['p_H'])
    Dz = fitmodel.depth_diameter(brg.D_IR, L['u'], k3)
    Kd = (prm.E * brg.d / (4.0 * abs(L['tau_H']))) * (1.0 / brg.D_eff ** 2 + 1.0 / Dz ** 2)
    q = 2.0 * L['X'] ** (-prm.e_weib)
    return dict(Kd=Kd, q=q, Dz=Dz, k3=k3, X=L['X'], tau_H=L['tau_H'])


def delta_max_life(brg: Bearing, Fr: float, L_req_h: float, n_rpm: float,
                   prm: Params = None) -> float:
    """
    Upper limit of the geometric interference set by the required life, µm — equation (3.43).

      · L10h < L_req_h : no interference meets the requirement → ValueError;
      · otherwise      : unique solution, equal to the thirteen-step chain to round-off.
    """
    prm = prm or Params()
    if L_req_h <= 0 or n_rpm <= 0:
        raise ValueError("Need L_req_h > 0 and n_rpm > 0.")
    L10h = life.million_rev_to_hours(life.L10_iso281(brg.C, Fr), n_rpm)
    if L10h < L_req_h:
        raise ValueError(f"L10h = {L10h:.0f} h < L_req = {L_req_h:.0f} h: no fit class can comply.")
    k = closed_form_constants(brg, Fr, prm)
    e, c, q = prm.e_weib, prm.c_exp, k['q']
    base = (1.0 + q) * (L10h / L_req_h) ** e - q
    return prm.roughness + (base ** (1.0 / (c * e)) - 1.0) / k['Kd'] * 1000.0


def ntn_min_interference(brg: Bearing, Fr: float) -> float:
    """
    Minimum effective interference against creep per NTN, µm — equations (3.40), (3.41).
    The two formulas switch at Fr = 0.3·C0r.
    """
    _check_load(Fr)
    if Fr <= 0.3 * brg.C0:
        return 0.08 * (brg.d * Fr / brg.B) ** 0.5
    return 0.02 * Fr / brg.B
