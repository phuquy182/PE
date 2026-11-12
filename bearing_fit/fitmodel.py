# -*- coding: utf-8 -*-
"""
fitmodel.py — Fit branch: from interference to fit pressure and hoop stress.

Based on Lamé's thick-walled cylinder solution, assuming a solid shaft of the same
material as the inner ring. Steps 2, 4 and 10 of Fig. 3.1,
equations (3.6), (3.9) – (3.11), (3.26) – (3.27).
"""


def effective_interference(delta_geo_um: float, roughness_um: float = 4.0) -> float:
    """
    Effective interference, mm — equation (3.6).

    delta_geo_um : geometric interference from ISO 286, µm
    roughness_um : loss of interference from asperity smoothing during press fitting,
                   for both mating surfaces together, µm

    When the geometric interference is below the roughness loss the fit has no
    effective interference; the function returns 0 and no fit pressure is produced.
    """
    return max(delta_geo_um - roughness_um, 0.0) / 1000.0


def fit_pressure(delta_eff: float, d: float, D_eff: float, E: float) -> float:
    """Contact pressure at the shaft – inner-ring interface, MPa — equation (3.9)."""
    if delta_eff <= 0.0:
        return 0.0
    return (E * delta_eff / d) * (D_eff ** 2 - d ** 2) / (2.0 * D_eff ** 2)


def hoop_stress(p_f: float, d: float, D_eff: float, D_z: float) -> float:
    """Tensile hoop stress in the inner ring at diameter D_z, MPa — equation (3.10)."""
    if p_f <= 0.0:
        return 0.0
    return p_f * d ** 2 / (D_eff ** 2 - d ** 2) * (1.0 + D_eff ** 2 / D_z ** 2)


def hoop_stress_at_bore(p_f: float, d: float, D_eff: float) -> float:
    """
    MAXIMUM hoop stress in the inner ring, at the bore surface D_z = d.

    Equation (3.11); reported in Table 4.3. The thesis does not check ring strength
    with this value (rolling contact fatigue only).

    Do not confuse with hoop_stress_at_depth(): that value, at the depth z0, is the
    one used in the fatigue correction of [9]. The two values differ and serve
    different purposes.

    Uses the equivalent ring D_eff of equation (3.8).
    """
    return hoop_stress(p_f, d, D_eff, d)


def k3_parameter(E: float, nu: float, D_IR: float, Dw: float, p_H: float) -> float:
    """
    Dimensionless parameter k3 locating the maximum shear stress — equation (3.27).

        k3 = E (R' + 1) / [4 (1 - nu^2) p_H]   with  R' = D_IR / Dw

    Original form of [9], step (5) of the worked example.
    NOTATION: [9] writes S_max for the maximum Hertz stress (p_H here) and d for the
    BALL DIAMETER; the thesis reserves d for the bore and uses Dw for the ball,
    so R' = D_IR/d of [9] becomes R' = D_IR/Dw.
    """
    R_prime = D_IR / Dw
    return E * (R_prime + 1.0) / (4.0 * (1.0 - nu * nu) * p_H)


def depth_diameter(D_IR: float, u: float, k3: float) -> float:
    """Diameter at the depth of maximum shear stress — equation (3.26): D_z = D_IR(1 - u/k3)."""
    return D_IR * (1.0 - u / k3)


def hoop_stress_at_depth(p_f: float, d: float, D_eff: float,
                         D_IR: float, u: float, k3: float) -> float:
    """
    Hoop stress at the location of maximum shear stress — equations (3.10), (3.26), (3.27).

    u  : dimensionless depth z0/b from equation (3.25)
    k3 : parameter from equation (3.27)
    """
    return hoop_stress(p_f, d, D_eff, depth_diameter(D_IR, u, k3))
