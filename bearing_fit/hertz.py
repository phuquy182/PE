# -*- coding: utf-8 -*-
"""
hertz.py — Elliptical Hertz contact between ball and raceway.

Uses the Brewe – Hamrock approximations for the contact ellipse and the curve fits
of [9] for the maximum subsurface shear stress.
Steps 5, 6, 8 and 9 of Fig. 3.1, equations (3.15) – (3.25).
"""
import math
from dataclasses import dataclass

Q_CONST = math.pi / 2.0 - 1.0   # constant q = π/2 − 1 in the Brewe – Hamrock approximation


@dataclass
class ContactGeom:
    """Geometric quantities of one elliptical contact."""
    Rx: float        # equivalent radius in the rolling direction, mm
    Ry: float        # equivalent radius across the groove, mm
    R: float         # combined equivalent radius, mm
    alpha_r: float   # radius ratio Ry / Rx
    k: float         # ellipticity a / b
    eps: float       # complete elliptic integral of the second kind (approx.)
    Im: float        # complete elliptic integral of the first kind (approx.)


def contact_geometry(Dw: float, d_race: float, f: float, inner: bool) -> ContactGeom:
    """
    Ball – raceway contact geometry — equations (3.15), (3.17).

    Dw     : ball diameter, mm
    d_race : raceway (groove bottom) diameter, mm
    f      : raceway conformity, groove radius divided by ball diameter
    inner  : True for the inner-ring contact (convex), False for the outer ring (concave)
    """
    r_ball = 0.5 * Dw
    sgn = 1.0 if inner else -1.0
    inv_Rx = 1.0 / r_ball + sgn * 1.0 / (0.5 * d_race)
    inv_Ry = 1.0 / r_ball - 1.0 / (f * Dw)          # concave groove, hence the minus sign
    Rx, Ry = 1.0 / inv_Rx, 1.0 / inv_Ry
    alpha_r = Ry / Rx
    k = alpha_r ** (2.0 / math.pi)                    # ellipticity
    eps = 1.0 + Q_CONST / alpha_r                     # elliptic integral, second kind
    Im = math.pi / 2.0 + Q_CONST * math.log(alpha_r)  # elliptic integral, first kind
    return ContactGeom(Rx, Ry, 1.0 / (inv_Rx + inv_Ry), alpha_r, k, eps, Im)


def equivalent_modulus(E: float, nu: float) -> float:
    """Equivalent modulus E′ of two bodies of the same material — equation (3.16)."""
    return E / (1.0 - nu * nu)


def contact_stiffness(cg: ContactGeom, Ep: float) -> float:
    """
    Contact stiffness constant K of one contact, Q = K·η^1.5 — equations (3.18), (3.19).
    Unit N/mm^1.5. Derived from the Hamrock – Dowson deflection relation.
    """
    return (math.pi * cg.k * Ep / cg.Im ** 1.5) * math.sqrt(2.0 * cg.eps * cg.R / 9.0)


def ellipse_axes(Q: float, cg: ContactGeom, Ep: float):
    """Semi-major axis a and semi-minor axis b of the contact ellipse, mm — equations (3.20), (3.21)."""
    base = 6.0 * cg.eps * Q * cg.R / (math.pi * Ep)
    a = (base * cg.k * cg.k) ** (1.0 / 3.0)
    b = (base / cg.k) ** (1.0 / 3.0)
    return a, b


def hertz_pressure(Q: float, a: float, b: float) -> float:
    """Maximum Hertz pressure at the centre of the ellipse, MPa — equation (3.22)."""
    return 3.0 * Q / (2.0 * math.pi * a * b)


def tau_ratio(b_over_a: float) -> float:
    """Ratio τ_max / p_H, curve fit of [9] — equation (3.23)."""
    x = b_over_a
    return 0.1002 * x ** 3 - 0.2431 * x ** 2 + 0.1639 * x + 0.30


def depth_ratio(b_over_a: float) -> float:
    """Dimensionless depth u = z0 / b, curve fit of [9] — equation (3.25)."""
    x = b_over_a
    return 0.3108 * x ** 3 - 0.5131 * x ** 2 - 0.1159 * x + 0.786


def subsurface_shear(Q: float, cg: ContactGeom, Ep: float):
    """
    Stress state at the contact — equations (3.20) – (3.25).
    Returns (p_H, τ_max, z0, a, b).
    τ_max is negative, following the sign convention of [9].
    """
    a, b = ellipse_axes(Q, cg, Ep)
    pH = hertz_pressure(Q, a, b)
    kappa = b / a
    tau = -tau_ratio(kappa) * pH
    z0 = depth_ratio(kappa) * b
    return pH, tau, z0, a, b
