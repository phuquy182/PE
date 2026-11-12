# -*- coding: utf-8 -*-
"""
life.py — Stress superposition, inner-race life ratio, bearing system life and
fit life factor f_fit.

Steps 11, 12 and 13 of Fig. 3.1, equations (3.28) – (3.35).
"""


def corrected_shear(tau_H: float, sigma_theta: float, sigma_res: float = 0.0) -> float:
    """
    Maximum shear stress corrected for the hoop stress — equation (3.28).

    Sign convention: tau_H is negative, sigma_theta is positive (tensile),
    sigma_res is negative when compressive. Subtracting a positive amount from a
    negative number increases the magnitude: larger shear amplitude, shorter life.

    Source: the form with sigma_res is the extended equation of [22] (form (2.20));
    with sigma_res = 0 it reduces to equation (1a)/(1b) of [9] — see section 3.5.2.
    The thesis always uses sigma_res = 0.
    """
    return tau_H - 0.5 * (sigma_theta + sigma_res)


def life_ratio(tau_H: float, tau_h: float, c: float = 9.0) -> float:
    """Inner-race life ratio L_R = (τ_H / τ_h)^c — equation (3.29).
    c is the shear-stress – life exponent."""
    return (tau_H / tau_h) ** c


def Qc_raceway(A: float, f: float, gamma: float, Dw: float, Z: int, inner: bool) -> float:
    """
    Basic dynamic capacity of one raceway, Lundberg – Palmgren theory.

    Equation (11.67), section 11.5 of [2]; (3.31), (3.32) in the thesis:

        Q_c = A (2f/(2f-1))^0.41 * (1∓γ)^1.39 / (1±γ)^(1/3) * (γ/cos α)^0.3
              * D_w^1.8 * Z^(-1/3)

    upper sign for the inner ring, lower sign for the outer ring; α = 0 for a
    radial deep groove ball bearing, so cos α = 1.

    A : Lundberg – Palmgren material constant for through-hardened 52100 steel,
        A = 98.1 with D_w in mm and Q_c in N (7450 in US units), Table A.2 of [2].
        A cancels in the ratio X and has no effect on f_fit.
    """
    s = -1.0 if inner else 1.0
    return (A * (2.0 * f / (2.0 * f - 1.0)) ** 0.41
            * (1.0 + s * gamma) ** 1.39
            / (1.0 - s * gamma) ** (1.0 / 3.0)
            * gamma ** 0.30 * Dw ** 1.8 * Z ** (-1.0 / 3.0))


def race_life_ratio_X(Qc_i: float, Qc_e: float, Qei: float, Qee: float) -> float:
    """Outer-to-inner race life ratio X = L_OR / L_IR — equations (3.30), (3.33)."""
    return (Qc_e / Qee) ** 3 / (Qc_i / Qei) ** 3


def f_fit_factor(LR: float, X: float, e: float = 1.11) -> float:
    """
    Fit life factor of the whole bearing — equations (3.34), (3.35).

    From strict-series reliability with Weibull distributions, Zaretsky's rule and
    the assumption that the ball-set life equals the outer-race life:

        f_fit = [ (1 + 2·X^(−e)) / (L_R^(−e) + 2·X^(−e)) ]^(1/e)

    The hoop stress acts on the inner ring only, so f_fit of the bearing is always
    larger than L_R.
    """
    xe = X ** (-e)
    return ((1.0 + 2.0 * xe) / (LR ** (-e) + 2.0 * xe)) ** (1.0 / e)


def L10_iso281(C: float, P: float) -> float:
    """Basic rating life per ISO 281, million revolutions. For pure radial load P = Fr."""
    return (C / P) ** 3


def million_rev_to_hours(L: float, n_rpm: float) -> float:
    """Convert life from million revolutions to operating hours."""
    return 1.0e6 * L / (60.0 * n_rpm)


def hours_to_million_rev(Lh: float, n_rpm: float) -> float:
    """Convert life from operating hours to million revolutions."""
    return 60.0 * n_rpm * Lh / 1.0e6
