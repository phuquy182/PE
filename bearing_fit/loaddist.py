# -*- coding: utf-8 -*-
"""
loaddist.py — Internal load distribution among the balls under pure radial load.

Follows ISO/TS 16281, reduced to Fa = 0 with no tilting moment.
Step 7 of Fig. 3.1, equations (3.12) – (3.14).
"""
import math


def total_stiffness(Ki: float, Ke: float) -> float:
    """
    Combined stiffness of the two contacts in series (ball – inner and ball – outer).
    Deflections add, so the inverse stiffnesses add with the exponent 2/3.
    """
    return (Ki ** (-2.0 / 3.0) + Ke ** (-2.0 / 3.0)) ** (-1.5)


def ball_loads(delta_r: float, Z: int, Pd: float, Kn: float):
    """
    Load on each ball for a radial displacement delta_r of the inner ring.
    Equations (3.12), (3.14). Balls outside the load zone carry 0.
    """
    loads = []
    for j in range(Z):
        psi = 2.0 * math.pi * j / Z                   # angular position of ball j
        dj = delta_r * math.cos(psi) - 0.5 * Pd       # elastic deflection of ball j
        loads.append(Kn * dj ** 1.5 if dj > 0.0 else 0.0)
    return loads


def radial_reaction(delta_r: float, Z: int, Pd: float, Kn: float) -> float:
    """Sum of ball-load projections on the load direction — right-hand side of (3.13)."""
    s = 0.0
    for j in range(Z):
        psi = 2.0 * math.pi * j / Z
        dj = delta_r * math.cos(psi) - 0.5 * Pd
        if dj > 0.0:
            s += Kn * dj ** 1.5 * math.cos(psi)
    return s


def solve_delta_r(Fr: float, Z: int, Pd: float, Kn: float,
                  tol: float = 1e-12, itmax: int = 200) -> float:
    """
    Solve the static equilibrium (3.13) by bisection for the radial
    displacement delta_r of the inner ring under the load Fr.
    """
    lo = max(0.5 * Pd, 0.0)          # below this no ball is loaded
    hi = lo + 1.0e-4
    while radial_reaction(hi, Z, Pd, Kn) < Fr:   # expand the bracket
        hi *= 2.0
        if hi > 100.0:
            raise RuntimeError("No bracket found for delta_r.")
    for _ in range(itmax):
        mid = 0.5 * (lo + hi)
        if radial_reaction(mid, Z, Pd, Kn) < Fr:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def equivalent_ball_loads(loads, Z: int):
    """
    Equivalent ball loads of Lundberg – Palmgren theory.
    Rotating ring (inner): cubic mean; stationary ring (outer): mean of power 10/3.
    Used for the race life ratio X in module life.
    """
    Qei = (sum(q ** 3 for q in loads) / Z) ** (1.0 / 3.0)
    Qee = (sum(q ** (10.0 / 3.0) for q in loads) / Z) ** (3.0 / 10.0)
    return Qei, Qee
