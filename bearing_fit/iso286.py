# -*- coding: utf-8 -*-
"""
iso286.py — Limit deviations of the shaft (ISO 286-2) and of the bearing bore (ISO 492).

All returned values are in micrometres (µm).
Step 1 of Fig. 3.1, equations (3.3) – (3.5) of the thesis.
"""

# Nominal size ranges, mm, as (lower bound exclusive, upper bound inclusive)
_RANGES = [(0, 3), (3, 6), (6, 10), (10, 18), (18, 30),
           (30, 50), (50, 80), (80, 120), (120, 180)]

# Standard tolerance grades IT4 – IT7, µm
_IT = {
    4: [3, 4, 4, 5, 6, 7, 8, 10, 12],
    5: [4, 5, 6, 8, 9, 11, 13, 15, 18],
    6: [6, 8, 9, 11, 13, 16, 19, 22, 25],
    7: [10, 12, 15, 18, 21, 25, 30, 35, 40],
}

# Fundamental deviations of the shaft, µm
#   g, h   : upper deviation es (≤ 0)
#   j      : lower deviation ei (separate table in ISO 286-2)
#   k … r  : lower deviation ei (≥ 0)
_FD_ES = {   # es — zones below the zero line
    'g': [-2, -4, -5, -6, -7, -9, -10, -12, -14],
    'h': [0, 0, 0, 0, 0, 0, 0, 0, 0],
}
_FD_EI = {   # ei — zones above the zero line
    'j': [-2, -2, -2, -3, -4, -5, -7, -9, -11],   # valid for grades 5 and 6 only
    'k': [0, 1, 1, 1, 2, 2, 2, 3, 3],
    'm': [2, 4, 6, 7, 8, 9, 11, 13, 15],
    'n': [4, 8, 10, 12, 15, 17, 20, 23, 27],
    'p': [6, 12, 15, 18, 22, 26, 32, 37, 43],
    'r': [10, 15, 19, 23, 28, 34, 41, 51, 63],
}
# Fundamental deviation r above 50 mm uses sub-ranges (ISO 286-2:2010, Table 27), µm
_R_SUB = [(50, 65, 41), (65, 80, 43), (80, 100, 51), (100, 120, 54),
          (120, 140, 63), (140, 160, 65), (160, 180, 68)]
_J_GRADES = (5, 6)   # grades covered by the j row above

# Bore deviation of the inner ring, Normal class (P0), ISO 492, µm.
# The bore tolerance lies below the zero line: ES = 0, EI < 0.
_BORE_TOL = [-8, -8, -8, -8, -10, -12, -15, -20, -25]


def _idx(d_mm: float) -> int:
    """Index of the nominal size range that contains d."""
    if d_mm <= 0:
        raise ValueError(f"Diameter must be positive, got {d_mm} mm.")
    for i, (lo, hi) in enumerate(_RANGES):
        if lo < d_mm <= hi:
            return i
    raise ValueError(f"Diameter {d_mm} mm is outside the table range (0 – 180 mm).")


def shaft_deviations(d_mm: float, fit: str):
    """
    Limit deviations of the shaft according to ISO 286-2.

    fit : tolerance class, e.g. 'k6', 'm5', 'js6', 'g6'.
    Returns (es, ei) in µm.
    """
    if not isinstance(fit, str) or len(fit.strip()) < 2 or not fit.strip()[-1].isdigit():
        raise ValueError(f"Invalid tolerance class: {fit!r} (valid examples: 'k5', 'm6').")
    fit = fit.strip()
    letter = fit[:-1]
    grade = int(fit[-1])
    if grade not in _IT:
        raise ValueError(f"No IT{grade} table; only IT4 – IT7 are available.")
    if letter == 'j' and grade not in _J_GRADES:
        # The j deviations of ISO 286-2 depend on the grade; the table above is
        # valid for j5 and j6 only. Raise instead of returning a wrong value.
        raise ValueError(f"The j row is available for grades 5 and 6 only, not for '{fit}'.")
    i = _idx(d_mm)
    it = _IT[grade][i]

    if letter == 'js':                       # symmetric about the zero line
        return it / 2.0, -it / 2.0
    if letter in _FD_ES:                     # es given
        es = _FD_ES[letter][i]
        return float(es), float(es - it)
    if letter in _FD_EI:                     # ei given
        ei = _FD_EI[letter][i]
        if letter == 'r' and d_mm > 50:      # ISO 286-2 Table 27 splits the ranges above 50 mm
            ei = next(v for lo, hi, v in _R_SUB if lo < d_mm <= hi)
        return float(ei + it), float(ei)
    raise ValueError(f"Fundamental deviation '{letter}' is not supported.")


def bore_deviations(d_mm: float):
    """Limit deviations of the inner-ring bore, Normal class. Returns (ES, EI) in µm."""
    return 0.0, float(_BORE_TOL[_idx(d_mm)])


def interference(d_mm: float, fit: str):
    """
    Diametral interference range between shaft and inner-ring bore — equations (3.3) – (3.5).

    Sign convention: positive = interference, negative = clearance.
    Returns a dict with dmax, dmin, dmean and the four limit deviations used.
    """
    es, ei = shaft_deviations(d_mm, fit)
    ES, EI = bore_deviations(d_mm)
    dmax = es - EI      # largest shaft in smallest bore
    dmin = ei - ES      # smallest shaft in largest bore
    return {'dmin': dmin, 'dmax': dmax, 'dmean': 0.5 * (dmax + dmin),
            'es': es, 'ei': ei, 'ES': ES, 'EI': EI}


# Fit classes used when scanning several classes
FIT_CLASSES = ['g6', 'h6', 'j6', 'k5', 'k6', 'm5', 'm6', 'n6', 'p6', 'r6']
