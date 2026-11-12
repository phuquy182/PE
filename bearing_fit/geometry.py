# -*- coding: utf-8 -*-
"""
geometry.py — Geometry of a single-row deep groove ball bearing.

Defines the Bearing data class and its derived quantities, most importantly the
effective outer diameter of the inner ring D_eff — equations (3.7) and (3.8),
step 3 of Fig. 3.1.
"""
import math
from dataclasses import dataclass

# Angle subtended by the raceway groove in the circular-segment area formula, rad
# (1.8543 rad ≈ 106.24°, shoulder height equal to 20 % of the ball diameter)
THETA_SECTOR = 1.8543


@dataclass(frozen=True)
class Bearing:
    """Geometry and load ratings of a single-row ball bearing."""
    name: str    # bearing designation
    d: float     # bore diameter, mm
    D: float     # outside diameter, mm
    B: float     # width (also inner-ring width), mm
    Dw: float    # ball diameter, mm
    Z: int       # number of balls
    fi: float    # inner raceway conformity
    fe: float    # outer raceway conformity
    Pd: float    # diametral operating clearance, mm
    C: float     # basic dynamic load rating, N
    C0: float    # basic static load rating, N

    def __post_init__(self):
        """Reject parameter sets that have no geometric meaning."""
        if not (self.d > 0 and self.D > self.d and self.B > 0):
            raise ValueError(f"Bearing {self.name}: need 0 < d < D and B > 0.")
        if not (0 < self.Dw < 0.5 * (self.D - self.d)):
            raise ValueError(f"Bearing {self.name}: ball diameter Dw must lie in (0, (D - d)/2).")
        if not (isinstance(self.Z, int) and self.Z >= 3):
            raise ValueError(f"Bearing {self.name}: number of balls Z must be an integer >= 3.")
        if self.Z * self.Dw >= math.pi * self.dm:
            raise ValueError(f"Bearing {self.name}: Z*Dw exceeds the pitch circumference.")
        if not (0.5 < self.fi < 1.0 and 0.5 < self.fe < 1.0):
            raise ValueError(f"Bearing {self.name}: conformities fi, fe must lie in (0.5, 1).")
        if self.Pd < 0 or self.C <= 0 or self.C0 <= 0:
            raise ValueError(f"Bearing {self.name}: need Pd >= 0, C > 0 and C0 > 0.")

    # ----- derived quantities -----
    @property
    def dm(self) -> float:
        """Pitch diameter."""
        return 0.5 * (self.d + self.D)

    @property
    def gamma(self) -> float:
        """Auxiliary parameter γ = Dw / dm — used in (3.31), (3.32)."""
        return self.Dw / self.dm

    @property
    def D_IR(self) -> float:
        """Inner raceway (groove bottom) diameter."""
        return self.dm - self.Dw

    @property
    def A_sector(self) -> float:
        """Area of the groove removed from the axial section, mm² — equation (3.7)."""
        r = 0.5 * self.Dw
        return 0.5 * r * r * (THETA_SECTOR - math.sin(THETA_SECTOR))

    @property
    def D_eff(self) -> float:
        """Effective outer diameter of the inner ring — equation (3.8)."""
        return self.D_IR + 0.2 * self.Dw - self.A_sector / self.B


# ---------------------------------------------------------------------------
# Bearings used in the thesis.
# Load ratings from the SKF catalogue; internal geometry per ISO/TS 16281.
# To change one value use dataclasses.replace(B6310, Pd=0.0); do not mutate.
# ---------------------------------------------------------------------------
B6210 = Bearing(name='6210', d=50.0, D=90.0, B=20.0, Dw=12.700, Z=10,
                fi=0.52, fe=0.52, Pd=0.013, C=37100.0, C0=23200.0)

B6310 = Bearing(name='6310', d=50.0, D=110.0, B=27.0, Dw=19.050, Z=8,
                fi=0.52, fe=0.52, Pd=0.013, C=65000.0, C0=38000.0)

# Paper 2: SKF 6306. C, C0 from the SKF catalogue [16]. Dw and Z are not published;
# the assumed values reproduce C within 0.1 % by ISO 281:2007 Eq. (1).
# P_d,op = 12.5 µm: middle of the Normal clearance group, 5 - 20 µm for d = 24 - 30 mm (SKF [16]).
B6306 = Bearing(name='6306', d=30.0, D=72.0, B=19.0, Dw=12.700, Z=8,
                fi=0.52, fe=0.52, Pd=0.0125, C=29600.0, C0=16000.0)

# Size-210 bearing of the worked example in [9]; used by run_verification.py (Table 3.4)
B210_REF = Bearing(name='210 (reference)', d=50.0, D=90.0, B=20.0, Dw=12.700, Z=10,
                   fi=0.52, fe=0.52, Pd=0.013, C=37100.0, C0=23200.0)
