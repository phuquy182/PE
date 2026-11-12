# -*- coding: utf-8 -*-
"""
bearing_fit — Effect of an interference fit on the fatigue life of a deep groove ball
bearing (the thirteen steps of Fig. 3.1 of the thesis).

Example:
    from bearing_fit import evaluate, B6310, delta_max_life
    r = evaluate(B6310, 'k5', 5200.0)      # at the mean interference of k5
    print(r.f_fit)                          # 0.8591
    print(delta_max_life(B6310, 5200.0, 18000.0, 1500.0))   # 15.70 µm
"""
import sys as _sys

# The scripts print Vietnamese text. On Windows, redirecting output to a file
# (python run_chapter4.py > log.txt) uses cp1252 by default and fails;
# switch the standard streams to UTF-8 to avoid this.
for _s in (_sys.stdout, _sys.stderr):
    if hasattr(_s, 'reconfigure') and (_s.encoding or '').lower().replace('-', '') != 'utf8':
        _s.reconfigure(encoding='utf-8')

from .geometry import Bearing, B6210, B6306, B6310, B210_REF
from .model import (Params, Result, evaluate, life_hours, delta_max_life,
                    closed_form_constants, ntn_min_interference)
from .iso286 import interference

__all__ = ['Bearing', 'B6210', 'B6306', 'B6310', 'B210_REF', 'Params', 'Result', 'evaluate',
           'life_hours', 'delta_max_life', 'closed_form_constants',
           'ntn_min_interference', 'interference']
