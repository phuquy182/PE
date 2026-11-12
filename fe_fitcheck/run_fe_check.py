# -*- coding: utf-8 -*-
"""
run_fe_check.py — Runs all FE cases of task 3 (Paper 2, direction 1) and writes
results/fe_results.json and results/fe_results.csv.

Cases (handover note 13/11, section 6):
  A  plain ring, outer diameter = D_eff of [9], shaft as long as the ring  (code check)
  B  real groove profile, shaft as long as the ring
  C  real groove profile, shaft longer than the ring by B on each side       (conclusion case)
  D  as C, shoulder depth 0.15 and 0.25 Dw                                     (geometry probe)
Plus: mesh convergence (6310, case C, three levels), linearity (6310, C, 30 µm),
shaft-length sensitivity (6310, C, extension 2B).

The original bearing_fit package is only imported, never modified.
"""
import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))          # bearing_fit
sys.path.insert(0, str(HERE))

import numpy as np

from bearing_fit import B6306, B6210, B6310, closed_form_constants, Params
from bearing_fit import fitmodel
import fe_axisym as F

OUT = HERE / 'results'
OUT.mkdir(exist_ok=True)

PRM = Params()
E, NU = PRM.E, PRM.nu
BEARINGS = [B6306, B6210, B6310]
DELTA = 0.010            # effective diametral interference, mm (10 µm)
LAMBDA = 0.08            # Fr / C for D_z
LEVEL_FINAL = 2          # mesh level used for all production runs


def deff_equal_area(brg):
    """Outer diameter of the plain ring with the same axial-section area, using the
    groove area of the program (0.1118 Dw², shoulder 0.2 Dw)."""
    return brg.D_IR + 0.4 * brg.Dw - 2.0 * brg.A_sector / brg.B


def mesh_sizes(brg, z0, level):
    s = 2 ** level
    return dict(h_glob=brg.B / 16 / s, h_fine=z0 / 2 / s, h_edge=brg.B / 64 / s), z0 / 2 / s


def formula(brg, D_eff, Dz, delta=DELTA):
    p = fitmodel.fit_pressure(delta, brg.d, D_eff, E)
    return dict(D_eff=D_eff, p_f=p,
                s_bore=fitmodel.hoop_stress_at_bore(p, brg.d, D_eff),
                s_z0=fitmodel.hoop_stress(p, brg.d, D_eff, Dz),
                k_cl=brg.d / D_eff,
                k_cl_lame_DIR=2.0 * float(F.lame_ur(p, brg.d, D_eff, 0.5 * brg.D_IR, E, NU)) / delta)


def run(brg, case, level=LEVEL_FINAL, delta=DELTA, shoulder=0.2, ext=1.0):
    k = closed_form_constants(brg, LAMBDA * brg.C, PRM)
    Dz = k['Dz']
    z0 = 0.5 * (brg.D_IR - Dz)
    sizes, ds_arc = mesh_sizes(brg, z0, level)
    hB = 0.5 * brg.B
    if case == 'A':
        sec = F.plain_ring_section(brg.d, brg.B, brg.D_eff, hB)
    else:
        Ls = hB if case == 'B' else hB + ext * brg.B
        sec = F.grooved_ring_section(brg.d, brg.B, brg.D_IR, brg.Dw, fi=brg.fi,
                                     shoulder=shoulder, L_shaft=Ls, ds_arc=ds_arc)
    t0 = time.time()
    r = F.solve(sec, delta, E=E, nu=NU, **sizes)
    dt = time.time() - t0
    s_z0 = r.sig_theta_at(0.5 * Dz)
    pz = r.p_nodes
    res = dict(bearing=brg.name, case=case, level=level, delta_um=delta * 1e3,
               shoulder=shoulder, shaft_ext_B=(0.0 if case in 'AB' else ext),
               n_dof=r.n_dof, t_s=round(dt, 2), z0=z0, Dz=Dz,
               h_fine=sizes['h_fine'], h_glob=sizes['h_glob'],
               p_f=r.p_mean, p_node_min=float(pz[:, 1].min()), p_node_max=float(pz[:, 1].max()),
               all_compressive=r.all_compressive,
               s_bore_mid=r.sig_bore_mid, s_bore_avg=r.sig_bore_avg, s_z0=s_z0,
               u_r_groove=r.u_r_groove, k_cl=2.0 * r.u_r_groove / delta,
               ring_area_half=F.ring_area(sec))
    # equivalent plain-ring diameter that reproduces the FE mean pressure by (3.9)
    x = 1.0 - 2.0 * r.p_mean * brg.d / (E * delta)
    res['D_eff_backcalc_p'] = brg.d / math.sqrt(x) if x > 0 else float('nan')
    res['D_eq_area_exact'] = 2.0 * (0.5 * brg.d + res['ring_area_half'] / hB)
    return res, r


def main():
    rows, extra = [], {}
    # ---- mesh convergence, 6310 case C
    conv = []
    for lev in range(3):
        res, _ = run(B6310, 'C', level=lev)
        conv.append(res)
        print('conv', lev, res['n_dof'], '%.5f %.5f %.5f %.6f' % (res['p_f'], res['s_bore_mid'],
                                                                  res['s_z0'], res['k_cl']))
    extra['convergence'] = conv
    # ---- main cases
    fe_fields = {}
    for brg in BEARINGS:
        for case in ('A', 'B', 'C'):
            res, r = run(brg, case)
            rows.append(res)
            if brg.name == '6310' and case == 'C':
                fe_fields['6310C'] = r
            print(brg.name, case, res['n_dof'], '%.4f %.4f %.4f %.5f %s' % (
                res['p_f'], res['s_bore_mid'], res['s_z0'], res['k_cl'], res['all_compressive']))
        for sh in (0.15, 0.25):
            res, _ = run(brg, 'D', shoulder=sh)
            rows.append(res)
            print(brg.name, 'D', sh, '%.4f %.4f %.4f %.5f' % (res['p_f'], res['s_bore_mid'],
                                                              res['s_z0'], res['k_cl']))
    # ---- linearity and shaft length, 6310 case C
    res30, _ = run(B6310, 'C', delta=0.030)
    extra['linearity_30um'] = res30
    res_ext, _ = run(B6310, 'C', ext=2.0)
    extra['shaft_ext_2B'] = res_ext
    # ---- formula values
    forms = {}
    for brg in BEARINGS:
        Dz = closed_form_constants(brg, LAMBDA * brg.C, PRM)['Dz']
        forms[brg.name] = dict(ref9=formula(brg, brg.D_eff, Dz),
                               equal_area=formula(brg, deff_equal_area(brg), Dz),
                               D_IR=brg.D_IR, Dw=brg.Dw, d=brg.d, B=brg.B, A_sector=brg.A_sector)
    out = dict(constants=dict(E=E, nu=NU, delta_mm=DELTA, Fr_over_C=LAMBDA),
               rows=rows, formula=forms, **extra)
    (OUT / 'fe_results.json').write_text(json.dumps(out, indent=1, default=float))
    keys = list(rows[0].keys())
    with open(OUT / 'fe_results.csv', 'w', encoding='utf-8') as f:
        f.write(','.join(keys) + '\n')
        for r_ in rows + conv + [res30, res_ext]:
            f.write(','.join(str(r_[k_]) for k_ in keys) + '\n')

    # ---- fields for the figures (6310, case C)
    r = fe_fields['6310C']
    rg = r.ring
    m = rg['mesh']
    np.savez(OUT / 'field_6310C.npz', p=m.p, t=m.t,
             sig_nodes=rg['sig_t'][rg['basis_s'].nodal_dofs[0]],
             p_nodes=r.p_nodes)
    b = B6310
    Dz = closed_form_constants(b, LAMBDA * b.C, PRM)['Dz']
    rr = np.linspace(0.5 * b.d, 0.5 * b.D_IR, 300)
    np.savez(OUT / 'radial_6310C.npz', r=rr, sig=r.sig_theta_at(rr, 0.0), Dz=Dz)
    print('done')


if __name__ == '__main__':
    main()
