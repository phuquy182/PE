# -*- coding: utf-8 -*-
"""
make_figures.py — Figures of the FE check (6310, case C) from results/*.npz.

Style follows bearing_fit/figstyle.py of the thesis (copied constants, the original
module is not imported so this folder stays independent).
"""
import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.ticker import FuncFormatter

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
from bearing_fit import B6310, fitmodel           # noqa: E402
import fe_axisym as F                              # noqa: E402

RES = HERE / 'results'
INK, INK2 = '#1a1a1a', '#555555'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#6250d6']
E, NU = 205878.0, 0.30


def apply(font_size=10):
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': font_size, 'axes.titlesize': font_size,
        'axes.labelsize': font_size, 'legend.fontsize': font_size - 1.5,
        'xtick.labelsize': font_size - 0.5, 'ytick.labelsize': font_size - 0.5, 'figure.dpi': 300,
        'axes.grid': True, 'grid.color': '#7a7a7a', 'grid.linewidth': 0.5,
        'grid.linestyle': (0, (2, 2)), 'grid.alpha': 0.45,
        'axes.edgecolor': '#4d4d4d', 'axes.linewidth': 0.8, 'axes.labelcolor': INK,
        'xtick.color': INK, 'ytick.color': INK, 'text.color': INK, 'legend.frameon': False,
    })


def vn_fmt(nd=None):
    def f(x, _):
        s = f'{x:.{nd}f}' if nd is not None else f'{x:g}'
        return s.replace('.', ',').replace('-', '−')
    return FuncFormatter(f)


def label_box():
    return dict(boxstyle='square,pad=0.15', fc='white', ec='none', alpha=0.85)


def fig_field():
    z = np.load(RES / 'field_6310C.npz')
    p, t, s = z['p'], z['t'], z['sig_nodes']
    b = B6310
    Dz = float(np.load(RES / 'radial_6310C.npz')['Dz'])
    tri = mtri.Triangulation(p[0], p[1], t.T)
    fig, (ax, axz) = plt.subplots(1, 2, figsize=(6.3, 3.3), gridspec_kw=dict(width_ratios=[1.55, 1]))
    levels = np.arange(np.floor(s.min()), np.ceil(s.max()) + 0.5, 1.0)
    for a_ in (ax, axz):
        cf = a_.tricontourf(tri, s, levels=levels, cmap='Blues')
        a_.set_aspect('equal')
        a_.grid(False)
    # shaft outline (part near the ring) for orientation
    a, hB = 0.5 * b.d, 0.5 * b.B
    ax.axvspan(a - 3.5, a, color='#e6e6e6', lw=0)
    ax.text(a - 1.75, 0.5 * hB, 'trục', ha='center', va='center', color=INK2, fontsize=8.5,
            rotation=90)
    pts = [(a, 0.0, '1  lỗ, $r=d/2$'), (0.5 * Dz, 0.0, '2  độ sâu $z_0$, $r=D_z/2$'),
           (0.5 * b.D_IR, 0.0, '3  đáy rãnh, $r=D_{IR}/2$')]
    mk = ['o', 's', '^']
    for (r_, z_, lab), m_ in zip(pts, mk):
        ax.plot(r_, z_, m_, ms=6, mfc='white', mec=INK, mew=1.0, clip_on=False, zorder=5)
        if r_ > 0.5 * b.D_IR - 1.2:
            axz.plot(r_, z_, m_, ms=6, mfc='white', mec=INK, mew=1.0, clip_on=False, zorder=5)
    ax.text(a + 0.5, 0.9, '1', fontsize=8.5, bbox=label_box())
    ax.text(0.5 * b.D_IR - 1.9, 0.9, '2, 3', fontsize=8.5, bbox=label_box())
    ax.set_xlim(a - 3.5, 0.5 * b.D_IR + 0.2 * b.Dw + 1.0)
    ax.set_ylim(0, hB + 0.5)
    ax.set_xlabel('r, mm')
    ax.set_ylabel('z, mm')
    # zoom on the groove bottom
    axz.set_xlim(0.5 * b.D_IR - 1.2, 0.5 * b.D_IR + 0.4)
    axz.set_ylim(0, 1.0)
    axz.set_xlabel('r, mm')
    for a_ in (ax, axz):
        a_.xaxis.set_major_formatter(vn_fmt())
        a_.yaxis.set_major_formatter(vn_fmt())
    axz.set_title('phóng to vùng đáy rãnh', fontsize=8.5, color=INK2)
    cb = fig.colorbar(cf, ax=[ax, axz], orientation='horizontal', fraction=0.06, pad=0.2, aspect=40)
    cb.set_label('σθ, MPa (δ_eff = 10 µm)')
    cb.set_ticks(levels[::2])
    cb.ax.xaxis.set_major_formatter(vn_fmt())
    handles = [plt.Line2D([], [], ls='none', marker=m_, mfc='white', mec=INK, label=lab)
               for (_, _, lab), m_ in zip(pts, mk)]
    fig.legend(handles=handles, loc='upper center', ncol=3, fontsize=7.8, bbox_to_anchor=(0.5, 1.0))
    fig.savefig(RES / 'fig_field_6310C.png', bbox_inches='tight', dpi=300)
    plt.close(fig)


def fig_radial():
    z = np.load(RES / 'radial_6310C.npz')
    r, s, Dz = z['r'], z['sig'], float(z['Dz'])
    d = json.loads((RES / 'fe_results.json').read_text())
    f = d['formula']['6310']
    b = B6310
    delta = d['constants']['delta_mm']
    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    ax.plot(r, s, color=SERIES[0], lw=2.0, ls='-', label='FE, ca C (z = 0)')
    for key, col, ls, lab in (('ref9', SERIES[1], (0, (5, 2)), 'Lamé, $D_{eff}$ theo [9] = %s mm'),
                              ('equal_area', SERIES[2], (0, (6, 2, 1.2, 2)),
                               'Lamé, $D_{eff}$ cùng diện tích = %s mm')):
        De = f[key]['D_eff']
        p = fitmodel.fit_pressure(delta, b.d, De, E)
        ax.plot(r, F.lame_hoop(p, b.d, De, r), color=col, lw=1.6, ls=ls,
                label=lab % f'{De:.2f}'.replace('.', ','))
    for x, lab, dx, ha in ((0.5 * Dz, '$D_z/2$', -0.05, 'right'), (0.5 * b.D_IR, '$D_{IR}/2$', 0.05, 'left')):
        ax.axvline(x, color=INK2, lw=0.7, ls=(0, (1, 1.5)))
        ax.text(x + dx, 34.6, lab, ha=ha, va='top', fontsize=8.5, color=INK2, bbox=label_box())
    ax.set_xlim(24.7, 31.3)
    ax.set_xlabel('r, mm')
    ax.set_ylabel('σθ, MPa')
    ax.set_ylim(24, 35)
    ax.xaxis.set_major_formatter(vn_fmt())
    ax.yaxis.set_major_formatter(vn_fmt())
    ax.legend(loc='lower left')
    fig.savefig(RES / 'fig_radial_6310C.png', bbox_inches='tight', dpi=300)
    plt.close(fig)


def fig_pressure():
    """Nodal contact pressure along the fit surface, 6310, cases B and C (normalised)."""
    import fe_axisym as F2
    from run_fe_check import run
    fig, ax = plt.subplots(figsize=(6.0, 2.9))
    for case, col, ls in (('B', SERIES[3], (0, (5, 2))), ('C', SERIES[0], '-')):
        res, r = run(B6310, case)
        pz = r.p_nodes
        o = np.argsort(pz[:, 0])
        ax.plot(pz[o, 0] / (0.5 * B6310.B), pz[o, 1] / res['p_f'], color=col, ls=ls, lw=1.6,
                label=f'ca {case}: p trung bình = {res["p_f"]:.2f} MPa'.replace('.', ','))
    ax.set_xlabel('z / (B/2)')
    ax.set_ylabel('p / p trung bình')
    ax.set_ylim(0, 3)
    ax.xaxis.set_major_formatter(vn_fmt())
    ax.yaxis.set_major_formatter(vn_fmt())
    ax.legend(loc='upper left')
    fig.savefig(RES / 'fig_pressure_6310.png', bbox_inches='tight', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    apply()
    fig_field()
    fig_radial()
    fig_pressure()
    print('figures written')
