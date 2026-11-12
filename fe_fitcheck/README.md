# fe_fitcheck — FE check of the fit branch (Paper 2, direction 1, task 3)

Axisymmetric linear-elastic FE model (scikit-fem + triangle) of a solid shaft and a
deep-groove-ball-bearing inner ring with interference, compared with the Lamé
formulas of `bearing_fit.fitmodel` using two effective diameters D_eff.

```bash
pip install scikit-fem==12.0.2 triangle scipy matplotlib reportlab
```

```bash
python run_fe_check.py
```

```bash
python make_figures.py && python make_report.py && python make_note.py
```

`run_fe_check.py` takes about 2 min. Outputs go to `results/`. The package
`../bearing_fit` is only imported, never modified (rebuilt from the Drive copy of
27/10/2026; regression: f_fit = 0.8591, delta_max_life = 15.70 µm).
