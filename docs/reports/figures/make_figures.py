"""
make_figures.py -- figures of docs/reports/dyer_supercharge_correction.tex

Every number below is copied from a result file of the repository:
  validation/waxman_2013_results.md                                  (Table 1, baseline bands)
  validation/kappa_correction/loocv_kappa_results.txt                (per fold, bands, limit check)
  validation/kappa_correction/loocv_functional_forms_results.txt     (forms, fitted parameters)
  validation/kappa_correction/scan_threshold_continuity_results.txt  (tank-pressure scan)
Run (needs only numpy + matplotlib):   python make_figures.py
Sizes are for a two-column paper: 3.4 in = one column, 7.0 in = full width.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.dirname(os.path.abspath(__file__))
BASE, FREE, PROT, ACC, GRID = "#6B7280", "#1F4E9C", "#D97706", "#2E7D6B", "#D1D5DB"
plt.rcParams.update({
    "font.family": "serif", "font.serif": ["STIXGeneral", "DejaVu Serif"], "mathtext.fontset": "stix",
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": .5, "grid.alpha": .8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.axisbelow": True,
    "legend.frameon": False, "legend.fontsize": 7, "savefig.bbox": "tight", "savefig.pad_inches": .03})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name)); plt.close(fig)


def labels(ax, bars, fmt="{:.1f}", dy=.25, size=6.5):
    for b in bars:
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + dy, fmt.format(b.get_height()),
                ha="center", va="bottom", fontsize=size)


# 1. baseline error by supercharge band (Table 1) -- one column ---------------------------------------------
fig, ax = plt.subplots(figsize=(3.3, 2.1))
b = ax.bar(["< 100 psi\n(n = 19)", "100-200 psi\n(n = 13)", ">= 200 psi\n(n = 32)"], [15.51, 6.24, 1.97],
           color=[PROT, "#E9A94B", BASE], width=.55)
labels(ax, b, "{:.2f} %", .3, 7)
ax.set_ylabel("baseline Dyer MAPE [%]"); ax.set_ylim(0, 18.5)
save(fig, "fig_baseline_bands.pdf")

# 2. shapes of the multiplier f(S) = kappa'/kappa -- full width -----------------------------------------------
S = np.linspace(0, 1000, 801)
fig, ax = plt.subplots(figsize=(7.0, 2.6))
ax.axvspan(41, 371, color="#E5E7EB", alpha=.9, lw=0)
ax.text(206, 1.54, "range of the data (41-371 psi)", ha="center", va="top", fontsize=7, color="#374151")
ax.axvline(94.7, color="k", lw=.7, ls=":"); ax.text(100, .03, "Part A (95 psi)", fontsize=6.5, va="bottom")
ax.plot(S, np.ones_like(S), color=BASE, lw=1.4, label="F0  baseline Dyer ($f=1$)")
ax.plot(S, np.where(S < 433, (S / 433) ** .9, 1), color=FREE, lw=1.6, label=r"F1  gated power ($\beta$=0.9, $S_{ref}$=433)")
ax.plot(S, np.minimum(1, S / 411), color=FREE, lw=1.1, ls="--", label=r"F2  linear, capped: min(1, $S/S_0$), $S_0$=411")
ax.plot(S, S / 411, color=FREE, lw=.9, ls=":", label=r"F2  linear, uncapped ($S/S_0$)")
ax.plot(S, S / (S + 236), color=ACC, lw=1.6, ls="-.", label=r"F6  soft gate $S/(S+S_1)$, $S_1$=236")
ax.plot(S, np.where(S < 90, (S / 90) ** 2.5, 1), color=PROT, lw=1.6, label=r"protected gate ($\beta$=2.5, $S_{ref}$=90)")
ax.set_xlim(0, 1000); ax.set_ylim(0, 1.6)
ax.set_xlabel("supercharge $S = P_{up}-P_{sat}(T_{up})$ [psi]"); ax.set_ylabel(r"multiplier  $f=\kappa'/\kappa$")
ax.legend(loc="center left", bbox_to_anchor=(1.01, .5))
save(fig, "fig_factor_shapes.pdf")

# 3. leave-one-curve-out, per fold -- full width --------------------------------------------------------------------
sup = [41, 79, 115, 169, 206, 269, 290, 329, 371]
base = [17.80, 13.67, 6.50, 5.93, 3.58, 1.08, 2.16, 1.60, 1.19]
free = [4.07, 1.31, 1.48, 0.75, 1.00, 0.84, 0.39, 0.55, 0.87]
gate = [3.65, 11.93, 6.50, 5.93, 3.58, 1.08, 2.16, 1.60, 1.19]
x = np.arange(len(sup)); w = .27
fig, ax = plt.subplots(figsize=(7.0, 2.3))
b0 = ax.bar(x - w, base, w, color=BASE, label="baseline Dyer")
b1 = ax.bar(x, free, w, color=FREE, label=r"gated $\kappa$, free fit")
b2 = ax.bar(x + w, gate, w, color=PROT, label=r"gated $\kappa$, $S_{ref}\leq$ 95 psi")
labels(ax, b0, "{:.1f}", .2, 6.5); labels(ax, b1, "{:.1f}", .2, 6.5)
ax.set_xticks(x); ax.set_xticklabels([str(s) for s in sup])
ax.set_xlabel("held-out curve: tank supercharge [psi]"); ax.set_ylabel("out-of-sample MAPE [%]")
ax.legend(); ax.set_ylim(0, 20.5)
save(fig, "fig_loocv_folds.pdf")

# 4. out-of-sample by supercharge band -- one column -----------------------------------------------------------------
vb, vf, vg = [15.62, 6.24, 1.95], [2.62, 1.14, 0.70], [8.01, 6.24, 1.95]
x = np.arange(3); fig, ax = plt.subplots(figsize=(3.3, 2.2))
bars = [ax.bar(x - w, vb, w, color=BASE, label="baseline Dyer"),
        ax.bar(x, vf, w, color=FREE, label=r"gated $\kappa$, free fit"),
        ax.bar(x + w, vg, w, color=PROT, label=r"gated $\kappa$, $S_{ref}\leq$ 95 psi")]
for bs in bars: labels(ax, bs, "{:.2f}", .2, 6)
ax.set_xticks(x); ax.set_xticklabels(["< 100 psi", "100-200 psi", ">= 200 psi"])
ax.set_ylabel("out-of-sample MAPE [%]"); ax.set_ylim(0, 19.5); ax.legend()
save(fig, "fig_loocv_bands.pdf")

# 5. forms comparison -- one column ---------------------------------------------------------------------------------------
names = ["F0 baseline Dyer (0)", "F1 gated power (2)", "F2 linear (1)", "F4 $\\sqrt{\\kappa^2-1}$ (0)",
         "F5 $\\kappa^2-1$ (0)", "F6 soft gate (1)", "F7 soft gate, exponent (2)"]
mape = [6.88, 1.30, 1.17, 3.66, 3.19, 2.25, 2.66]
cols = [BASE, FREE, FREE, ACC, ACC, PROT, PROT]
fig, ax = plt.subplots(figsize=(3.3, 2.3))
y = np.arange(len(names))[::-1]
ax.barh(y, mape, color=cols, height=.62)
for yy, m in zip(y, mape): ax.text(m + .08, yy, f"{m:.2f} %", va="center", fontsize=6.5)
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=7); ax.set_xlim(0, 8.4)
ax.set_xlabel("pooled out-of-sample MAPE [%] (mode A)")
save(fig, "fig_forms.pdf")

# 6. limit check -- one column -------------------------------------------------------------------------------------------------
Sx = [.05, .5, 2, 5, 10, 20, 41, 95, 200]
dy = [42.60/36.78, 42.70/36.86, 43.03/37.13, 43.68/37.66, 44.76/38.53, 46.84/40.20, 50.97/43.51, 60.42/51.04, 75.79/63.17]
c1 = [1.000, 1.000, 1.002, 1.004, 1.008, 1.016, 1.034, 1.078, 1.146]
c2 = [1.000, 1.000, 1.000, 1.000, 1.001, 1.008, 1.045, 1.184, 1.200]
fig, ax = plt.subplots(figsize=(3.3, 2.4))
ax.semilogx(Sx, dy, "o--", ms=3, color=BASE, label="Dyer, uncorrected")
ax.semilogx(Sx, c1, "s-", ms=3, color=FREE, label=r"gated, $\beta$=1, $S_{ref}$=400 psi")
ax.semilogx(Sx, c2, "^:", ms=3.5, color=PROT, label=r"gated, $\beta$=2.5, $S_{ref}$=90 psi")
ax.axhline(1, color="k", lw=.7, ls=":")
ax.text(.06, 1.095, "uncorrected: $(1+\\dot m_{SPI}/\\dot m_{HEM})/2\\approx1.16$", fontsize=6.5, color="#374151")
ax.set_ylim(.99, 1.32)
ax.set_xlabel("supercharge $S$ [psi]"); ax.set_ylabel("prediction / HEM")
ax.legend(loc="upper left")
save(fig, "fig_limit.pdf")

# 7. tank-pressure scan (needle valve, 8 mm) -- one column ---------------------------------------------------------------------
P = np.arange(52.0, 58.0001, .25); nan = np.nan
bs_ = [186.7, 187.9, 189.1, 190.3, 191.5, 192.7, 193.9] + [nan]*6 + [340.3, 341.7, 343.1, 344.4, 345.8, 347.2, 348.5, 349.8, 351.1, 352.5, 353.8, 355.1]
co_ = [186.7, 187.9, 189.1, 190.3, 191.5, 192.7, 193.9, 196.1, 199.1, 202.1, 205.1, 208.0, 211.0, 213.9, 216.9, 219.8, 222.7, 225.6, 228.4, 231.3, 234.1, 237.0, 239.8, 242.6, 245.4]
pr_ = [186.7, 187.9, 189.1, 190.3, 191.5, 192.7, 193.9, 194.8, 195.7, 197.0, 198.8, 201.2, 204.1, 207.6, 211.6, 216.2, 221.2, 226.6, 232.5, 238.6, 245.0, 251.6, 258.4, 265.4, 272.4]
assert len(bs_) == len(co_) == len(pr_) == len(P)
fig, ax = plt.subplots(figsize=(3.3, 2.5))
ax.axvspan(53.6, 55.1, color="#F3D5D5", alpha=.8, lw=0)
ax.text(54.35, 262, "baseline:\nno convergence", ha="center", va="center", fontsize=6.5, color="#7F1D1D")
ax.annotate("", xy=(55.25, 340.3), xytext=(53.5, 193.9), arrowprops=dict(arrowstyle="->", color=BASE, ls="--", lw=.8))
ax.text(54.95, 232, "+75 %", fontsize=7, color="#374151")
ax.plot(P, bs_, "o--", ms=2.8, color=BASE, label="baseline")
ax.plot(P, co_, "s-", ms=2.5, color=FREE, label=r"gated $\kappa$ ($\beta$=1, 400 psi)")
ax.plot(P, pr_, "^:", ms=3, color=PROT, label=r"gated $\kappa$ ($\beta$=2.5, 90 psi)")
ax.set_xlabel("tank pressure [bar]"); ax.set_ylabel("converged mass flow [g/s]")
ax.legend(loc="upper left"); ax.set_ylim(170, 375)
save(fig, "fig_scan.pdf")
print("figures written to", OUT)
