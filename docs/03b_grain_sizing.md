# 3b. Fuel Grain Sizing: From Oxidiser Mass Flow to Grain Geometry

## 3b.1 Why a hybrid motor's fuel side is a different problem from a solid motor's

In a solid rocket motor, fuel and oxidiser are pre-mixed throughout the propellant grain; the burn rate is governed by chamber pressure, via Saint-Robert's law, $\dot r = a P^n$. In a hybrid motor, the grain is **fuel only** — paraffin, HTPB, or another solid polymer/wax — with the oxidiser (here, N₂O) flowing through a central channel (the *port*) cast or machined into the grain.

The flame does not sit at the fuel surface itself; it sits inside a turbulent diffusion boundary layer a short distance above it, where fuel vapour (driven off the hot surface) meets oxidiser diffusing in from the core flow. It is **convective heat transfer from this boundary layer flame back to the fuel surface** that drives regression (surface recession) — not chamber pressure directly. This is the physical reason the fuel side of a hybrid motor needs a fundamentally different model from a solid motor's Saint-Robert's law.

## 3b.2 The Marxman regression rate correlation

Marxman's classical turbulent-boundary-layer analysis, and decades of subsequent empirical work, converge on a power-law correlation in the **oxidiser mass flux** through the port (not the absolute oxidiser mass flow):

$$\dot r = a\, G_o^n, \qquad G_o = \frac{\dot m_{ox}}{A_{port}}$$

- $\dot r$ — regression rate (radial recession speed of the fuel surface)
- $G_o$ — oxidiser mass flux, kg/(m²·s) — mass flow **per unit port cross-sectional area**, not the total flow
- $a$, $n$ — empirical coefficients, fit from test data for a specific fuel/oxidiser/injector/scale combination

Marxman's own boundary-layer derivation predicts $n \approx 0.8$. Real fitted values for most propellant combinations cluster closer to $n \approx 0.5$; the gap reflects physics the simple analysis does not fully capture — blowing (fuel vapour injection into the boundary layer reduces the heat transfer coefficient) and, for paraffin specifically, an additional mechanism where a thin melted liquid layer on the surface is mechanically entrained into the flow by shear, adding a second regression pathway beyond pure vaporisation and pushing paraffin's regression rate well above what a purely vapour-phase model like HTPB's would predict at the same $G_o$.

**Why this project treats $a$ and $n$ as required inputs, not defaults.** Unlike a saturation pressure correlation or a material density, $a$ and $n$ are not properties of "paraffin" or "HTPB" in the abstract — they are fit parameters tied to a specific test article: injector design (which sets the turbulence and mixing at the port inlet), motor scale, and often chamber pressure range. Published values for nominally the same paraffin/N₂O combination differ by a factor of 2–3× between independent studies (see `references.md`). Shipping a single literature $(a, n)$ pair as a tool default would imply a precision the underlying data does not support, and would be exactly the kind of unsourced "magic number" this project's conventions exist to avoid. `grain_sizing.py`'s `FUEL_PROPERTIES` table therefore provides **density only** (a genuine, well-defined material property) plus a literature *reference range* for $a, n$, explicitly labelled as orientation, not a design value — the actual $a, n$ used in a calculation must come from the team's own test data or a source they have specifically vetted for their scale and injector.

## 3b.3 From target O/F to fuel mass flow

Given the oxidiser mass flow $\dot m_{ox}$ (already computed by the rest of this tool) and a target oxidiser-to-fuel mass ratio $OF$ (from thermochemical sizing — CEA, RPA, or equivalent, external to this project), the required fuel mass flow follows directly from the definition of $OF$:

$$\dot m_{fuel} = \frac{\dot m_{ox}}{OF}$$

If the grain has $N_{ports}$ identical circular ports sharing the oxidiser flow equally, both $\dot m_{ox}$ and the target $\dot m_{fuel}$ are simply divided by $N_{ports}$ for the per-port calculation below.

## 3b.4 From fuel mass flow to port radius

The fuel mass flow from one circular port of radius $r$ and length $L$ is the volumetric burn rate of a thin annular shell at the surface, converted to mass via density:

$$\dot m_{fuel} = \rho_{fuel} \cdot \underbrace{(2\pi r)}_{\text{perimeter}} \cdot L \cdot \dot r = \rho_{fuel}\, 2\pi r L\, a\left(\frac{\dot m_{ox}}{\pi r^2}\right)^{\!n}$$

Collecting the powers of $r$:

$$\dot m_{fuel}(r) = K\, r^{\,1-2n}, \qquad K = 2\pi \rho_{fuel} L\, a\, \pi^{-n}\, \dot m_{ox}^{\,n}$$

### The $n=0.5$ special case

When $n = 0.5$, the exponent $1-2n$ is exactly zero: **the port radius drops out of the equation entirely.** Fuel mass flow — and therefore O/F — becomes independent of port geometry. This is a well-known, practically important property of hybrid motors with $n \approx 0.5$ fuels (paraffin among them): the O/F stays comparatively stable as the port regresses over the burn, without needing careful geometric tuning. It is part of why $n \approx 0.5$ fuels are attractive from a design-robustness standpoint, independent of their absolute regression rate.

### The direction-reversal subtlety for $n \neq 0.5$

For $n < 0.5$ (exponent positive), increasing the radius increases fuel flow — the intuitive direction: a bigger hole, more burning perimeter, more fuel. For $n > 0.5$ (exponent **negative**), the reverse happens: increasing the radius *decreases* fuel flow. The reason is that $G_o$ falls as $1/r^2$ while the burning perimeter only grows as $r$; for $n>0.5$, the mass-flux term's decline dominates the perimeter term's growth, and net regression rate — hence mass flow — falls. This is easy to get backwards, and is verified explicitly (both directions, not assumed) in `test_grain_sizing.py`.

For $n \neq 0.5$, the equation inverts in closed form:

$$r_0 = \left(\frac{\dot m_{fuel}}{K}\right)^{\!\frac{1}{1-2n}}$$

`grain_sizing.py` solves this by bisection regardless of $n$ (the same convention this project already uses for other transcendental relationships — saturation temperature inversion, the Henry-Fauske critical flow), which needs no special-casing at $n=0.5$ and gives a natural cross-check against the closed form in testing.

## 3b.5 Multiple ports: a real, simple design lever

For a **fixed total port area**, splitting a single port into $N$ identical circular ports increases total burning perimeter by a factor of $\sqrt N$ (pure geometry — a circle of area $A$ has perimeter $2\sqrt{\pi A}$; $N$ circles of total area $A$, each of area $A/N$, have total perimeter $N \cdot 2\sqrt{\pi A/N} = 2\sqrt{N\pi A}$). More burning perimeter, for the same regression rate, means more fuel mass flow. This is a genuine, widely used hybrid motor design technique, and is simple to model exactly — unlike non-circular single-port shapes (star, wagon-wheel), which change perimeter-to-area relationship as they regress in ways that do not stay self-similar and need dedicated published geometric solutions this project does not yet implement (see `future_work.md`).

## 3b.6 What this module deliberately does not do

- **No transient burn simulation.** The port radius, $G_o$, and regression rate all evolve as the grain burns back; this project computes only the **initial** ($t=0$) sizing point, plus an explicitly-flagged first-order, deliberately conservative estimate of how far the port grows over a given burn duration (holding the *initial* regression rate constant, which over-estimates the true burnback since real $\dot r$ falls as the port grows, for $n>0$ regardless of which side of 0.5 it's on relative to *radius growth* — the point is that $G_o$ itself always falls as the port opens up, for fixed $\dot m_{ox}$). A full transient model — tracking radius, $G_o$, $\dot r$, and hence O/F drift over the whole burn — is `future_work.md` Priority 8.
- **No specific impulse estimate.** $I_{sp}$ requires a chemical equilibrium combustion code (CEA, RPA, or similar); this project does not implement or wrap one. Obtain $I_{sp}$ at the design O/F from CEA/RPA directly.
- **No non-circular port shapes.** See Section 3b.5.
- **Grain length is a required input, not a derived output**, deliberately — a length derived from an unsourced L/D heuristic would be exactly the kind of unjustified default this project avoids; most teams already know their available case length as a hard constraint.

---
*Related documents: [03_two_phase_flow.md](03_two_phase_flow.md) — the injector-side two-phase models this grain sizing step chains onto · [04_implementation.md](04_implementation.md), Section 4.7 — computational implementation of `grain_sizing.py`.*
