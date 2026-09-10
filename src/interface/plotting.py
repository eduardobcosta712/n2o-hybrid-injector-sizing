"""
plotting.py

Interactive Plotly figures for the N2O injector sizing tool (app.py).
All figures support native pan, zoom, and hover. Kept separate from the
Streamlit interface: this module only builds figures from already-computed
results, with no knowledge of Streamlit or of how those results were obtained.

Units: SI throughout (Pa, K), converted to bar/degC for display only.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from n2o_properties import P_sat, T_sat as T_sat_func, T_MIN, T_MAX

# ---------------------------------------------------------------------------
# Colour system — dark theme, production quality
# Backgrounds are transparent so figures blend with the app (no "island" effect)
# Data colours are deliberate and consistent across all figures
# ---------------------------------------------------------------------------
C_LIQUID  = "#60A5FA"    # blue-400  — liquid / subcooled state
C_VAPOR   = "#F87171"    # red-400   — vapour / flashing / warning
C_SAT     = "#94A3B8"    # slate-400 — saturation curve (reference, not data)
C_TANK    = "#60A5FA"    # same as liquid — tank operating point
C_INLET   = "#34D399"    # emerald-400 — injector inlet / safe margin
C_CHAMBER = "#F87171"    # same as vapour — chamber (downstream)
C_HF      = "#FB923C"    # orange-400 — Henry-Fauske non-equilibrium ceiling (diagnostic)
C_GRID    = "rgba(148,163,184,0.12)"   # barely-there grid
C_TEXT    = "#E2E8F0"    # slate-200 — primary text
C_TEXT_2  = "#64748B"    # slate-500 — secondary text / axis labels

PLOTLY_LAYOUT = dict(
    paper_bgcolor="rgba(0,0,0,0)",   # fully transparent — no island
    plot_bgcolor="rgba(15,23,42,0.6)",  # very dark slate, semi-transparent
    font=dict(family="Inter, sans-serif", color=C_TEXT, size=12),
    margin=dict(l=58, r=20, t=42, b=50),
    legend=dict(
        bgcolor="rgba(15,23,42,0.85)",
        bordercolor="rgba(148,163,184,0.15)",
        borderwidth=1,
        font=dict(size=11, color=C_TEXT),
    ),
    xaxis=dict(
        gridcolor=C_GRID,
        zerolinecolor="rgba(148,163,184,0.2)",
        linecolor="rgba(148,163,184,0.15)",
        tickfont=dict(color=C_TEXT_2, size=11),
        title_font=dict(color=C_TEXT_2, size=12),
        showgrid=True,
        zeroline=False,
    ),
    yaxis=dict(
        gridcolor=C_GRID,
        zerolinecolor="rgba(148,163,184,0.2)",
        linecolor="rgba(148,163,184,0.15)",
        tickfont=dict(color=C_TEXT_2, size=11),
        title_font=dict(color=C_TEXT_2, size=12),
        showgrid=True,
        zeroline=False,
    ),
)


def plot_pressure_along_line(trace, P_tank, T_tank, segments):
    """
    Interactive Plotly figure: pressure at each point along the feed line,
    alongside P_sat(T_tank). X-axis is cumulative pipe length in metres;
    fittings are point-losses shown at the same x as the preceding segment.
    """
    x_positions = [0.0]
    x_current = 0.0
    for seg in segments:
        if seg["type"] == "pipe":
            x_current += seg["L"]
        x_positions.append(x_current)

    P_bar = [P_tank / 1e5] + [s["pressure_after_Pa"] / 1e5 for s in trace]
    P_sat_val = P_sat(T_tank) / 1e5

    labels = ["Tank exit"] + [
        f"After seg {s['segment_index']} ({s['segment_type']})" for s in trace
    ]
    hover = [
        f"<b>{lab}</b><br>P = {p:.3f} bar<br>ΔT_sub = {s['delta_T_sub_K']:.2f} K"
        if i > 0 else f"<b>{lab}</b><br>P = {p:.3f} bar"
        for i, (lab, p, s) in enumerate(
            zip(labels, P_bar, [None] + list(trace)))
    ]

    fig = go.Figure()

    # Subcooling margin shading
    fig.add_trace(go.Scatter(
        x=x_positions, y=P_bar,
        fill="tonexty", fillcolor="rgba(96,165,250,0.08)",
        line=dict(color="rgba(0,0,0,0)"), showlegend=False, hoverinfo="skip"
    ))

    # P_sat reference line
    fig.add_hline(y=P_sat_val, line=dict(color=C_CHAMBER, dash="dash", width=1.5),
                  annotation_text=f"P_sat({T_tank-273.15:.1f}°C) = {P_sat_val:.2f} bar",
                  annotation_font_color=C_CHAMBER, annotation_font_size=11)

    # Fluid pressure trace
    fig.add_trace(go.Scatter(
        x=x_positions, y=P_bar,
        mode="lines+markers",
        line=dict(color=C_TANK, width=2.5),
        marker=dict(size=7, color=C_TANK, line=dict(color="white", width=1)),
        name="Fluid pressure",
        hovertemplate="%{customdata}<extra></extra>",
        customdata=hover,
    ))

    fig.update_layout(**PLOTLY_LAYOUT,
        title=dict(text="Pressure along the feed line", font=dict(size=13, color=C_TEXT)),
        xaxis_title="Position along line (m)",
        yaxis_title="Pressure (bar)",
    )
    return fig


def plot_PT_diagram(T_tank, P_tank, P_injector_inlet, P_chamber):
    """
    Interactive Plotly P-T diagram: saturation curve with liquid/vapour
    regions shaded, and the three operating points marked. The chamber
    point is plotted at T_sat(P_chamber), since in the HEM/Dyer model
    the exiting mixture is in equilibrium on the saturation curve.
    Native Plotly zoom/pan replaces the static zoom checkbox.
    """
    T_curve = np.linspace(T_MIN, T_MAX, 300)
    P_curve_bar = [P_sat(T) / 1e5 for T in T_curve]
    T_curve_C = T_curve - 273.15

    T_tank_C = T_tank - 273.15
    T_chamber_C = T_sat_func(P_chamber) - 273.15

    fig = go.Figure()

    # Liquid region fill
    fig.add_trace(go.Scatter(
        x=T_curve_C, y=[p + (80 - p) for p in P_curve_bar],
        fill=None, line=dict(color="rgba(0,0,0,0)"),
        showlegend=False, hoverinfo="skip"
    ))
    fig.add_trace(go.Scatter(
        x=T_curve_C, y=P_curve_bar,
        fill="tonexty", fillcolor="rgba(96,165,250,0.08)",
        line=dict(color="rgba(0,0,0,0)"), name="Liquid region", hoverinfo="skip"
    ))

    # Vapour region fill
    fig.add_trace(go.Scatter(
        x=T_curve_C, y=P_curve_bar,
        fill="tozeroy", fillcolor="rgba(232,137,74,0.10)",
        line=dict(color="rgba(0,0,0,0)"), name="Vapour / two-phase", hoverinfo="skip"
    ))

    # Saturation curve
    fig.add_trace(go.Scatter(
        x=T_curve_C, y=P_curve_bar,
        mode="lines", line=dict(color="white", width=2),
        name="Saturation curve",
        hovertemplate="T = %{x:.1f}°C<br>P_sat = %{y:.2f} bar<extra></extra>"
    ))

    # Arrow (dashed line tank→chamber)
    fig.add_trace(go.Scatter(
        x=[T_tank_C, T_chamber_C], y=[P_tank/1e5, P_chamber/1e5],
        mode="lines", line=dict(color="rgba(200,200,200,0.35)", dash="dot", width=1.5),
        showlegend=False, hoverinfo="skip"
    ))

    # Operating points
    points = [
        (T_tank_C,    P_tank/1e5,           C_TANK,    "square",   "Tank",
         f"T = {T_tank_C:.1f}°C<br>P = {P_tank/1e5:.1f} bar"),
        (T_tank_C,    P_injector_inlet/1e5,  C_INLET,   "triangle-up","Injector inlet",
         f"T = {T_tank_C:.1f}°C<br>P = {P_injector_inlet/1e5:.2f} bar"),
        (T_chamber_C, P_chamber/1e5,         C_CHAMBER, "triangle-down","Chamber (T_sat)",
         f"T_sat = {T_chamber_C:.1f}°C<br>P = {P_chamber/1e5:.1f} bar"),
    ]
    for tx, py, col, sym, name, tip in points:
        fig.add_trace(go.Scatter(
            x=[tx], y=[py], mode="markers",
            marker=dict(symbol=sym, size=13, color=col,
                        line=dict(color="white", width=1.5)),
            name=name,
            hovertemplate=f"<b>{name}</b><br>{tip}<extra></extra>"
        ))

    fig.update_layout(**{k: v for k, v in PLOTLY_LAYOUT.items() if k not in ("xaxis","yaxis")},
        title=dict(text="P-T diagram — operating points vs. saturation curve",
                   font=dict(size=13, color=C_TEXT)),
        xaxis=dict(title="Temperature (°C)", gridcolor=C_GRID, zerolinecolor=C_GRID),
        yaxis=dict(title="Pressure (bar)", range=[0, 80], gridcolor=C_GRID,
                   zerolinecolor=C_GRID),
    )
    return fig


def plot_model_comparison(m_spi, m_hem, m_dyer, m_target, m_dot_crit_HF=None,
                          choked=False):
    """
    Bar chart comparing SPI, HEM, and Dyer mass flow predictions against
    the design target. Specific to Design mode: shows the magnitude of
    the two-phase correction and the spread between model limits.

    Parameters
    ----------
    m_spi, m_hem, m_dyer, m_target : float
        Mass flow rates, kg/s, as before.
    m_dot_crit_HF : float or None, optional
        Henry-Fauske non-equilibrium critical flow ceiling, kg/s (added
        September 2026, docs/future_work.md Priority 1). When provided,
        drawn as a reference line -- DIAGNOSTIC ONLY, since it is not
        applied to m_dyer itself (see injector_two_phase.dyer_mass_flow
        docstring for why). Omit or pass None to hide it (e.g. when the
        ceiling could not be computed -- see "HF_unavailable_reason" in
        the Dyer result).
    choked : bool, optional
        Whether m_dyer exceeds m_dot_crit_HF. Only affects styling (the
        ceiling line is drawn in a warning colour when True). Ignored if
        m_dot_crit_HF is None.
    """
    labels = ["SPI", "HEM", "Dyer (adopted)", "Target"]
    values = [m_spi * 1000, m_hem * 1000, m_dyer * 1000, m_target * 1000]
    colours = ["#3d5a80", "#2a4a3a", "#2d6a4f", "rgba(200,200,200,0.15)"]
    text = [f"{v:.1f} g/s" for v in values]

    fig = go.Figure()
    for lab, val, col, txt in zip(labels, values, colours, text):
        fig.add_trace(go.Bar(
            name=lab, x=[lab], y=[val],
            marker_color=col,
            marker_line=dict(color="rgba(255,255,255,0.08)", width=1),
            text=txt, textposition="outside",
            textfont=dict(color="#8090a0", size=11,
                          family="JetBrains Mono, monospace"),
            hovertemplate=f"<b>{lab}</b><br>{val:.1f} g/s<extra></extra>",
        ))

    # Target line
    fig.add_hline(y=m_target * 1000,
                  line=dict(color="rgba(180,180,180,0.3)", dash="dot", width=1.5),
                  annotation_text=f"Target: {m_target*1000:.1f} g/s",
                  annotation_font_color="rgba(150,150,150,0.7)",
                  annotation_font_size=10)

    # Henry-Fauske non-equilibrium ceiling (diagnostic only -- see
    # docstring). Drawn as a distinct dashed line, orange normally, red
    # if the Dyer bar actually exceeds it (choked=True), so the chart
    # itself flags the same condition the text warning describes.
    if m_dot_crit_HF is not None:
        hf_colour = C_CHAMBER if choked else C_HF
        fig.add_hline(
            y=m_dot_crit_HF * 1000,
            line=dict(color=hf_colour, dash="dashdot", width=1.5),
            annotation_text=(
                f"Henry-Fauske ceiling (diagnostic): {m_dot_crit_HF*1000:.1f} g/s"
                + (" — EXCEEDED" if choked else "")
            ),
            annotation_font_color=hf_colour,
            annotation_font_size=10,
            annotation_position="top left",
        )

    fig.update_layout(
        **{k: v for k, v in PLOTLY_LAYOUT.items() if k not in ("xaxis", "yaxis")},
        title=dict(text="Model comparison — mass flow predictions",
                   font=dict(size=13, color=C_TEXT)),
        xaxis=dict(title="", gridcolor=C_GRID, zerolinecolor=C_GRID),
        yaxis=dict(title="Mass flow (g/s)", gridcolor=C_GRID,
                   zerolinecolor=C_GRID),
        showlegend=False,
        bargap=0.35,
    )
    return fig


def plot_line_profile(segments_ui):
    """
    SVG-style Plotly diagram showing the feed line in profile:
    pipes as horizontal rectangles (length proportional to L),
    fittings as standard schematic symbols.

    Parameters
    ----------
    segments_ui : list of dict
        The raw UI segment list from st.session_state.segments
        (with L_mm, D_mm, fitting_name, K keys).

    Returns
    -------
    plotly.graph_objects.Figure
    """
    # Layout constants
    PIPE_H   = 0.18   # pipe height in y-units
    FIT_W    = 0.4    # fitting width in x-units (fixed)
    Y_CL     = 0.5    # centreline y
    SCALE    = 4.0    # x-scale: 1 m of pipe = SCALE x-units

    shapes = []
    annotations = []
    x_cur = 0.0

    for seg in segments_ui:
        if seg["type"] == "pipe":
            L = seg.get("L_mm", 1000.0) / 1000.0
            w = L * SCALE
            # Pipe rectangle
            shapes.append(dict(
                type="rect",
                x0=x_cur, x1=x_cur + w,
                y0=Y_CL - PIPE_H / 2, y1=Y_CL + PIPE_H / 2,
                fillcolor="rgba(40,70,110,0.7)",
                line=dict(color="rgba(74,144,217,0.5)", width=1),
            ))
            # Label: length
            annotations.append(dict(
                x=x_cur + w / 2, y=Y_CL + PIPE_H / 2 + 0.07,
                text=f"{seg['L_mm']:.0f} mm",
                showarrow=False,
                font=dict(size=9, color="#4a6fa5",
                          family="JetBrains Mono, monospace"),
                xanchor="center",
            ))
            # Label: diameter
            annotations.append(dict(
                x=x_cur + w / 2, y=Y_CL - PIPE_H / 2 - 0.07,
                text=f"ID {seg['D_mm']:.1f} mm",
                showarrow=False,
                font=dict(size=8, color="#3d4a5c",
                          family="JetBrains Mono, monospace"),
                xanchor="center",
            ))
            x_cur += w

        else:
            fname = seg.get("fitting_name", "fitting")
            K = seg.get("K", 0)
            # Fitting symbol: diamond shape
            cx = x_cur + FIT_W / 2
            shapes.append(dict(
                type="path",
                path=(f"M {x_cur},{Y_CL} "
                      f"L {cx},{Y_CL + PIPE_H * 0.8} "
                      f"L {x_cur + FIT_W},{Y_CL} "
                      f"L {cx},{Y_CL - PIPE_H * 0.8} Z"),
                fillcolor="rgba(60,40,80,0.75)",
                line=dict(color="rgba(130,90,180,0.6)", width=1),
            ))
            # Short name for label
            short = fname.split("(")[0].strip()
            if len(short) > 14:
                short = short[:14]
            annotations.append(dict(
                x=cx, y=Y_CL + PIPE_H * 0.8 + 0.07,
                text=short,
                showarrow=False,
                font=dict(size=8, color="#7060a0",
                          family="Inter, sans-serif"),
                xanchor="center",
            ))
            annotations.append(dict(
                x=cx, y=Y_CL - PIPE_H * 0.8 - 0.07,
                text=f"K={K:.2f}",
                showarrow=False,
                font=dict(size=8, color="#3d4a5c",
                          family="JetBrains Mono, monospace"),
                xanchor="center",
            ))
            x_cur += FIT_W

    # Tank symbol on the left
    shapes.append(dict(
        type="rect",
        x0=-0.35, x1=0.0,
        y0=Y_CL - 0.3, y1=Y_CL + 0.3,
        fillcolor="rgba(30,50,80,0.8)",
        line=dict(color="rgba(74,144,217,0.4)", width=1.5),
    ))
    annotations.append(dict(
        x=-0.175, y=Y_CL,
        text="Tank", showarrow=False,
        font=dict(size=9, color="#4a6fa5", family="Inter, sans-serif"),
        xanchor="center",
    ))

    # Injector symbol on the right
    shapes.append(dict(
        type="path",
        path=(f"M {x_cur},{Y_CL+0.22} "
              f"L {x_cur+0.28},{Y_CL+0.08} "
              f"L {x_cur+0.28},{Y_CL-0.08} "
              f"L {x_cur},{Y_CL-0.22} Z"),
        fillcolor="rgba(40,80,50,0.75)",
        line=dict(color="rgba(80,160,100,0.5)", width=1.5),
    ))
    annotations.append(dict(
        x=x_cur + 0.14, y=Y_CL + 0.3,
        text="Injector", showarrow=False,
        font=dict(size=9, color="#2d6a4f", family="Inter, sans-serif"),
        xanchor="center",
    ))

    x_max = x_cur + 0.5
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers",
                             showlegend=False, hoverinfo="skip"))
    fig.update_layout(
        **{k: v for k, v in PLOTLY_LAYOUT.items()
           if k not in ("xaxis", "yaxis", "margin")},
        title=dict(text="Feed line schematic", font=dict(size=13, color=C_TEXT)),
        shapes=shapes,
        annotations=annotations,
        xaxis=dict(range=[-0.5, x_max], showgrid=False,
                   zeroline=False, showticklabels=False,
                   gridcolor=C_GRID),
        yaxis=dict(range=[0.0, 1.0], showgrid=False,
                   zeroline=False, showticklabels=False,
                   scaleanchor="x", scaleratio=1,
                   gridcolor=C_GRID),
        margin=dict(l=10, r=10, t=40, b=10),
        height=200,
    )
    return fig


def plot_subcooling_margin(trace, P_tank, T_tank, segments):
    """
    Plot the subcooling margin (delta_T_sub) at each point along the
    feed line -- how many degrees of margin remain before the fluid
    reaches saturation. Complements the pressure chart by showing
    the thermal safety margin directly.
    """
    x_positions = [0.0]
    x_cur = 0.0
    for seg in segments:
        if seg["type"] == "pipe":
            x_cur += seg["L"]
        x_positions.append(x_cur)

    # delta_T_sub at tank exit
    from n2o_properties import T_sat as T_sat_fn
    dT_tank = T_sat_fn(P_tank) - T_tank
    dT_vals = [dT_tank] + [s["delta_T_sub_K"] for s in trace]
    labels = ["Tank exit"] + [
        f"After seg {s['segment_index']} ({s['segment_type']})" for s in trace
    ]

    fig = go.Figure()

    # Zero line (saturation threshold)
    fig.add_hline(y=0, line=dict(color=C_CHAMBER, dash="dash", width=1.5),
                  annotation_text="Saturation (flashing threshold)",
                  annotation_font_color=C_CHAMBER, annotation_font_size=10)

    # Shading: positive = subcooled (good), negative = flashing
    fig.add_trace(go.Scatter(
        x=x_positions, y=[max(0, v) for v in dT_vals],
        fill="tozeroy", fillcolor="rgba(45,106,79,0.12)",
        line=dict(color="rgba(0,0,0,0)"), showlegend=False, hoverinfo="skip"))
    fig.add_trace(go.Scatter(
        x=x_positions, y=[min(0, v) for v in dT_vals],
        fill="tozeroy", fillcolor="rgba(122,32,32,0.18)",
        line=dict(color="rgba(0,0,0,0)"), showlegend=False, hoverinfo="skip"))

    hover = [f"<b>{lab}</b><br>dT_sub = {v:.2f} K"
             for lab, v in zip(labels, dT_vals)]
    fig.add_trace(go.Scatter(
        x=x_positions, y=dT_vals,
        mode="lines+markers",
        line=dict(color=C_INLET, width=2.5),
        marker=dict(size=7, color=C_INLET, line=dict(color="white", width=1)),
        name="Subcooling margin",
        hovertemplate="%{customdata}<extra></extra>",
        customdata=hover,
    ))

    fig.update_layout(**PLOTLY_LAYOUT,
        title=dict(text="Subcooling margin along the feed line",
                   font=dict(size=13, color=C_TEXT)),
        xaxis_title="Position along line (m)",
        yaxis_title="Delta T_sub (K)",
    )
    return fig


def plot_sensitivity(T_tank, P_tank, P_chamber, model_segments, Cd, A_injector,
                     roughness=1.5e-6):
    """
    Sensitivity chart: how does the real mass flow change as tank
    temperature varies +/-5 deg C? Shows the design point and the
    sensitivity envelope, helping the user understand launch-condition risk.
    """
    import sys, os
    _MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model")
    sys.path.insert(0, os.path.abspath(_MODEL_DIR))
    from full_system import evaluate_full_system
    from n2o_properties import P_sat as P_sat_fn

    T_range = np.linspace(T_tank - 5, T_tank + 5, 22)
    m_vals, valid, flash_temps = [], [], []

    for T in T_range:
        if T < 182.33 + 273.15 - 273.15 or T > 309.52:
            continue
        try:
            r = evaluate_full_system(0.3, T, P_tank, model_segments,
                                     Cd, A_injector, P_chamber, roughness)
            m = r.get("m_dot_real")
            if r["feed_line_result"]["flashing_detected"]:
                flash_temps.append(T - 273.15)
            # Use m_dot_real whether or not flashing — model now always returns a value
            m_vals.append(m * 1000 if m is not None else None)
            valid.append(T - 273.15)
        except Exception:
            pass

    T_vals_C = [T - 273.15 for T in T_range
                if 182.33 <= T <= 309.52]

    fig = go.Figure()

    # Main sensitivity curve
    m_clean = [m for m in m_vals if m is not None]
    T_clean = [t for t, m in zip(valid, m_vals) if m is not None]

    if T_clean:
        fig.add_trace(go.Scatter(
            x=T_clean, y=m_clean,
            mode="lines+markers",
            line=dict(color=C_TANK, width=2),
            marker=dict(size=5),
            name="Real mass flow",
            hovertemplate="T = %{x:.1f}°C<br>m_dot = %{y:.1f} g/s<extra></extra>",
        ))

    # Design point marker
    try:
        r_dp = evaluate_full_system(0.3, T_tank, P_tank, model_segments,
                                     Cd, A_injector, P_chamber, roughness)
        if r_dp["m_dot_real"]:
            fig.add_trace(go.Scatter(
                x=[T_tank - 273.15], y=[r_dp["m_dot_real"] * 1000],
                mode="markers",
                marker=dict(symbol="star", size=14, color=C_CHAMBER,
                            line=dict(color="white", width=1.5)),
                name="Design point",
                hovertemplate=f"Design point<br>T = {T_tank-273.15:.1f}°C<br>"
                              f"m_dot = {r_dp['m_dot_real']*1000:.1f} g/s<extra></extra>",
            ))
    except Exception:
        pass

    # Flashing zone shading
    if flash_temps:
        fig.add_vrect(x0=min(flash_temps), x1=max(flash_temps) + 0.5,
                      fillcolor="rgba(122,32,32,0.12)",
                      line_width=0,
                      annotation_text="Flashing zone",
                      annotation_font_color=C_CHAMBER,
                      annotation_font_size=9)

    fig.update_layout(**PLOTLY_LAYOUT,
        title=dict(text="Sensitivity: mass flow vs tank temperature (+/-5 deg C)",
                   font=dict(size=13, color=C_TEXT)),
        xaxis_title="Tank temperature (deg C)",
        yaxis_title="Real mass flow (g/s)",
    )
    return fig


def plot_segment_losses(trace, segments):
    """
    Horizontal bar chart showing the pressure drop contributed by each
    segment, making it immediately clear where the feed line losses
    are concentrated.
    """
    labels, drops, types = [], [], []
    for s in trace:
        seg = segments[s["segment_index"]] if s["segment_index"] < len(segments) else {}
        t = s["segment_type"]
        if t == "pipe":
            lbl = f"Pipe {s['segment_index']+1} ({seg.get('L',0)*1000:.0f}mm)"
        else:
            lbl = f"Fitting {s['segment_index']+1}"
        labels.append(lbl)
        drops.append(s["pressure_drop_Pa"] / 1e5)
        types.append(t)

    colours = [C_TANK if t == "pipe" else "rgba(130,90,180,0.8)"
               for t in types]

    fig = go.Figure(go.Bar(
        y=labels, x=drops,
        orientation="h",
        marker_color=colours,
        marker_line=dict(color="rgba(255,255,255,0.08)", width=0.5),
        hovertemplate="%{y}<br>dP = %{x:.4f} bar<extra></extra>",
    ))
    fig.update_layout(
        **{k: v for k, v in PLOTLY_LAYOUT.items() if k not in ("xaxis","yaxis")},
        title=dict(text="Pressure drop by segment", font=dict(size=13, color=C_TEXT)),
        xaxis=dict(title="Pressure drop (bar)", gridcolor=C_GRID,
                   zerolinecolor=C_GRID),
        yaxis=dict(title="", autorange="reversed",
                   gridcolor=C_GRID, zerolinecolor=C_GRID),
        bargap=0.25,
    )
    return fig


def plot_tornado(T_tank, P_tank, P_chamber, model_segments, Cd, A_injector,
                 roughness=1.5e-6):
    """
    Tornado plot: one-at-a-time sensitivity of real mass flow to each
    uncertain input, at the current design point.

    Inputs varied:
        T_tank     : +/-5 K         -- ambient day variation
        P_tank     : +/-5%          -- regulator / gauge uncertainty
        Cd         : +/-0.05 (abs)  -- inlet geometry uncertainty
        D_orifice  : +/-0.05 mm     -- machining tolerance
        Line length: +/-20%         -- installation uncertainty
    """
    import sys, os, math
    _DIR = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.abspath(os.path.join(_DIR, '..', 'model')))
    from full_system import evaluate_full_system

    def run(T, P, cd, A, segs):
        try:
            r = evaluate_full_system(0.3, T, P, segs, cd, A, P_chamber, roughness)
            m = r.get('m_dot_real')
            return float(m) * 1000 if m is not None else None
        except Exception:
            return None

    m_nom = run(T_tank, P_tank, Cd, A_injector, model_segments)
    if m_nom is None:
        fig = go.Figure()
        fig.add_annotation(text="Nominal point unavailable",
                           xref="paper", yref="paper", x=0.5, y=0.5,
                           showarrow=False, font=dict(color=C_TEXT))
        fig.update_layout(**PLOTLY_LAYOUT,
                          title=dict(text="Sensitivity tornado",
                                     font=dict(size=13, color=C_TEXT)))
        return fig

    D_equiv = math.sqrt(4.0 * A_injector / math.pi)

    def scale_pipes(segs, f):
        return [{**s, 'L': s['L'] * f} if s.get('type') == 'pipe' else dict(s)
                for s in segs]

    # (label, m_low, m_high) — note: temp bug fix: use +/-5 K not +/-5%
    candidates = [
        ("Tank temperature (±5 K)",
         run(T_tank - 5.0, P_tank, Cd, A_injector, model_segments),
         run(T_tank + 5.0, P_tank, Cd, A_injector, model_segments)),
        ("Tank pressure (±5%)",
         run(T_tank, P_tank * 0.95, Cd, A_injector, model_segments),
         run(T_tank, P_tank * 1.05, Cd, A_injector, model_segments)),
        ("Discharge coeff. Cd (±0.05)",
         run(T_tank, P_tank, max(0.3, Cd - 0.05), A_injector, model_segments),
         run(T_tank, P_tank, min(1.0, Cd + 0.05), A_injector, model_segments)),
        ("Orifice diameter (±0.05 mm)",
         run(T_tank, P_tank, Cd,
             math.pi * (max(1e-4, D_equiv - 0.05e-3) / 2) ** 2, model_segments),
         run(T_tank, P_tank, Cd,
             math.pi * ((D_equiv + 0.05e-3) / 2) ** 2, model_segments)),
        ("Feed line length (±20%)",
         run(T_tank, P_tank, Cd, A_injector, scale_pipes(model_segments, 0.80)),
         run(T_tank, P_tank, Cd, A_injector, scale_pipes(model_segments, 1.20))),
    ]

    rows = [(lbl, lo, hi) for lbl, lo, hi in candidates
            if lo is not None and hi is not None]
    if not rows:
        fig = go.Figure()
        fig.update_layout(**PLOTLY_LAYOUT)
        return fig

    # Sort by total swing
    rows.sort(key=lambda r: abs(r[2] - r[1]), reverse=True)
    labels = [r[0] for r in rows]
    m_lo   = [r[1] for r in rows]
    m_hi   = [r[2] for r in rows]

    # Deltas from nominal — pre-formatted as strings to avoid Plotly float
    # formatting issues in go.Bar hovertemplate (%{x} does not reliably
    # apply format strings to floats with many decimal places)
    d_lo = [lo - m_nom for lo in m_lo]
    d_hi = [hi - m_nom for hi in m_hi]

    # customdata: [delta_str, m_dot_str] for each bar
    cd_lo = [[f"{d:+.1f}", f"{m:.0f}"] for d, m in zip(d_lo, m_lo)]
    cd_hi = [[f"{d:+.1f}", f"{m:.0f}"] for d, m in zip(d_hi, m_hi)]

    fig = go.Figure()

    # Low-side bar (to the left of nominal)
    fig.add_trace(go.Bar(
        name="Lower bound",
        y=labels, x=d_lo,
        orientation="h",
        marker=dict(
            color="rgba(96,165,250,0.75)",
            line=dict(color="rgba(96,165,250,0.3)", width=0.5),
        ),
        hovertemplate="<b>%{y}</b><br>Δm = %{customdata[0]} g/s<br>m_dot = %{customdata[1]} g/s<extra></extra>",
        customdata=cd_lo,
    ))

    # High-side bar (to the right of nominal)
    fig.add_trace(go.Bar(
        name="Upper bound",
        y=labels, x=d_hi,
        orientation="h",
        marker=dict(
            color="rgba(248,113,113,0.75)",
            line=dict(color="rgba(248,113,113,0.3)", width=0.5),
        ),
        hovertemplate="<b>%{y}</b><br>Δm = %{customdata[0]} g/s<br>m_dot = %{customdata[1]} g/s<extra></extra>",
        customdata=cd_hi,
    ))

    fig.add_vline(x=0, line=dict(color="rgba(226,232,240,0.5)", width=1.2, dash="dot"))
    fig.add_annotation(
        x=0, y=1.07, xref="x", yref="paper",
        text=f"Nominal: {m_nom:.0f} g/s",
        showarrow=False, xanchor="center",
        font=dict(size=11, color=C_TEXT_2),
    )

    base = {k: v for k, v in PLOTLY_LAYOUT.items()
            if k not in ("xaxis", "yaxis", "legend")}
    fig.update_layout(
        **base,
        title=dict(text="Sensitivity tornado — Δm_dot per input",
                   font=dict(size=13, color=C_TEXT)),
        barmode="overlay",
        bargap=0.30,
        xaxis=dict(
            title="Δ mass flow from nominal (g/s)",
            gridcolor=C_GRID,
            zerolinecolor="rgba(226,232,240,0.3)",
            tickfont=dict(color=C_TEXT_2, size=11),
            title_font=dict(color=C_TEXT_2, size=12),
        ),
        yaxis=dict(
            title="",
            autorange="reversed",
            gridcolor=C_GRID,
            tickfont=dict(color=C_TEXT, size=11),
        ),
        legend=dict(
            bgcolor="rgba(15,23,42,0.85)",
            bordercolor="rgba(148,163,184,0.15)",
            borderwidth=1,
            font=dict(size=11, color=C_TEXT),
            orientation="h", x=0.5, y=-0.18, xanchor="center",
        ),
        showlegend=True,
    )
    return fig
