"""
app.py — N2O Hybrid Rocket Injector Sizing Tool (Streamlit interface).

Two modes:
  Sizing — known orifice area -> predicted real mass flow (SPI / Dyer)
  Design — target mass flow   -> required orifice area (SPI + Dyer correction)
"""

import os, sys, math
_MODEL_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "model")
sys.path.insert(0, os.path.abspath(_MODEL_DIR))

import streamlit as st
import json
from n2o_properties import P_sat, T_sat, rho_liquid_sat, nu_vapor_sat
from feed_line import evaluate_feed_line
from injector_spi import spi_mass_flow, spi_sufficient, orifice_area_from_target_flow
from injector_two_phase import dyer_mass_flow
from full_system import evaluate_full_system
from plotting import (plot_pressure_along_line, plot_PT_diagram, plot_model_comparison,
                      plot_line_profile, plot_subcooling_margin, plot_sensitivity,
                      plot_segment_losses)
from export import generate_pdf

M_N2O = 44.013

st.set_page_config(page_title="N2O Injector Sizing", page_icon=None, layout="wide")

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600&family=JetBrains+Mono:wght@400;500&display=swap');
html,body,[class*="css"]{font-family:'Inter',sans-serif;}
.stApp{background:#0d1117;}
.hero-title{font-size:1.75rem;font-weight:600;color:#e0e6f0;margin-bottom:0.1rem;}
.hero-byline{color:#3d4a5c;font-size:0.75rem;letter-spacing:.1em;text-transform:uppercase;margin-bottom:1.5rem;}
.hero-desc{color:#7a8494;font-size:.9rem;line-height:1.7;max-width:700px;margin-bottom:1.5rem;}
.mode-card{background:#111827;border:1px solid #1f2937;border-radius:5px;padding:1.3rem 1.5rem;margin-bottom:.5rem;}
.mode-card h3{color:#b0bac8;margin:0 0 .3rem;font-size:.9rem;font-weight:600;letter-spacing:.03em;}
.mode-card p{color:#5a6475;margin:0;font-size:.82rem;line-height:1.55;}
.assumption-box{background:#0f0f18;border:1px solid #1e1a30;border-left:3px solid #3d3560;border-radius:3px;padding:.9rem 1.1rem;margin-top:1.2rem;}
.assumption-box h4{color:#7060a0;margin:0 0 .45rem;font-size:.65rem;letter-spacing:.12em;text-transform:uppercase;}
.assumption-box ul{color:#5a6475;margin:0;padding-left:1rem;font-size:.8rem;line-height:1.9;}
.disclaimer-box{background:#12080a;border:1px solid #2a1018;border-left:3px solid #6a2030;border-radius:3px;padding:.9rem 1.1rem;margin-top:.8rem;}
.disclaimer-box h4{color:#8a3040;margin:0 0 .4rem;font-size:.65rem;letter-spacing:.12em;text-transform:uppercase;}
.disclaimer-box p{color:#5a4048;margin:0;font-size:.8rem;line-height:1.7;}
.section-label{font-size:.62rem;font-weight:600;letter-spacing:.14em;color:#3d5a80;text-transform:uppercase;margin:1rem 0 .2rem;border-bottom:1px solid #1a2030;padding-bottom:.2rem;}
.param-help{font-size:.74rem;color:#3d4a5c;margin-top:-.3rem;margin-bottom:.45rem;line-height:1.5;}
.result-card{background:#0f1520;border:1px solid #1a2235;border-radius:4px;padding:.85rem .9rem;text-align:center;}
.result-card .label{font-size:.62rem;color:#3d4a5c;letter-spacing:.1em;text-transform:uppercase;}
.result-card .value{font-size:1.4rem;font-weight:500;color:#c0cad8;font-family:'JetBrains Mono',monospace;margin:.1rem 0;}
.result-card .sub{font-size:.72rem;color:#3d4a5c;}
.badge-ok{display:inline-block;background:#0a1a12;color:#2d8a5a;border:1px solid #1a3028;border-radius:3px;padding:.12rem .55rem;font-size:.75rem;font-family:'JetBrains Mono',monospace;}
.badge-warn{display:inline-block;background:#1a0d0d;color:#9a3030;border:1px solid #301515;border-radius:3px;padding:.12rem .55rem;font-size:.75rem;font-family:'JetBrains Mono',monospace;}
.badge-info{display:inline-block;background:#0a1525;color:#3d5a80;border:1px solid #152035;border-radius:3px;padding:.12rem .55rem;font-size:.75rem;font-family:'JetBrains Mono',monospace;}
.diag-box{background:#100a0a;border:1px solid #201010;border-left:3px solid #602020;border-radius:3px;padding:.9rem 1.2rem;margin-top:.9rem;}
.diag-box h4{color:#803030;margin:0 0 .45rem;font-size:.62rem;letter-spacing:.12em;text-transform:uppercase;}
.diag-box ul{color:#7a8494;margin:0;padding-left:1rem;font-size:.82rem;line-height:1.9;}
.estimate-box{background:#0f100a;border:1px solid #202010;border-left:3px solid #505020;border-radius:3px;padding:.9rem 1.2rem;margin-top:.6rem;}
.estimate-box h4{color:#707030;margin:0 0 .35rem;font-size:.62rem;letter-spacing:.12em;text-transform:uppercase;}
.estimate-box p{color:#7a8494;margin:0;font-size:.8rem;line-height:1.65;}
.insert-row{text-align:center;padding:.1rem 0;opacity:.4;transition:opacity .15s;}
.insert-row:hover{opacity:1;}
.divider{border:none;border-top:1px solid #1a2030;margin:1.3rem 0;}
.block-container{padding-top:1.8rem;padding-bottom:1.8rem;}
div[data-testid="stSidebarContent"]{background:#0a0e17;border-right:1px solid #151c28;}
</style>
""", unsafe_allow_html=True)

# ── Fitting K table ───────────────────────────────────────────────────────────
FITTING_K = {
    "Ball valve (fully open)":    0.05,
    "Ball valve (1/2 open)":      5.5,
    "Needle valve (fully open)":  2.0,
    "Globe valve (fully open)":   10.0,
    "Check valve":                2.5,
    "Elbow 90 deg (standard)":    0.9,
    "Elbow 90 deg (long radius)": 0.4,
    "Elbow 45 deg":               0.4,
    "Tee (flow-through)":         0.4,
    "Tee (branch)":               1.5,
    "Union / coupling":           0.04,
    "Custom (enter K manually)":  None,
}

# ── Session state ─────────────────────────────────────────────────────────────
if "page" not in st.session_state:
    st.session_state.page = "landing"
if "segments" not in st.session_state:
    st.session_state.segments = [
        {"type": "pipe",    "L_mm": 1000.0, "D_mm": 8.0},
        {"type": "fitting", "D_mm": 8.0,
         "fitting_name": "Ball valve (fully open)", "K": 0.05},
        {"type": "pipe",    "L_mm": 500.0,  "D_mm": 8.0},
    ]

# Presets — generic design scenarios, no real-world team references
PRESETS = {
    "Custom (current settings)": None,
    "Compact lab motor (low pressure)": {
        "T_tank_C": 15.0, "P_tank_bar": 45.0, "Cd": 0.65,
        "P_chamber_bar": 12.0,
        "segments": [
            {"type": "pipe", "L_mm": 500.0, "D_mm": 8.0},
            {"type": "fitting", "D_mm": 8.0,
             "fitting_name": "Ball valve (fully open)", "K": 0.05},
        ]
    },
    "Mid-scale hybrid (moderate pressure)": {
        "T_tank_C": 20.0, "P_tank_bar": 55.0, "Cd": 0.65,
        "P_chamber_bar": 20.0,
        "segments": [
            {"type": "pipe", "L_mm": 1000.0, "D_mm": 10.0},
            {"type": "fitting", "D_mm": 10.0,
             "fitting_name": "Ball valve (fully open)", "K": 0.05},
            {"type": "pipe", "L_mm": 500.0, "D_mm": 10.0},
            {"type": "fitting", "D_mm": 10.0,
             "fitting_name": "Elbow 90 deg (standard)", "K": 0.9},
        ]
    },
    "High-pressure research motor": {
        "T_tank_C": 10.0, "P_tank_bar": 65.0, "Cd": 0.70,
        "P_chamber_bar": 30.0,
        "segments": [
            {"type": "pipe", "L_mm": 800.0, "D_mm": 12.0},
            {"type": "fitting", "D_mm": 12.0,
             "fitting_name": "Ball valve (fully open)", "K": 0.05},
            {"type": "pipe", "L_mm": 400.0, "D_mm": 12.0},
        ]
    },
    "Near-critical conditions (risk case)": {
        "T_tank_C": 28.0, "P_tank_bar": 60.0, "Cd": 0.62,
        "P_chamber_bar": 18.0,
        "segments": [
            {"type": "pipe", "L_mm": 2000.0, "D_mm": 8.0},
            {"type": "fitting", "D_mm": 8.0,
             "fitting_name": "Needle valve (fully open)", "K": 2.0},
            {"type": "fitting", "D_mm": 8.0,
             "fitting_name": "Elbow 90 deg (standard)", "K": 0.9},
        ]
    },
}

if "preset_applied" not in st.session_state:
    st.session_state.preset_applied = "Custom (current settings)"

# ── Caching ───────────────────────────────────────────────────────────────────
@st.cache_data
def _run_full_system(m_dot, T_tank, P_tank, seg_t, Cd, A, P_chamber, roughness):
    return evaluate_full_system(m_dot, T_tank, P_tank,
                                [dict(s) for s in seg_t], Cd, A, P_chamber, roughness)

@st.cache_data
def _run_feed_line(m_dot, T_tank, P_tank, seg_t, roughness):
    return evaluate_feed_line(m_dot, T_tank, P_tank,
                              [dict(s) for s in seg_t], roughness)

@st.cache_data
def _run_dyer(Cd, A, T_tank, P_inlet, P_chamber, rho_l_up, rho_l_down, rho_v_down):
    return dyer_mass_flow(Cd, A, T_tank, P_inlet, P_chamber,
                          rho_l_up, rho_l_down, rho_v_down)

def _to_tuple(segments):
    return tuple(tuple(sorted(s.items())) for s in segments)

# ── Sidebar ───────────────────────────────────────────────────────────────────
def render_sidebar():
    with st.sidebar:
        if st.button("Back to home"):
            st.session_state.page = "landing"; st.rerun()
        st.markdown("---")

        # Presets
        st.markdown('<div class="section-label">Quick start</div>',
                    unsafe_allow_html=True)
        chosen_preset = st.selectbox(
            "Load preset configuration",
            list(PRESETS.keys()),
            index=0,
            help="Load a pre-configured scenario as a starting point.")
        if chosen_preset != "Custom (current settings)" and chosen_preset != st.session_state.preset_applied:
            p = PRESETS[chosen_preset]
            st.session_state.segments = [dict(s) for s in p["segments"]]
            st.session_state.preset_applied = chosen_preset
            st.session_state["_preset_T"] = p["T_tank_C"]
            st.session_state["_preset_P"] = p["P_tank_bar"]
            st.session_state["_preset_Cd"] = p["Cd"]
            st.session_state["_preset_Pc"] = p["P_chamber_bar"]
            st.rerun()
        if chosen_preset == "Custom (current settings)":
            st.session_state.preset_applied = "Custom (current settings)"

        st.markdown("---")
        st.markdown('<div class="section-label">Tank conditions</div>',
                    unsafe_allow_html=True)
        T_default = st.session_state.get("_preset_T", 20.0)
        P_default = st.session_state.get("_preset_P", 55.0)
        T_C = st.number_input("Tank temperature (deg C)",
                               min_value=-10.0, max_value=35.0,
                               value=float(T_default), step=0.5)
        st.markdown('<div class="param-help">N2O temperature at tank exit. '
                    'Critical point: 36.4 deg C.</div>', unsafe_allow_html=True)
        P_bar = st.number_input("Tank pressure (bar)",
                                 min_value=5.0, max_value=71.0,
                                 value=float(P_default), step=0.5)
        st.markdown('<div class="param-help">Must exceed P_sat(T_tank) for '
                    'liquid at the outlet. '
                    'P_sat(20 deg C) = 51.4 bar.</div>', unsafe_allow_html=True)

        st.markdown('<div class="section-label">Injector</div>',
                    unsafe_allow_html=True)
        Cd_default = st.session_state.get("_preset_Cd", 0.65)
        Pc_default = st.session_state.get("_preset_Pc", 20.0)
        Cd = st.number_input("Discharge coefficient Cd",
                              min_value=0.3, max_value=0.99,
                              value=float(Cd_default), step=0.01)
        st.markdown('<div class="param-help">Sharp-edged ~0.61, '
                    'well-rounded ~0.82, typical injectors 0.60-0.75.</div>',
                    unsafe_allow_html=True)
        Pc_bar = st.number_input("Chamber pressure (bar)",
                                  min_value=1.0, max_value=65.0,
                                  value=float(Pc_default), step=0.5)
        st.markdown('<div class="param-help">Typically 15-30% below tank '
                    'pressure for combustion stability.</div>',
                    unsafe_allow_html=True)

        st.markdown('<div class="section-label">Advanced</div>',
                    unsafe_allow_html=True)
        with st.expander("Pipe roughness"):
            rough_um = st.number_input("Absolute roughness (um)",
                                        min_value=0.1, max_value=500.0,
                                        value=1.5, step=0.5)
            st.markdown('<div class="param-help">Smooth SS tubing ~1.5 um, '
                        'commercial steel ~46 um.</div>', unsafe_allow_html=True)

        T_tank = T_C + 273.15
        P_tank = P_bar * 1e5
        P_chamber = Pc_bar * 1e5
        roughness = rough_um * 1e-6

        # Subcooling margin badge
        Psat = P_sat(T_tank)
        margin = P_tank - Psat
        if margin > 0:
            st.markdown(f'<div class="badge-ok">&#10003; {margin/1e5:.2f} bar subcooling</div>',
                        unsafe_allow_html=True)
        else:
            st.markdown(f'<div class="badge-warn">&#9888; Below P_sat '
                        f'({Psat/1e5:.2f} bar)</div>', unsafe_allow_html=True)

        # Combustion stability check
        dP_inj = P_tank - P_chamber
        stab_pct = dP_inj / P_chamber * 100 if P_chamber > 0 else 0
        st.markdown("---")
        st.markdown('<div class="section-label">Combustion stability</div>',
                    unsafe_allow_html=True)
        if stab_pct >= 20:
            st.markdown(
                f'<div style="background:#0a1a12;color:#2d8a5a;border:1px solid #1a3028;'
                f'border-radius:3px;padding:.35rem .6rem;font-size:.75rem;'
                f'font-family:JetBrains Mono,monospace;line-height:1.4;">'
                f'&#10003; dP/Pc = {stab_pct:.0f}% &nbsp;(stable)</div>',
                unsafe_allow_html=True)
        elif stab_pct >= 15:
            st.markdown(
                f'<div style="background:#0a1525;color:#3d5a80;border:1px solid #152035;'
                f'border-radius:3px;padding:.35rem .6rem;font-size:.75rem;'
                f'font-family:JetBrains Mono,monospace;line-height:1.4;">'
                f'&#9888; dP/Pc = {stab_pct:.0f}% &nbsp;(marginal)</div>',
                unsafe_allow_html=True)
        else:
            st.markdown(
                f'<div style="background:#1a0d0d;color:#9a3030;border:1px solid #301515;'
                f'border-radius:3px;padding:.35rem .6rem;font-size:.75rem;'
                f'font-family:JetBrains Mono,monospace;line-height:1.4;">'
                f'&#9888; dP/Pc = {stab_pct:.0f}%<br>below 15% — instability risk</div>',
                unsafe_allow_html=True)
        st.markdown('<div class="param-help">Injector pressure drop / chamber '
                    'pressure. Min 15% recommended for combustion stability.</div>',
                    unsafe_allow_html=True)

        return T_tank, P_tank, Cd, P_chamber, roughness

# ── Segment list with insert-between ─────────────────────────────────────────
def render_segments(mode_key):
    from plotting import plot_line_profile
    st.markdown('<div class="section-label">Feed line geometry</div>',
                unsafe_allow_html=True)
    st.caption("Tank exit to injector inlet. Use Insert buttons to add between components.")

    def _add_pipe_at(pos):
        st.session_state.segments.insert(
            pos, {"type": "pipe", "L_mm": 500.0, "D_mm": 8.0})
        st.rerun()

    def _add_fitting_at(pos):
        st.session_state.segments.insert(
            pos, {"type": "fitting", "D_mm": 8.0,
                  "fitting_name": "Ball valve (fully open)", "K": 0.05})
        st.rerun()

    # Insert at the very beginning
    ci1, ci2, ci3 = st.columns([2, 2, 6])
    with ci1:
        if st.button("+ Pipe", key=f"{mode_key}_ins_pipe_0",
                     help="Insert pipe at start"):
            _add_pipe_at(0)
    with ci2:
        if st.button("+ Fitting", key=f"{mode_key}_ins_fit_0",
                     help="Insert fitting at start"):
            _add_fitting_at(0)

    for i, seg in enumerate(st.session_state.segments):
        # ── Segment row ──
        c1, c2, c3, c4, c_rm = st.columns([1.5, 2.0, 1.5, 1.8, 0.4])
        with c1:
            seg["type"] = st.selectbox(
                "Type", ["pipe", "fitting"],
                index=0 if seg["type"] == "pipe" else 1,
                key=f"{mode_key}_type_{i}",
                label_visibility="collapsed")
        if seg["type"] == "pipe":
            with c2:
                seg["L_mm"] = st.number_input(
                    "Length mm", value=float(seg.get("L_mm", 1000.0)),
                    min_value=1.0, step=50.0,
                    key=f"{mode_key}_L_{i}",
                    label_visibility="collapsed",
                    help="Pipe length in mm")
            with c3:
                seg["D_mm"] = st.number_input(
                    "ID mm", value=float(seg.get("D_mm", 8.0)),
                    min_value=0.5, step=0.5,
                    key=f"{mode_key}_D_{i}",
                    label_visibility="collapsed",
                    help="Inner diameter in mm")
            with c4:
                total_L = sum(s.get("L_mm", 0) for s in st.session_state.segments
                              if s["type"] == "pipe")
                st.markdown(
                    f'<div style="padding-top:.45rem;color:#3d4a5c;'
                    f'font-size:.75rem;">'
                    f'L={seg["L_mm"]:.0f} mm &nbsp;·&nbsp; '
                    f'ID={seg["D_mm"]:.1f} mm</div>',
                    unsafe_allow_html=True)
        else:
            with c2:
                fname = seg.get("fitting_name", "Ball valve (fully open)")
                options = list(FITTING_K.keys())
                idx = options.index(fname) if fname in options else 0
                chosen = st.selectbox(
                    "Fitting type", options, index=idx,
                    key=f"{mode_key}_fname_{i}",
                    label_visibility="collapsed")
                seg["fitting_name"] = chosen
                if FITTING_K[chosen] is not None:
                    seg["K"] = FITTING_K[chosen]
            with c3:
                seg["D_mm"] = st.number_input(
                    "ID mm", value=float(seg.get("D_mm", 8.0)),
                    min_value=0.5, step=0.5,
                    key=f"{mode_key}_D_{i}",
                    label_visibility="collapsed",
                    help="Inner diameter of adjacent pipe")
            with c4:
                if FITTING_K[chosen] is None:
                    seg["K"] = st.number_input(
                        "K", value=float(seg.get("K", 1.0)),
                        min_value=0.0, step=0.1,
                        key=f"{mode_key}_K_{i}",
                        label_visibility="collapsed",
                        help="Loss coefficient K")
                else:
                    kv = seg.get("K", 0)
                    st.markdown(
                        f'<div style="padding-top:.45rem;color:#3d4a5c;'
                        f'font-size:.75rem;">K = {kv:.2f}</div>',
                        unsafe_allow_html=True)
        with c_rm:
            if st.button("x", key=f"{mode_key}_rm_{i}", help="Remove segment"):
                st.session_state.segments.pop(i)
                st.rerun()

        # ── Insert after this segment ──
        ci1, ci2, ci3 = st.columns([2, 2, 6])
        with ci1:
            if st.button("+ Pipe", key=f"{mode_key}_ins_pipe_{i+1}",
                         help=f"Insert pipe after segment {i+1}"):
                _add_pipe_at(i + 1)
        with ci2:
            if st.button("+ Fitting", key=f"{mode_key}_ins_fit_{i+1}",
                         help=f"Insert fitting after segment {i+1}"):
                _add_fitting_at(i + 1)

    # Live schematic
    if st.session_state.segments:
        st.markdown('<div class="section-label" style="margin-top:.8rem">'
                    'Line schematic</div>', unsafe_allow_html=True)
        st.plotly_chart(
            plot_line_profile(st.session_state.segments),
            use_container_width=True,
            config={"staticPlot": True})

    model_segments = []
    for seg in st.session_state.segments:
        if seg["type"] == "pipe":
            model_segments.append({"type": "pipe",
                "L": seg["L_mm"] / 1000.0, "D": seg["D_mm"] / 1000.0})
        else:
            model_segments.append({"type": "fitting",
                "D": seg["D_mm"] / 1000.0,
                "K": float(seg.get("K", 0.05))})
    return model_segments

# ── Diagnostics ───────────────────────────────────────────────────────────────
def render_diagnostics(T_tank, P_tank, P_chamber, Cd, A_injector=None):
    Psat = P_sat(T_tank)
    suggestions = []
    pipe_segs = [s for s in st.session_state.segments if s["type"] == "pipe"]
    total_L = sum(s["L_mm"] for s in pipe_segs) / 1000.0
    min_D = min((s["D_mm"] for s in pipe_segs), default=8.0)
    high_k = [s for s in st.session_state.segments
               if s["type"] == "fitting" and s.get("K", 0) > 1.0]
    if P_tank <= Psat:
        suggestions.append(
            f"<b>Raise tank pressure</b> to at least {Psat/1e5+1:.1f} bar "
            f"(P_sat at {T_tank-273.15:.1f} deg C = {Psat/1e5:.2f} bar). "
            f"Current: {P_tank/1e5:.1f} bar.")
    if total_L > 1.5:
        suggestions.append(
            f"<b>Shorten the feed line</b>: {total_L:.1f} m total. "
            "Target under 1.5 m where possible.")
    if min_D < 10.0:
        suggestions.append(
            f"<b>Increase pipe inner diameter</b> (min {min_D:.0f} mm). "
            "Pressure drop scales as 1/D^4.")
    if high_k:
        suggestions.append(
            f"<b>Review {len(high_k)} high-K fitting(s)</b> (K > 1.0). "
            "Fully-open ball valves (K = 0.05) are strongly preferred.")
    if T_tank - 273.15 > 25.0:
        suggestions.append(
            f"<b>Pre-cool the oxidiser</b>: at {T_tank-273.15:.1f} deg C, "
            f"P_sat = {Psat/1e5:.2f} bar — margin is critically small.")
    st.markdown(
        '<div class="diag-box"><h4>Flashing detected — injector models not evaluated</h4>'
        "<ul>" + "".join(f"<li>{s}</li>" for s in suggestions) +
        ("" if suggestions else "<li>Review line geometry and tank pressure.</li>") +
        "</ul></div>", unsafe_allow_html=True)
    if A_injector and A_injector > 0:
        try:
            rho_sat = rho_liquid_sat(T_tank)
            dP = max(P_tank - P_chamber, 0.0)
            m_spi_ref = Cd * A_injector * math.sqrt(2 * rho_sat * dP)
            st.markdown(
                f'<div class="estimate-box"><h4>Reference — SPI prediction (not valid here)</h4>'
                f"<p>The SPI model (single-phase liquid, no flashing) would predict "
                f"<strong>{m_spi_ref*1000:.1f} g/s</strong>. "
                f"Because flashing has begun in the feed line, the real mass flow will be "
                f"<strong>substantially lower</strong> — documented cases show reductions "
                f"of 3 to 10 times below the SPI value under these conditions. "
                f"No reliable numerical estimate can be given without modelling two-phase "
                f"flow in the feed line itself (outside the current scope). "
                f"Fix the feed line flashing first, then re-evaluate.</p>"
                f"</div>", unsafe_allow_html=True)
        except Exception:
            pass

# ── Result cards ──────────────────────────────────────────────────────────────
def render_result_cards(result):
    if result["feed_line_result"]["flashing_detected"]:
        r1, r2, r3 = st.columns(3)
        with r1:
            st.markdown('<div class="result-card"><div class="label">Flashing in line</div>'
                        '<div class="value" style="color:#8a2020">YES</div></div>',
                        unsafe_allow_html=True)
        with r2:
            st.markdown('<div class="result-card"><div class="label">SPI sufficient</div>'
                        '<div class="value" style="color:#2a3040">N/A</div></div>',
                        unsafe_allow_html=True)
        with r3:
            st.markdown('<div class="result-card"><div class="label">Real mass flow</div>'
                        '<div class="value" style="color:#2a3040">N/A</div></div>',
                        unsafe_allow_html=True)
        return False
    spi_ok = result["spi_sufficient"]
    m_dot = result["m_dot_real"]
    ir = result["injector_result"]
    r1, r2, r3, r4 = st.columns(4)
    with r1:
        st.markdown('<div class="result-card"><div class="label">Flashing in line</div>'
                    '<div class="value" style="color:#2d6a4f">NO</div></div>',
                    unsafe_allow_html=True)
    with r2:
        col = "#2d6a4f" if spi_ok else "#7a6020"
        st.markdown(f'<div class="result-card"><div class="label">SPI sufficient</div>'
                    f'<div class="value" style="color:{col}">{"YES" if spi_ok else "NO"}</div>'
                    f'<div class="sub">{"SPI used" if spi_ok else "Dyer used"}</div></div>',
                    unsafe_allow_html=True)
    with r3:
        st.markdown(f'<div class="result-card"><div class="label">Real mass flow</div>'
                    f'<div class="value">{m_dot*1000:.1f}</div>'
                    f'<div class="sub">g/s</div></div>', unsafe_allow_html=True)
    with r4:
        if ir:
            over = (ir["m_dot_SPI"] - m_dot) / ir["m_dot_SPI"] * 100
            st.markdown(f'<div class="result-card"><div class="label">SPI over-prediction</div>'
                        f'<div class="value" style="color:#802020">{over:.1f}%</div>'
                        f'<div class="sub">SPI: {ir["m_dot_SPI"]*1000:.1f} g/s</div></div>',
                        unsafe_allow_html=True)
        else:
            st.markdown('<div class="result-card"><div class="label">SPI over-prediction</div>'
                        '<div class="value" style="color:#2d6a4f">0%</div>'
                        '<div class="sub">SPI valid here</div></div>',
                        unsafe_allow_html=True)
    if ir:
        st.caption(f"Dyer: kappa = {ir['kappa']:.3f}  |  "
                   f"exit vapour quality x = {ir['x_exit']:.3f}  |  "
                   f"HEM prediction: {ir['m_dot_HEM']*1000:.1f} g/s")
    return True

# ══════════════════════════════════════════════════════════════════════════════
# LANDING
# ══════════════════════════════════════════════════════════════════════════════
if st.session_state.page == "landing":
    st.markdown('<div class="hero-title">N2O Hybrid Rocket Injector Sizing</div>',
                unsafe_allow_html=True)
    st.markdown('<div class="hero-byline">Eduardo Costa &nbsp;·&nbsp; '
                'Instituto Superior Tecnico &nbsp;·&nbsp; '
                'Hybrid Rocket Propulsion</div>', unsafe_allow_html=True)
    st.markdown("""<div class="hero-desc">
        This tool predicts the real oxidiser mass flow through N2O hybrid rocket
        injectors, accounting for premature vaporisation (flashing) that causes the
        standard SPI model to over-predict flow by a factor of several times when
        the fluid approaches its saturation pressure along the feed path.
        The model chains three physics blocks: feed line pressure drop
        (Darcy-Weisbach), injector sufficiency check (SPI), and two-phase
        correction (Dyer/NHNE).
    </div>""", unsafe_allow_html=True)

    col1, col2 = st.columns(2, gap="large")
    with col1:
        st.markdown("""<div class="mode-card">
            <h3>Sizing mode</h3>
            <p>Known orifice geometry — predict the real mass flow the system
            will deliver, accounting for two-phase effects along the feed path.</p>
        </div>""", unsafe_allow_html=True)
        if st.button("Open — Sizing", use_container_width=True):
            st.session_state.page = "sizing"; st.rerun()
    with col2:
        st.markdown("""<div class="mode-card">
            <h3>Design mode</h3>
            <p>Target mass flow — find the orifice area required, with the
            Dyer correction applied so SPI under-sizing is avoided.</p>
        </div>""", unsafe_allow_html=True)
        if st.button("Open — Design", use_container_width=True):
            st.session_state.page = "design"; st.rerun()

    st.markdown("""<div class="assumption-box">
        <h4>Model assumptions and scope</h4>
        <ul>
            <li>Feed line treated as <strong>adiabatic</strong>: fluid temperature
                constant at T_tank along the line.</li>
            <li>Liquid density treated as <strong>incompressible</strong> along
                the line (function of T only).</li>
            <li><strong>Steady-state</strong> flow only — transient start-up
                effects not modelled.</li>
            <li>Tank temperature is a <strong>direct user input</strong>: thermal
                balance with the environment not modelled.</li>
            <li>Dyer model uses <strong>literature reference values</strong> for
                Cd and weighting, not team-calibrated data.</li>
            <li>When flashing is detected in the feed line, injector models are
                <strong>not evaluated</strong>. A conservative upper-bound
                estimate is provided instead.</li>
        </ul>
    </div>""", unsafe_allow_html=True)

    st.markdown("""<div class="disclaimer-box">
        <h4>Disclaimer</h4>
        <p>This tool is provided for <strong>predictive and educational purposes
        only</strong>. It is an academic project under active development and has
        not been independently validated against a comprehensive experimental
        dataset. Results should not be used as the sole basis for engineering
        decisions or hardware fabrication. The author accepts no responsibility
        for any damages, losses, or injuries arising from the use of this tool.
        Always verify critical design parameters with qualified engineers and
        appropriate experimental testing.</p>
    </div>""", unsafe_allow_html=True)

# ══════════════════════════════════════════════════════════════════════════════
# SIZING MODE
# ══════════════════════════════════════════════════════════════════════════════
if st.session_state.page == "sizing":
    T_tank, P_tank, Cd, P_chamber, roughness = render_sidebar()

    st.markdown("## Sizing mode")
    st.caption("Known orifice geometry — predicted real mass flow")

    with st.sidebar:
        st.markdown('<div class="section-label">Orifice geometry</div>',
                    unsafe_allow_html=True)
        N_holes = st.number_input("Number of orifices",
                                   min_value=1, max_value=500, value=1, step=1)
        d_mm = st.number_input("Orifice diameter (mm)",
                                min_value=0.1, max_value=30.0, value=2.19, step=0.01)
        st.markdown('<div class="param-help">Total area = N x pi x (d/2)^2. '
                    'Individual hole diameter.</div>', unsafe_allow_html=True)
        A_total = N_holes * math.pi * (d_mm / 2e3) ** 2
        st.markdown(f'<div class="badge-info">Total area: {A_total*1e6:.4f} mm2</div>',
                    unsafe_allow_html=True)

    m_dot_line = st.number_input(
        "Design mass flow for line evaluation (g/s)",
        min_value=10.0, max_value=5000.0, value=500.0, step=10.0,
        help="Used only to compute the feed line pressure drop. "
             "The injector result is independent of this value.")
    st.markdown('<div class="param-help">This value is used to calculate '
                'the velocity (and friction losses) in the feed line. '
                'It does not affect the injector model directly.</div>',
                unsafe_allow_html=True)

    model_segments = render_segments("siz")

    if not model_segments:
        st.info("Add at least one feed line segment using the buttons above.")
    else:
        try:
            seg_t = _to_tuple(model_segments)
            result = _run_full_system(m_dot_line / 1000.0, T_tank, P_tank,
                                      seg_t, Cd, A_total, P_chamber, roughness)
            st.markdown('<hr class="divider">', unsafe_allow_html=True)
            ok = render_result_cards(result)
            if not ok:
                render_diagnostics(T_tank, P_tank, P_chamber, Cd,
                                    A_injector=A_total)
            else:
                ir = result.get("injector_result")
                # Combustion stability check on actual injector dP
                P_inlet = result["P_injector_inlet"]
                actual_dP = P_inlet - P_chamber
                stab = actual_dP / P_chamber * 100
                if stab < 15:
                    st.warning(f"Injector dP/Pc = {stab:.1f}% — below the 15% "
                               "minimum recommended for combustion stability. "
                               "Consider increasing the pressure drop (smaller orifice "
                               "or higher tank pressure).")

                # Segment losses table
                with st.expander("Pressure drop by segment"):
                    st.plotly_chart(
                        plot_segment_losses(
                            result["feed_line_result"]["trace"],
                            model_segments),
                        use_container_width=True,
                        config={"scrollZoom": False})

                # Main diagrams
                st.markdown('<div class="section-label" '
                            'style="margin-top:1.4rem">Diagrams</div>',
                            unsafe_allow_html=True)
                pc1, pc2 = st.columns(2)
                with pc1:
                    st.plotly_chart(
                        plot_pressure_along_line(
                            result["feed_line_result"]["trace"],
                            P_tank, T_tank, model_segments),
                        use_container_width=True,
                        config={"scrollZoom": True})
                with pc2:
                    st.plotly_chart(
                        plot_PT_diagram(T_tank, P_tank,
                                        result["P_injector_inlet"], P_chamber),
                        use_container_width=True,
                        config={"scrollZoom": True})

                # Subcooling margin chart
                pc3, pc4 = st.columns(2)
                with pc3:
                    st.plotly_chart(
                        plot_subcooling_margin(
                            result["feed_line_result"]["trace"],
                            P_tank, T_tank, model_segments),
                        use_container_width=True,
                        config={"scrollZoom": True})
                with pc4:
                    st.plotly_chart(
                        plot_sensitivity(T_tank, P_tank, P_chamber,
                                         model_segments, Cd, A_total, roughness),
                        use_container_width=True,
                        config={"scrollZoom": True})

                # PDF export
                st.markdown('<div class="section-label" '
                            'style="margin-top:1rem">Export</div>',
                            unsafe_allow_html=True)
                result_copy = dict(result)
                result_copy["_T_tank"] = T_tank
                result_copy["_P_tank"] = P_tank
                result_copy["_P_chamber"] = P_chamber
                pdf_bytes = generate_pdf(
                    mode="sizing",
                    inputs={"T_tank_C": T_tank - 273.15,
                            "P_tank_bar": P_tank / 1e5,
                            "Cd": Cd,
                            "P_chamber_bar": P_chamber / 1e5,
                            "N_holes": N_holes,
                            "d_mm": d_mm,
                            "A_total": A_total,
                            "roughness_um": roughness * 1e6},
                    results=result_copy,
                    segments_ui=st.session_state.segments,
                    model_segments=model_segments)
                st.download_button(
                    label="Download PDF report",
                    data=pdf_bytes,
                    file_name="n2o_injector_sizing_report.pdf",
                    mime="application/pdf")
        except ValueError as e:
            st.error(f"Model error: {e}")

# ══════════════════════════════════════════════════════════════════════════════
# DESIGN MODE
# ══════════════════════════════════════════════════════════════════════════════
elif st.session_state.page == "design":
    T_tank, P_tank, Cd, P_chamber, roughness = render_sidebar()

    st.markdown("## Design mode")
    st.caption("Target mass flow — required orifice area with Dyer correction")

    with st.sidebar:
        st.markdown('<div class="section-label">Target</div>',
                    unsafe_allow_html=True)
        m_dot_target_gs = st.number_input(
            "Target oxidiser mass flow (g/s)",
            min_value=10.0, max_value=5000.0, value=300.0, step=10.0)
        st.markdown('<div class="param-help">Required flow at the design '
                    'O/F ratio and chamber pressure.</div>', unsafe_allow_html=True)
        N_holes = st.number_input("Number of orifices",
                                   min_value=1, max_value=500, value=4, step=1)
        st.markdown('<div class="param-help">Splits the total area across '
                    'N identical holes. More holes = smaller diameter each.</div>',
                    unsafe_allow_html=True)

    model_segments = render_segments("des")

    if not model_segments:
        st.info("Add at least one feed line segment using the buttons above.")
    else:
        try:
            m_dot_target = m_dot_target_gs / 1000.0
            seg_t = _to_tuple(model_segments)
            fl = _run_feed_line(m_dot_target, T_tank, P_tank, seg_t, roughness)
            P_inlet = fl["P_final"]

            st.markdown('<hr class="divider">', unsafe_allow_html=True)

            if fl["flashing_detected"]:
                r1, r2, r3 = st.columns(3)
                with r1:
                    st.markdown('<div class="result-card">'
                                '<div class="label">Flashing in line</div>'
                                '<div class="value" style="color:#8a2020">YES</div>'
                                '</div>', unsafe_allow_html=True)
                with r2:
                    st.markdown('<div class="result-card">'
                                '<div class="label">SPI area</div>'
                                '<div class="value" style="color:#2a3040">N/A</div>'
                                '</div>', unsafe_allow_html=True)
                with r3:
                    st.markdown('<div class="result-card">'
                                '<div class="label">Dyer area</div>'
                                '<div class="value" style="color:#2a3040">N/A</div>'
                                '</div>', unsafe_allow_html=True)
                render_diagnostics(T_tank, P_tank, P_chamber, Cd)
            else:
                dP_inj = P_inlet - P_chamber
                rho_l_up = rho_liquid_sat(T_tank)
                A_spi = orifice_area_from_target_flow(
                    m_dot_target, Cd, rho_l_up, dP_inj)
                A_iter = A_spi
                T_down = T_sat(P_chamber)
                rho_l_down = rho_liquid_sat(T_down)
                rho_v_down = M_N2O / nu_vapor_sat(T_down)
                for _ in range(40):
                    dr = _run_dyer(Cd, A_iter, T_tank, P_inlet, P_chamber,
                                   rho_l_up, rho_l_down, rho_v_down)
                    if abs(dr["m_dot_Dyer"] - m_dot_target) / m_dot_target < 1e-7:
                        break
                    A_iter *= m_dot_target / dr["m_dot_Dyer"]
                A_dyer = A_iter
                d_spi  = math.sqrt(4 * A_spi  / (N_holes * math.pi)) * 1000
                d_dyer = math.sqrt(4 * A_dyer / (N_holes * math.pi)) * 1000
                pct = (A_dyer - A_spi) / A_spi * 100

                r1, r2, r3, r4 = st.columns(4)
                with r1:
                    st.markdown(f'<div class="result-card">'
                                f'<div class="label">SPI area</div>'
                                f'<div class="value">{A_spi*1e6:.3f}</div>'
                                f'<div class="sub">mm2 total</div></div>',
                                unsafe_allow_html=True)
                with r2:
                    st.markdown(f'<div class="result-card">'
                                f'<div class="label">SPI hole diameter</div>'
                                f'<div class="value">{d_spi:.3f}</div>'
                                f'<div class="sub">mm x {N_holes}</div></div>',
                                unsafe_allow_html=True)
                with r3:
                    st.markdown(f'<div class="result-card">'
                                f'<div class="label">Dyer area</div>'
                                f'<div class="value" style="color:#2d6a4f">'
                                f'{A_dyer*1e6:.3f}</div>'
                                f'<div class="sub">mm2 total</div></div>',
                                unsafe_allow_html=True)
                with r4:
                    st.markdown(f'<div class="result-card">'
                                f'<div class="label">Dyer hole diameter</div>'
                                f'<div class="value" style="color:#2d6a4f">'
                                f'{d_dyer:.3f}</div>'
                                f'<div class="sub">mm x {N_holes}</div></div>',
                                unsafe_allow_html=True)

                st.caption(
                    f"Dyer area is {pct:.1f}% larger than SPI — the additional area "
                    "compensates for two-phase flow reduction so the system delivers "
                    "the target mass flow.")

                # Combustion stability check
                actual_dP = P_inlet - P_chamber
                stab = actual_dP / P_chamber * 100
                if stab < 15:
                    st.warning(f"Injector dP/Pc = {stab:.1f}% — below the 15% "
                               "minimum recommended for combustion stability.")

                st.markdown('<div class="section-label" '
                            'style="margin-top:1.4rem">Diagrams</div>',
                            unsafe_allow_html=True)
                pc1, pc2 = st.columns(2)
                with pc1:
                    st.plotly_chart(
                        plot_model_comparison(
                            dr["m_dot_SPI"], dr["m_dot_HEM"],
                            dr["m_dot_Dyer"], m_dot_target),
                        use_container_width=True,
                        config={"scrollZoom": False})
                with pc2:
                    st.plotly_chart(
                        plot_PT_diagram(T_tank, P_tank, P_inlet, P_chamber),
                        use_container_width=True,
                        config={"scrollZoom": True})

                # Sensitivity for the Dyer area
                pc3, pc4 = st.columns(2)
                with pc3:
                    st.plotly_chart(
                        plot_sensitivity(T_tank, P_tank, P_chamber,
                                         model_segments, Cd, A_dyer, roughness),
                        use_container_width=True,
                        config={"scrollZoom": True})
                with pc4:
                    result_check = _run_full_system(
                        m_dot_target, T_tank, P_tank, seg_t,
                        Cd, A_dyer, P_chamber, roughness)
                    st.plotly_chart(
                        plot_subcooling_margin(
                            result_check["feed_line_result"]["trace"],
                            P_tank, T_tank, model_segments),
                        use_container_width=True,
                        config={"scrollZoom": True})

                # PDF export
                st.markdown('<div class="section-label" '
                            'style="margin-top:1rem">Export</div>',
                            unsafe_allow_html=True)
                export_results = {
                    "feed_line_result": result_check["feed_line_result"],
                    "P_injector_inlet": P_inlet,
                    "_T_tank": T_tank, "_P_tank": P_tank, "_P_chamber": P_chamber,
                    "A_spi": A_spi, "A_dyer": A_dyer,
                    "d_spi": d_spi, "d_dyer": d_dyer, "pct": pct,
                    "m_spi": dr["m_dot_SPI"], "m_hem": dr["m_dot_HEM"],
                    "m_dyer": dr["m_dot_Dyer"], "m_target": m_dot_target,
                }
                pdf_bytes = generate_pdf(
                    mode="design",
                    inputs={"T_tank_C": T_tank - 273.15,
                            "P_tank_bar": P_tank / 1e5,
                            "Cd": Cd,
                            "P_chamber_bar": P_chamber / 1e5,
                            "m_dot_target_gs": m_dot_target_gs,
                            "N_holes": N_holes,
                            "roughness_um": roughness * 1e6},
                    results=export_results,
                    segments_ui=st.session_state.segments,
                    model_segments=model_segments)
                st.download_button(
                    label="Download PDF report",
                    data=pdf_bytes,
                    file_name="n2o_injector_design_report.pdf",
                    mime="application/pdf")

        except ValueError as e:
            st.error(f"Model error: {e}")
