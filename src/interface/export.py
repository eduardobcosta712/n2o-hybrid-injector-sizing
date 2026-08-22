"""
export.py

PDF report generation for the N2O injector sizing tool.
Uses ReportLab for layout and Matplotlib for chart rendering
(avoids the Chrome dependency of Plotly/kaleido).

The report includes: input summary, key results, pressure-along-line
chart, and P-T diagram -- one page, ready to attach to a project report.
"""

import io
import math
import datetime
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from reportlab.lib.units import mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                 Table, TableStyle, Image as RLImage,
                                 HRFlowable)
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_RIGHT

# Colours matching the app's dark theme, mapped to light paper
C_BLUE   = colors.HexColor("#2c5f8a")
C_GREEN  = colors.HexColor("#2d6a4f")
C_RED    = colors.HexColor("#7a2020")
C_DARK   = colors.HexColor("#1a2030")
C_LIGHT  = colors.HexColor("#f0f4f8")
C_MID    = colors.HexColor("#d0d8e4")
C_ACCENT = colors.HexColor("#4a6fa5")
C_WARN   = colors.HexColor("#b05020")


def _fig_to_image(fig, width_mm=85, height_mm=55):
    """Render a Matplotlib figure to a ReportLab Image object."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=180, bbox_inches="tight",
                facecolor=fig.get_facecolor())
    buf.seek(0)
    plt.close(fig)
    return RLImage(buf, width=width_mm * mm, height=height_mm * mm)


def _make_pressure_chart(trace, P_tank, T_tank, segments):
    """Matplotlib version of the pressure-along-line chart for PDF export."""
    from n2o_properties import P_sat
    fig, ax = plt.subplots(figsize=(4.5, 2.8))
    fig.patch.set_facecolor("#f8fafc")
    ax.set_facecolor("#f0f4f8")

    x_positions = [0.0]
    x_cur = 0.0
    for seg in segments:
        if seg["type"] == "pipe":
            x_cur += seg["L"]
        x_positions.append(x_cur)

    P_bar = [P_tank / 1e5] + [s["pressure_after_Pa"] / 1e5 for s in trace]
    P_sat_val = P_sat(T_tank) / 1e5

    ax.plot(x_positions, P_bar, marker="o", color="#2c5f8a",
            linewidth=1.8, markersize=5, zorder=3)
    ax.axhline(P_sat_val, color="#7a2020", linestyle="--",
               linewidth=1.2, label=f"P_sat = {P_sat_val:.1f} bar")
    ax.fill_between(x_positions,
                    [P_sat_val] * len(x_positions), P_bar,
                    where=[p > P_sat_val for p in P_bar],
                    color="#2c5f8a", alpha=0.08)
    ax.set_xlabel("Position along line (m)", fontsize=8)
    ax.set_ylabel("Pressure (bar)", fontsize=8)
    ax.set_title("Pressure along feed line", fontsize=9, fontweight="bold")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    ax.tick_params(labelsize=7)
    fig.tight_layout(pad=0.5)
    return fig


def _make_pt_chart(T_tank, P_tank, P_inlet, P_chamber):
    """Matplotlib version of the P-T diagram for PDF export."""
    from n2o_properties import P_sat, T_sat, T_MIN, T_MAX
    fig, ax = plt.subplots(figsize=(4.5, 2.8))
    fig.patch.set_facecolor("#f8fafc")
    ax.set_facecolor("#f0f4f8")

    T_curve = np.linspace(T_MIN, T_MAX, 200)
    P_curve = [P_sat(T) / 1e5 for T in T_curve]
    T_c = T_curve - 273.15

    ax.plot(T_c, P_curve, color="#1a2030", linewidth=1.8, label="Saturation curve")
    ax.fill_between(T_c, P_curve, 80, color="#2c5f8a", alpha=0.07)
    ax.fill_between(T_c, 0, P_curve, color="#c06020", alpha=0.06)

    T_tc = T_tank - 273.15
    T_chc = T_sat(P_chamber) - 273.15
    ax.scatter([T_tc], [P_tank / 1e5], color="#2c5f8a", s=50,
               zorder=5, label=f"Tank ({P_tank/1e5:.1f} bar)", marker="s")
    ax.scatter([T_tc], [P_inlet / 1e5], color="#2d6a4f", s=50,
               zorder=5, label=f"Injector inlet ({P_inlet/1e5:.1f} bar)", marker="^")
    ax.scatter([T_chc], [P_chamber / 1e5], color="#7a2020", s=50,
               zorder=5, label=f"Chamber ({P_chamber/1e5:.1f} bar)", marker="v")

    ax.set_ylim(0, 80)
    ax.set_xlabel("Temperature (deg C)", fontsize=8)
    ax.set_ylabel("Pressure (bar)", fontsize=8)
    ax.set_title("P-T diagram", fontsize=9, fontweight="bold")
    ax.legend(fontsize=6)
    ax.grid(alpha=0.3)
    ax.tick_params(labelsize=7)
    fig.tight_layout(pad=0.5)
    return fig


def _make_comparison_chart(m_spi, m_hem, m_dyer, m_target):
    """Matplotlib bar chart for model comparison in Design mode PDF."""
    fig, ax = plt.subplots(figsize=(4.5, 2.8))
    fig.patch.set_facecolor("#f8fafc")
    ax.set_facecolor("#f0f4f8")

    labels = ["SPI", "HEM", "Dyer", "Target"]
    vals = [m_spi * 1000, m_hem * 1000, m_dyer * 1000, m_target * 1000]
    clrs = ["#4a6fa5", "#2a4a3a", "#2d6a4f", "#888888"]
    bars = ax.bar(labels, vals, color=clrs, edgecolor="white",
                  linewidth=0.5, width=0.5)
    ax.axhline(m_target * 1000, color="#888888", linestyle=":",
               linewidth=1.2)
    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v + max(vals) * 0.02,
                f"{v:.0f}", ha="center", va="bottom", fontsize=7)
    ax.set_ylabel("Mass flow (g/s)", fontsize=8)
    ax.set_title("Model comparison", fontsize=9, fontweight="bold")
    ax.grid(alpha=0.3, axis="y")
    ax.tick_params(labelsize=7)
    fig.tight_layout(pad=0.5)
    return fig


def generate_pdf(mode, inputs, results, segments_ui, model_segments):
    """
    Generate a one-page PDF summary report.

    Parameters
    ----------
    mode : str
        "sizing" or "design"
    inputs : dict
        All user inputs (T_tank_C, P_tank_bar, Cd, P_chamber_bar, etc.)
    results : dict
        Full result dict from evaluate_full_system or design mode calculation
    segments_ui : list
        Raw UI segment list (for schematic labels)
    model_segments : list
        SI segment list (for chart generation)

    Returns
    -------
    bytes
        PDF file content as bytes, ready for st.download_button.
    """
    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4,
                            leftMargin=18*mm, rightMargin=18*mm,
                            topMargin=16*mm, bottomMargin=16*mm)

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle("title", fontSize=14, textColor=C_DARK,
                                  fontName="Helvetica-Bold", spaceAfter=2)
    sub_style = ParagraphStyle("sub", fontSize=8, textColor=C_ACCENT,
                                fontName="Helvetica", spaceAfter=8)
    section_style = ParagraphStyle("section", fontSize=9, textColor=C_BLUE,
                                    fontName="Helvetica-Bold", spaceBefore=8,
                                    spaceAfter=4)
    body_style = ParagraphStyle("body", fontSize=8, textColor=C_DARK,
                                 fontName="Helvetica", spaceAfter=4,
                                 leading=12)
    warn_style = ParagraphStyle("warn", fontSize=7, textColor=C_WARN,
                                 fontName="Helvetica-Oblique", spaceAfter=6,
                                 leading=10)

    story = []
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    mode_label = "Sizing mode" if mode == "sizing" else "Design mode"

    # ── Header ───────────────────────────────────────────────────────────────
    story.append(Paragraph("N2O Hybrid Rocket Injector Sizing", title_style))
    story.append(Paragraph(f"{mode_label} &nbsp;·&nbsp; Generated {now}",
                           sub_style))
    story.append(Paragraph("Eduardo Costa &nbsp;·&nbsp; Instituto Superior Tecnico",
                           sub_style))
    story.append(HRFlowable(width="100%", thickness=1,
                             color=C_ACCENT, spaceAfter=6))

    # ── Disclaimer ───────────────────────────────────────────────────────────
    story.append(Paragraph(
        "<b>Disclaimer:</b> This report is generated by an academic tool under "
        "active development. Results are for predictive purposes only and must "
        "not be used as the sole basis for engineering decisions or hardware "
        "fabrication. The author accepts no responsibility for any damages "
        "arising from the use of this tool.", warn_style))
    story.append(HRFlowable(width="100%", thickness=0.5,
                             color=C_MID, spaceAfter=4))

    # ── Input summary table ───────────────────────────────────────────────────
    # Page usable width = 210 - 36 = 174mm. Split: label 52, value 35, label 52, value 35
    story.append(Paragraph("Input Parameters", section_style))
    input_data = [
        ["Parameter", "Value", "Parameter", "Value"],
        ["Tank temperature", f"{inputs['T_tank_C']:.1f} deg C",
         "Chamber pressure", f"{inputs['P_chamber_bar']:.1f} bar"],
        ["Tank pressure", f"{inputs['P_tank_bar']:.1f} bar",
         "Discharge coeff Cd", f"{inputs['Cd']:.2f}"],
    ]
    if mode == "sizing":
        input_data.append([
            "No. orifices", str(inputs.get("N_holes", "-")),
            "Orifice diameter", f"{inputs.get('d_mm', 0):.2f} mm"])
        input_data.append([
            "Total area", f"{inputs.get('A_total', 0)*1e6:.4f} mm2",
            "Pipe roughness", f"{inputs.get('roughness_um', 1.5):.1f} um"])
    else:
        input_data.append([
            "Target flow", f"{inputs.get('m_dot_target_gs', 0):.0f} g/s",
            "No. orifices", str(inputs.get("N_holes", "-"))])
        input_data.append([
            "Pipe roughness", f"{inputs.get('roughness_um', 1.5):.1f} um",
            "", ""])

    t = Table(input_data, colWidths=[52*mm, 35*mm, 52*mm, 35*mm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), C_DARK),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 7.5),
        ("BACKGROUND", (0, 1), (-1, -1), C_LIGHT),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [C_LIGHT, colors.white]),
        ("GRID", (0, 0), (-1, -1), 0.5, C_MID),
        ("FONTNAME", (0, 1), (0, -1), "Helvetica-Bold"),
        ("FONTNAME", (2, 1), (2, -1), "Helvetica-Bold"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("WORDWRAP", (0, 0), (-1, -1), True),
    ]))
    story.append(t)

    # ── Feed line segments table ──────────────────────────────────────────────
    story.append(Paragraph("Feed Line Geometry", section_style))
    seg_data = [["#", "Type", "Details", "Loss (approx.)"]]
    for i, seg in enumerate(segments_ui):
        if seg["type"] == "pipe":
            detail = f"L = {seg['L_mm']:.0f} mm, ID = {seg['D_mm']:.1f} mm"
            loss = "friction (Darcy-Weisbach)"
        else:
            detail = f"{seg.get('fitting_name','fitting')}, ID = {seg['D_mm']:.1f} mm"
            loss = f"K = {seg.get('K', 0):.2f}"
        seg_data.append([str(i + 1), seg["type"].capitalize(), detail, loss])

    ts = Table(seg_data, colWidths=[8*mm, 16*mm, 80*mm, 48*mm])
    ts.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), C_DARK),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 7.5),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [C_LIGHT, colors.white]),
        ("GRID", (0, 0), (-1, -1), 0.4, C_MID),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 5),
        ("WORDWRAP", (0, 0), (-1, -1), True),
    ]))
    story.append(ts)

    # ── Results ───────────────────────────────────────────────────────────────
    story.append(Paragraph("Results", section_style))
    fl = results.get("feed_line_result", {})
    flashing = fl.get("flashing_detected", False)

    if mode == "sizing":
        res_data = [["Metric", "Value"]]
        res_data.append(["Flashing in feed line",
                          "YES — injector models not evaluated" if flashing else "NO"])
        if not flashing:
            res_data.append(["Injector inlet pressure",
                              f"{results.get('P_injector_inlet', 0)/1e5:.2f} bar"])
            res_data.append(["Model used",
                              "SPI" if results.get("spi_sufficient") else "Dyer (two-phase)"])
            res_data.append(["Real mass flow",
                              f"{results.get('m_dot_real', 0)*1000:.1f} g/s"])
            ir = results.get("injector_result")
            if ir:
                over = (ir["m_dot_SPI"] - results["m_dot_real"]) / ir["m_dot_SPI"] * 100
                res_data.append(["SPI prediction (reference)",
                                  f"{ir['m_dot_SPI']*1000:.1f} g/s"])
                res_data.append(["SPI over-prediction",
                                  f"{over:.1f}%"])
                res_data.append(["Dyer kappa",
                                  f"{ir['kappa']:.3f}"])
                res_data.append(["Exit vapour quality x",
                                  f"{ir['x_exit']:.3f}"])
    else:
        res_data = [["Metric", "Value"]]
        res_data.append(["Flashing in feed line",
                          "YES — design not valid" if flashing else "NO"])
        if not flashing:
            res_data.append(["SPI total area",
                              f"{results.get('A_spi', 0)*1e6:.4f} mm2"])
            res_data.append(["SPI hole diameter",
                              f"{results.get('d_spi', 0):.4f} mm"])
            res_data.append(["Dyer total area (recommended)",
                              f"{results.get('A_dyer', 0)*1e6:.4f} mm2"])
            res_data.append(["Dyer hole diameter (recommended)",
                              f"{results.get('d_dyer', 0):.4f} mm"])
            res_data.append(["Dyer vs SPI area increase",
                              f"{results.get('pct', 0):.1f}%"])

    tr = Table(res_data, colWidths=[90*mm, 58*mm])
    tr.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), C_DARK),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE", (0, 0), (-1, -1), 8),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [C_LIGHT, colors.white]),
        ("GRID", (0, 0), (-1, -1), 0.4, C_MID),
        ("FONTNAME", (0, 1), (0, -1), "Helvetica-Bold"),
        ("TOPPADDING", (0, 0), (-1, -1), 3),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
    ]))
    story.append(tr)

    # ── Charts ────────────────────────────────────────────────────────────────
    if not flashing and fl.get("trace"):
        story.append(Paragraph("Diagrams", section_style))
        fig1 = _make_pressure_chart(
            fl["trace"], results.get("_P_tank", 55e5),
            results.get("_T_tank", 293.15), model_segments)
        fig2 = _make_pt_chart(
            results.get("_T_tank", 293.15),
            results.get("_P_tank", 55e5),
            results.get("P_injector_inlet", 50e5),
            results.get("_P_chamber", 20e5))

        img1 = _fig_to_image(fig1, width_mm=87, height_mm=54)
        img2 = _fig_to_image(fig2, width_mm=87, height_mm=54)
        chart_table = Table([[img1, img2]], colWidths=[90*mm, 90*mm])
        chart_table.setStyle(TableStyle([
            ("ALIGN", (0, 0), (-1, -1), "CENTER"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ]))
        story.append(chart_table)

        if mode == "design" and "m_spi" in results:
            fig3 = _make_comparison_chart(
                results["m_spi"], results["m_hem"],
                results["m_dyer"], results["m_target"])
            img3 = _fig_to_image(fig3, width_mm=87, height_mm=54)
            story.append(Spacer(1, 4*mm))
            story.append(img3)

    # ── Footer ────────────────────────────────────────────────────────────────
    story.append(Spacer(1, 4*mm))
    story.append(HRFlowable(width="100%", thickness=0.5, color=C_MID))
    story.append(Paragraph(
        "Generated by N2O Hybrid Rocket Injector Sizing Tool — Eduardo Costa, "
        "Instituto Superior Tecnico. "
        "This document is provided for informational purposes only.",
        ParagraphStyle("footer", fontSize=6.5, textColor=colors.grey,
                        fontName="Helvetica-Oblique")))

    doc.build(story)
    buf.seek(0)
    return buf.read()
