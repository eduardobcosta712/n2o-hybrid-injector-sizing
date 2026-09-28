"""
explore_supercharge_correction.py

*** EXPLORATORIO -- NAO altera injector_two_phase.py nem full_system.py ***

Testa a hipotese de que uma correcao dependente de supercharge no parametro
kappa de Dyer,

    kappa' = kappa * (supercharge / supercharge_ref) ** beta

reduz o MAPE observado a baixo supercharge (docs/future_work.md, Priority 1,
ponto 1), SEM degradar a banda ja validada (dP = 8-14 bar, alto supercharge).

COMO CORRER (no teu repo local, com CoolProp instalado):
  1. Aplica o PATCH_injector_two_phase.py: cola a funcao
     dyer_mass_flow_corrected() em src/model/injector_two_phase.py, logo a
     seguir a dyer_mass_flow() existente (nao mexe em mais nada).
  2. Copia este ficheiro para validation/explore_supercharge_correction.py
  3. python validation/explore_supercharge_correction.py
     (corre de qualquer diretoria, tal como o waxman_2013_validation.py)

Reutiliza directamente as funcoes de carregamento e calibracao de
validation/waxman_2013_validation.py (load_multiseries, curve_state,
spi_cd_samples) -- nao ha logica duplicada para o parsing dos CSVs
digitalizados.

Metodologia (espelha run_part_b() de waxman_2013_validation.py):
  1. Cd calibrado UMA vez, pooled sobre a janela de fase unica de cada
     curva (30 psi <= dP <= supercharge) -- kappa nao entra aqui, e por
     isso a calibracao de Cd nao e afetada pela correcao.
  2. Para cada (beta, supercharge_ref) numa grelha, recalcula-se
     m_dot_Dyer_corrected em TODOS os pontos de duas fases (dP > supercharge)
     e mede-se o MAPE global e por banda de supercharge/dP.
  3. Restricao dura: a banda "sagrada" (Parte A original, os 4 pontos
     Nino&Razavi, dP 8-14 bar, injector 2 geometry, Cd=0.65) nao pode
     degradar significativamente -- IMPORTANTE: usa sempre CD_PART_A=0.65,
     NAO o Cd pooled da Fig.13 (injector 3) -- sao geometrias diferentes.
"""

import math
import os
import sys

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO_ROOT, "src", "model"))

from n2o_properties import P_sat, T_sat, rho_liquid_sat, nu_vapor_sat, M_N2O
from injector_two_phase import dyer_mass_flow, dyer_mass_flow_corrected

# Reutiliza o loader e a calibracao ja escritos e testados em
# waxman_2013_validation.py -- evita duplicar logica de parsing.
sys.path.insert(0, _REPO_ROOT)
from validation.waxman_2013_validation import (
    load_multiseries, curve_state, spi_cd_samples, A_INJ3, MIN_DP_PSI,
    PSI_PA, BAR_PER_PSI, pct_error, summarize,
)

# --- Parte A original (Nino & Razavi 2019, injector 2 geometry) ------------
# Geometria DIFERENTE do injector 3 usado na Fig.13 -- Cd proprio, nao
# o cd_pooled calibrado mais abaixo.
T1_A = 280.0
P1_A = 4.36e6
D_A = 0.0015
A_A = math.pi * (D_A / 2.0) ** 2
CD_PART_A = 0.65   # valor generico ja assumido pelo repo para o injector 2
CD_PART_A_FIG15 = 0.681   # Cd MEDIDO para o injector 2 (Waxman Fig.15, digitised) --
                          # ver waxman_2013_validation.py, print_part_d() /
                          # a sensibilidade ja feita la: MAPE baseline sobe de
                          # 2.76% para 4.41% so por trocar 0.65 -> 0.681 (sem
                          # nenhuma correcao de kappa envolvida). Usado abaixo
                          # para isolar quanto da degradacao com a correcao
                          # "agressiva" e de facto efeito de kappa, e quanto
                          # e so o Cd assumido ja nao ser o correto.
CASES_A = [
    ("Pre-critical",  0.84, 44.0),
    ("Critical",      0.98, 46.5),
    ("Post-critical1", 1.09, 47.5),
    ("Post-critical2", 1.37, 48.0),
]


def _rho_v(T):
    return M_N2O / nu_vapor_sat(T)


def two_phase_points(curves):
    """Todos os pontos de duas fases (dP > supercharge) do injector 3,
    com o estado a montante -- mesma selecao que dy_rows em run_part_b()."""
    pts = []
    for c in curves:
        T, P1_pa = curve_state(c)
        rho_l_up = rho_liquid_sat(T)
        for dP, m_exp in zip(c["dP_psi"], c["y"]):
            if dP <= c["super_psi"] or dP < MIN_DP_PSI:
                continue
            P2 = P1_pa - dP * PSI_PA
            if P2 <= 1e4:
                continue
            pts.append({"super_psi": c["super_psi"], "dP_psi": dP,
                        "dP_bar": dP * BAR_PER_PSI,
                        "m_exp": m_exp, "T": T, "P1": P1_pa, "P2": P2,
                        "rho_l_up": rho_l_up})
    return pts


def predict_baseline(pts, cd):
    """Dyer original (beta=0), sem linha a montante -- negligivel, tal
    como assumido no script original (LINE ~= 5 cm de 25.4 mm ID)."""
    out = []
    for p in pts:
        T_d = T_sat(p["P2"])
        rho_l_d = rho_liquid_sat(T_d)
        rho_v_d = _rho_v(T_d)
        r = dyer_mass_flow(cd, A_INJ3, p["T"], p["P1"], p["P2"],
                           p["rho_l_up"], rho_l_d, rho_v_d)
        out.append({**p, "m_model": r["m_dot_Dyer"],
                    "err": pct_error(r["m_dot_Dyer"], p["m_exp"])})
    return out


def predict_corrected(pts, cd, beta, supercharge_ref_pa):
    out = []
    for p in pts:
        T_d = T_sat(p["P2"])
        rho_l_d = rho_liquid_sat(T_d)
        rho_v_d = _rho_v(T_d)
        r = dyer_mass_flow_corrected(cd, A_INJ3, p["T"], p["P1"], p["P2"],
                                     p["rho_l_up"], rho_l_d, rho_v_d,
                                     beta=beta, supercharge_ref_Pa=supercharge_ref_pa)
        m = r["m_dot_Dyer_corrected"]
        out.append({**p, "m_model": m, "err": pct_error(m, p["m_exp"]),
                    "kappa": r["kappa"], "kappa_corr": r["kappa_corrected"]})
    return out


def fmt(s):
    return f"n={s['n']:3d}  mean={s['mean']:+6.2f}%  MAPE={s['mape']:5.2f}%  max|err|={s['max']:5.1f}%"


def run_part_a_corrected(beta, supercharge_ref_pa, cd=None):
    """
    Verifica a banda sagrada (8-14 bar, Nino & Razavi) com a correcao
    aplicada.

    cd : float ou None
        Se None (default), usa CD_PART_A (0.65 -- o valor generico ja
        assumido pelo repo). Passa CD_PART_A_FIG15 (0.681, o valor
        MEDIDO para esta geometria) para testar se a degradacao vista
        com a correcao "agressiva" e efeito genuino da correcao de
        kappa, ou (parcial/totalmente) um artefacto do Cd assumido ja
        estar errado antes de qualquer correcao entrar em jogo.
    """
    if cd is None:
        cd = CD_PART_A
    rows = []
    rho_l_up = rho_liquid_sat(T1_A)
    for label, dP_MPa, m_exp_gs in CASES_A:
        P2 = P1_A - dP_MPa * 1e6
        T_d = T_sat(P2)
        rho_l_d = rho_liquid_sat(T_d)
        rho_v_d = _rho_v(T_d)
        r = dyer_mass_flow_corrected(cd, A_A, T1_A, P1_A, P2,
                                     rho_l_up, rho_l_d, rho_v_d,
                                     beta=beta, supercharge_ref_Pa=supercharge_ref_pa)
        m = r["m_dot_Dyer_corrected"]
        m_exp_kgs = m_exp_gs * 1e-3
        rows.append({"label": label, "err": pct_error(m, m_exp_kgs),
                     "m_model": m, "m_exp": m_exp_kgs})
    return rows


def main():
    print("=" * 78)
    print("EXPLORACAO: correcao GATED de kappa dependente do supercharge")
    print("(CoolProp real -- confirmar que 'import CoolProp' funciona antes)")
    print("=" * 78)

    curves, _n_dropped = load_multiseries("waxman_fig13_mdot_vs_dP_by_supercharge.csv")
    samples = [v for c in curves for v in spi_cd_samples(c)]
    cd_pooled = sum(samples) / len(samples)
    print(f"\nCd pooled (janela monofasica, injector 3, n={len(samples)}): {cd_pooled:.4f}")

    pts = two_phase_points(curves)
    print(f"Pontos de duas fases (dP > supercharge): {len(pts)}")

    supercharge_A_psi = (P1_A - P_sat(T1_A)) / PSI_PA
    print(f"\nATENCAO: o supercharge dos 4 pontos da Parte A e "
          f"{supercharge_A_psi:.1f} psi -- QUALQUER supercharge_ref >= "
          f"{supercharge_A_psi:.0f} psi faz a correcao atuar TAMBEM na "
          f"Parte A (nao fica protegida so por a correcao ser 'gated').")

    # --- Baseline (beta=0, Dyer original) ---
    base = predict_baseline(pts, cd_pooled)
    print(f"\nBASELINE (Dyer original) -- deve bater com waxman_2013_results.md:")
    print("  Global:", fmt(summarize([r["err"] for r in base])))
    for lo, hi, name in [(0, 100, "supercharge < 100 psi"),
                         (100, 200, "100 <= supercharge < 200 psi"),
                         (200, 1000, "supercharge >= 200 psi")]:
        b = [r["err"] for r in base if lo <= r["super_psi"] < hi]
        print(f"  {name:<32}", fmt(summarize(b)))

    part_a_base = run_part_a_corrected(beta=0.0, supercharge_ref_pa=200 * PSI_PA)
    print("  Parte A (4 pontos, 8-14 bar, deve dar ~3.51%):",
          fmt(summarize([r["err"] for r in part_a_base])))

    # --- Grid search (GATED: factor=1 fora da zona corrigida) ---
    print(f"\n{'='*78}\nGRID SEARCH (GATED): beta x supercharge_ref\n{'='*78}")
    betas = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0]
    # refs_psi inclui valores ABAIXO e ACIMA do supercharge da Parte A
    # (~82 psi), para tornar visivel o trade-off no relatorio.
    refs_psi = [50, 60, 70, 80, 100, 150, 200, 250, 300]

    results = []
    for ref_psi in refs_psi:
        ref_pa = ref_psi * PSI_PA
        for beta in betas:
            corr = predict_corrected(pts, cd_pooled, beta, ref_pa)
            s_global = summarize([r["err"] for r in corr])
            part_a = run_part_a_corrected(beta, ref_pa)
            s_a = summarize([r["err"] for r in part_a])
            results.append({"beta": beta, "ref_psi": ref_psi,
                            "mape_global": s_global["mape"], "mean_global": s_global["mean"],
                            "mape_a": s_a["mape"], "mean_a": s_a["mean"],
                            "protects_A": ref_psi <= supercharge_A_psi})

    print(f"\n{'ref[psi]':>8} {'beta':>5} {'global MAPE':>12} {'global mean':>12} "
          f"{'PartA MAPE':>11} {'PartA mean':>11}  {'ProtegeA?':>10}")
    for r in results:
        flag = "SIM (gate)" if r["protects_A"] else ("degrada" if r["mape_a"] > 3.6 else "-")
        print(f"{r['ref_psi']:8.0f} {r['beta']:5.1f} {r['mape_global']:12.2f} "
              f"{r['mean_global']:+12.2f} {r['mape_a']:11.2f} {r['mean_a']:+11.2f}  {flag:>10}")

    # --- Duas leituras separadas, para deixar o trade-off explicito ---
    print(f"\n{'='*78}")
    protected = [r for r in results if r["protects_A"]]
    if protected:
        best_protected = min(protected, key=lambda r: r["mape_global"])
        print("OPCAO CONSERVADORA -- Parte A protegida por CONSTRUCAO "
              f"(supercharge_ref <= {supercharge_A_psi:.0f} psi):")
        print(f"  beta = {best_protected['beta']}, "
              f"supercharge_ref = {best_protected['ref_psi']} psi")
        print(f"  MAPE global: {best_protected['mape_global']:.2f}% "
              f"(baseline: {summarize([r['err'] for r in base])['mape']:.2f}%)")
        print(f"  MAPE Parte A: {best_protected['mape_a']:.2f}% "
              f"(identico ao baseline, por construcao)")

    print()
    unrestricted = [r for r in results if r["mape_a"] <= 5.0]
    if unrestricted:
        best_unrestricted = min(unrestricted, key=lambda r: r["mape_global"])
        print("OPCAO AGRESSIVA -- so exige MAPE_A <= 5% (Parte A pode "
              "degradar um pouco):")
        print(f"  beta = {best_unrestricted['beta']}, "
              f"supercharge_ref = {best_unrestricted['ref_psi']} psi")
        print(f"  MAPE global: {best_unrestricted['mape_global']:.2f}% "
              f"(baseline: {summarize([r['err'] for r in base])['mape']:.2f}%)")
        print(f"  MAPE Parte A: {best_unrestricted['mape_a']:.2f}% "
              f"(baseline: {summarize([r['err'] for r in part_a_base])['mape']:.2f}%)")

        best_corr = predict_corrected(pts, cd_pooled, best_unrestricted["beta"],
                                       best_unrestricted["ref_psi"] * PSI_PA)
        print(f"\n  Breakdown por supercharge (opcao agressiva):")
        for lo, hi, name in [(0, 100, "< 100 psi"), (100, 200, "100-200 psi"),
                             (200, 1000, ">= 200 psi")]:
            b = [r["err"] for r in best_corr if lo <= r["super_psi"] < hi]
            print(f"    {name:<15}", fmt(summarize(b)))

        # -------------------------------------------------------------
        # TESTE: quanto da degradacao da Parte A e efeito de kappa vs.
        # efeito do Cd assumido (0.65) ja nao ser o correto para comecar?
        # Cd_FIG15 (0.681) e o valor MEDIDO diretamente para este injector
        # (Waxman Fig.15, digitised) -- ja se sabia, mesmo antes desta
        # correcao, que so trocar o Cd fazia o MAPE baseline subir de
        # 2.76% para 4.41% (ver waxman_2013_validation.py print_part_d()
        # e a sensibilidade ja feita la).
        # -------------------------------------------------------------
        print(f"\n{'='*78}")
        print("ISOLAR EFEITO DE Cd vs. EFEITO DE KAPPA na Parte A")
        print(f"{'='*78}")

        beta_best = best_unrestricted["beta"]
        ref_best_pa = best_unrestricted["ref_psi"] * PSI_PA

        combos = [
            ("Cd=0.65 (assumido), beta=0 (sem correcao)", CD_PART_A, 0.0, 200 * PSI_PA),
            ("Cd=0.65 (assumido), beta=best (com correcao)", CD_PART_A, beta_best, ref_best_pa),
            ("Cd=0.681 (medido, Fig.15), beta=0 (sem correcao)", CD_PART_A_FIG15, 0.0, 200 * PSI_PA),
            ("Cd=0.681 (medido, Fig.15), beta=best (com correcao)", CD_PART_A_FIG15, beta_best, ref_best_pa),
        ]
        for label, cd_test, beta_test, ref_test in combos:
            rows = run_part_a_corrected(beta_test, ref_test, cd=cd_test)
            s = summarize([r["err"] for r in rows])
            print(f"  {label:<52}", fmt(s))

        print(f"\n  Leitura: compara a 1a linha com a 3a (efeito de trocar so o Cd,")
        print(f"  SEM correcao de kappa) contra a 2a vs a 4a (mesmo efeito, MAS com")
        print(f"  a correcao ja aplicada). Se a diferenca 3a-1a for parecida com a")
        print(f"  diferenca 4a-2a, o Cd explica a maior parte do 'aumento' de MAPE")
        print(f"  que atribuimos a correcao de kappa -- nao a correcao em si.")


if __name__ == "__main__":
    main()
