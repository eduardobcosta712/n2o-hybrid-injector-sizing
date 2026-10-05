"""
audit_repo.py -- consistency audit of the repository (keep in the repo: tools/audit_repo.py)

Run from the repository root, after any change to physics, results or documents:
    python tools/audit_repo.py

Checks
  1. README tree  vs  the files really tracked by git (both directions)
  2. expected files present (report, figures, kappa_correction folder, option tests)
  3. test counts: pytest collection  vs  the numbers quoted in the README banner
  4. stale statements (old test count, moved scripts, retracted claims, ...)
  5. paths quoted in backticks in the .md files that do not exist
  6. code hooks of the kappa-correction feature (solver option, app warning, PDF note)
  7. version: CITATION.cff vs CHANGELOG

Exit code 1 if any FAIL. WARN lines need a human look (they can be legitimate).
"""
import os
import re
import subprocess
import sys
from collections import Counter

ROOT = os.getcwd()
fails, warns = [], []


def ok(msg): print("  ok    ", msg)
def fail(msg): fails.append(msg); print("  FAIL  ", msg)
def warn(msg): warns.append(msg); print("  WARN  ", msg)
def head(t): print("\n" + t)


def read(path):
    with open(os.path.join(ROOT, path), encoding="utf-8", errors="ignore") as f:
        return f.read()


def tracked_files():
    try:
        out = subprocess.run(["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True).stdout
        files = [l.strip() for l in out.splitlines() if l.strip()]
        if files: return files
    except Exception:
        pass
    skip = {".git", "__pycache__", "venv", ".venv", ".pytest_cache", ".streamlit", "node_modules"}
    out = []
    for d, dirs, fs in os.walk(ROOT):
        dirs[:] = [x for x in dirs if x not in skip]
        for f in fs:
            if not f.endswith((".pyc", ".aux", ".log", ".out", ".toc")):
                out.append(os.path.relpath(os.path.join(d, f), ROOT).replace("\\", "/"))
    return out


FILES = tracked_files()
TEXT = [f for f in FILES if f.endswith((".md", ".py", ".txt", ".yml", ".cff", ".tex")) and f != "tools/audit_repo.py"]

# ------------------------------------------------------------------ 1. README tree
head("1. README tree vs tracked files")
tree_files, tree_dirs = set(), set()
try:
    lines = read("README.md").splitlines()
    start = next(i for i, l in enumerate(lines) if l.strip().startswith("n2o-hybrid-injector-sizing/"))
    stack = []
    for l in lines[start + 1:]:
        if l.strip().startswith("```"): break
        m = re.search(r"(├── |└── )", l)
        if not m: continue
        depth = m.start() // 4 + 1
        name = l[m.end():].split()[0]
        is_dir = name.endswith("/")
        name = name.rstrip("/")
        stack = stack[:depth - 1] + [name]
        (tree_dirs if is_dir else tree_files).add("/".join(stack))
    missing = sorted(p for p in tree_files if p not in FILES)
    unlisted = sorted(p for p in FILES if p not in tree_files
                      and not any(p == d or p.startswith(d + "/") for d in ())
                      and p != "README.md")
    unlisted = [p for p in unlisted if p not in tree_files]
    if missing: fail("in the README tree but not tracked/existing: " + ", ".join(missing))
    else: ok(f"all {len(tree_files)} files listed in the tree exist")
    if unlisted: warn("tracked but not in the README tree: " + ", ".join(unlisted))
    else: ok("every tracked file is in the tree")
except StopIteration:
    fail("could not find the tree block in README.md")

# ------------------------------------------------------------------ 2. expected files
head("2. expected files")
EXPECT = ["docs/reports/dyer_supercharge_correction.tex", "docs/reports/dyer_supercharge_correction.pdf",
          "docs/reports/figures/make_figures.py",
          "validation/kappa_correction/README.md", "validation/kappa_correction/loocv_kappa_correction.py",
          "validation/kappa_correction/loocv_functional_forms.py", "validation/kappa_correction/scan_threshold_continuity.py",
          "validation/kappa_correction/loocv_kappa_results.txt", "validation/kappa_correction/loocv_functional_forms_results.txt",
          "validation/kappa_correction/scan_threshold_continuity_results.txt",
          "tests/test_kappa_correction_option.py", "tools/audit_repo.py"]
EXPECT += [f"docs/reports/figures/{n}.pdf" for n in ("fig_baseline_bands", "fig_factor_shapes", "fig_forms",
                                                     "fig_limit", "fig_loocv_bands", "fig_loocv_folds", "fig_scan")]
bad = [p for p in EXPECT if p not in FILES]
if bad: fail("missing: " + ", ".join(bad))
else: ok(f"all {len(EXPECT)} expected files are tracked")
for gone in ("validation/waxman_2013_experimental_data.csv",
             "apply_fixes_v1_0_1.py", "apply_kappa_option.py", "apply_docs_kappa.py", "apply_app_notice.py",
             "apply_readme_tree.py", "bump_release.py"):
    if gone in FILES: fail(f"should not be in the repo: {gone}")

# ------------------------------------------------------------------ 3. test counts
head("3. test counts: pytest vs README banner")
try:
    r = subprocess.run([sys.executable, "-m", "pytest", "--collect-only", "-q", "tests"], cwd=ROOT,
                       capture_output=True, text=True)
    per = Counter(m.group(1) for l in r.stdout.splitlines() if (m := re.match(r"tests[\\/](test_\w+)\.py::", l)))
    total = sum(per.values())
    if not total: warn("pytest collected nothing (is CoolProp installed? " + (r.stderr or r.stdout)[-200:].strip() + ")")
    else:
        banner = read("README.md")
        mt = re.search(r"\*\*(\d+) automated tests\*\*", banner)
        if not mt: warn("no '**N automated tests**' banner in README.md")
        elif int(mt.group(1)) != total: fail(f"README says {mt.group(1)} tests, pytest collects {total}")
        else: ok(f"README banner and pytest agree: {total} tests")
        for mod, n in re.findall(r"`(test_\w+)` (\d+)", banner):
            if per.get(mod) != int(n): fail(f"README says {mod} = {n}, pytest collects {per.get(mod)}")
        for mod in per:
            if f"`{mod}`" not in banner: warn(f"{mod} ({per[mod]} tests) is not in the README banner")
except Exception as e:
    warn(f"could not run pytest: {e}")

# ------------------------------------------------------------------ 4. stale statements
head("4. stale statements")
STALE = [
    (r"\b262[ -](automated )?tests?\b", "old test count (262)", ("CHANGELOG.md", "examples/", "validation/kappa_correction/")),
    (r"roughly half (of )?the Dyer flow", "'HEM = half of Dyer' (only true if SPI >> HEM)", ()),
    (r"not called by\s+`?full_system", "outdated: full_system.py now accepts kappa_correction", ()),
    (r"validation/(loocv_|scan_threshold)", "script moved to validation/kappa_correction/", ()),
    (r"waxman_2013_experimental_data", "removed file still referenced", ("docs/references.md", "CHANGELOG.md")),
    (r"python examples/generate_example_03", "wrong path (the script is in tests/)", ()),
]
found = False
for pat, why, skip in STALE:
    for f in TEXT:
        if any(f.startswith(s) for s in skip): continue
        for i, line in enumerate(read(f).splitlines(), 1):
            if re.search(pat, line):
                found = True; warn(f"{f}:{i}: {why}")
if not found: ok("none found")

# ------------------------------------------------------------------ 5. backticked paths
head("5. paths quoted in .md files that do not exist")
bn = {os.path.basename(f) for f in FILES}
miss = set()
for f in [x for x in TEXT if x.endswith(".md")]:
    for p in re.findall(r"`([A-Za-z0-9_./\\-]+\.(?:py|md|csv|txt|pdf|tex|yml|png|cff))`", read(f)):
        p = p.replace("\\", "/")
        if "*" in p or p.startswith(("results/", ".streamlit")): continue
        cands = [p, os.path.normpath(os.path.join(os.path.dirname(f), p)).replace("\\", "/")]
        if any(c in FILES for c in cands): continue
        if "/" not in p and p in bn: continue
        miss.add((f, p))
if miss:
    for f, p in sorted(miss): warn(f"{f}: `{p}` not found")
else: ok("all quoted paths resolve")

# ------------------------------------------------------------------ 6. code hooks
head("6. kappa-correction hooks in the code")
HOOKS = [("src/model/full_system.py", "kappa_correction", "solver option"),
         ("src/model/full_system.py", "_validate_kappa_correction", "argument validation"),
         ("src/interface/app.py", "render_kappa_notice_sizing", "app warning (Sizing)"),
         ("src/interface/app.py", "render_kappa_notice_design", "app warning (Design)"),
         ("src/interface/export.py", "Exploratory note", "PDF note")]
for path, needle, what in HOOKS:
    if path not in FILES: fail(f"{path} missing"); continue
    (ok if needle in read(path) else fail)(f"{what}: '{needle}' in {path}")

# ------------------------------------------------------------------ 7. version
head("7. version")
try:
    cff = re.search(r'^version:\s*"?([\d.]+)"?', read("CITATION.cff"), re.M).group(1)
    rel = re.findall(r"^## \[(\d[\d.]*)\]", read("CHANGELOG.md"), re.M)
    if "## [Unreleased]" in read("CHANGELOG.md"): warn("CHANGELOG still has an [Unreleased] section (run bump_release.py before tagging)")
    if rel and rel[0] != cff: fail(f"CITATION.cff version {cff} != latest CHANGELOG release {rel[0]}")
    else: ok(f"CITATION.cff version {cff} matches the CHANGELOG" if rel else f"CITATION.cff version {cff}")
except Exception as e:
    warn(f"could not compare versions: {e}")

print(f"\n{len(fails)} FAIL, {len(warns)} WARN")
sys.exit(1 if fails else 0)
