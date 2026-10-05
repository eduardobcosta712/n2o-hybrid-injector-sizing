"""
apply_readme_tree.py   (one-off -- do not commit)

Completes the README repository tree: docs/reports/ (tex + figures) and tools/.
    python apply_readme_tree.py --dry-run | python apply_readme_tree.py
Same safety design as the other patch scripts.
"""
import os, re, sys

EDITS = [
    ("README.md", r"│       └── dyer_supercharge_correction.pdf",
     "│       ├── dyer_supercharge_correction.tex\n"
     "│       ├── dyer_supercharge_correction.pdf\n"
     "│       └── figures/\n"
     "│           ├── make_figures.py\n"
     "│           ├── fig_baseline_bands.pdf\n"
     "│           ├── fig_factor_shapes.pdf\n"
     "│           ├── fig_forms.pdf\n"
     "│           ├── fig_limit.pdf\n"
     "│           ├── fig_loocv_bands.pdf\n"
     "│           ├── fig_loocv_folds.pdf\n"
     "│           └── fig_scan.pdf"),
    ("README.md", r"├── tests/",
     "├── tools/\n│   └── audit_repo.py\n├── tests/"),
]

def flex(t): return r"\s+".join(re.escape(x) for x in t.split())

args = [a for a in sys.argv[1:] if not a.startswith("--")]
root = os.path.abspath(args[0] if args else ".")
text = {}; nl = {}; errs = []
for path, old, new in EDITS:
    full = os.path.join(root, path)
    if path not in text:
        if not os.path.exists(full): errs.append(f"{path}: not found"); text[path] = None; continue
        with open(full, encoding="utf-8", newline="") as f: s = f.read()
        nl[path] = "\r\n" if "\r\n" in s else "\n"; text[path] = s
    if text[path] is None: continue
    m = list(re.finditer(flex(old), text[path]))
    if len(m) != 1: errs.append(f"{path}: expected exactly 1 match, found {len(m)} for: {old[:50]}"); continue
    m = m[0]; text[path] = text[path][:m.start()] + new.replace("\n", nl[path]) + text[path][m.end():]
if errs:
    print("NOTHING WRITTEN:"); [print("  -", e) for e in errs]; sys.exit(1)
if "--dry-run" in sys.argv: print("Dry run OK."); sys.exit(0)
for path, s in text.items():
    with open(os.path.join(root, path), "w", encoding="utf-8", newline="") as f: f.write(s)
    print("patched:", path)
