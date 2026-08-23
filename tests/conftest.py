"""
conftest.py

pytest configuration file. Automatically adds src/model/ to sys.path
so that all test files can import the model modules (n2o_properties,
feed_line, etc.) without needing to modify any existing source file
or install the package.

This file is loaded automatically by pytest before any test runs.
"""

import sys
import os

# Resolve the path to src/model/ relative to this file's location.
# conftest.py lives in tests/, so we go up one level to the repo root,
# then down into src/model/.
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_MODEL_DIR = os.path.join(_REPO_ROOT, "src", "model")

if _MODEL_DIR not in sys.path:
    sys.path.insert(0, _MODEL_DIR)
