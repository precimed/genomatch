from pathlib import Path
import sys

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent
SRC_MODULE_DIR = REPO_ROOT / "src" / "genomatch"

for path in (TESTS_DIR, SRC_MODULE_DIR):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)
