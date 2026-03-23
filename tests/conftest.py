from pathlib import Path
import sys

TESTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = TESTS_DIR.parent
SRC_ROOT = REPO_ROOT / "src"

for path in (TESTS_DIR, SRC_ROOT):
    path_str = str(path)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)
