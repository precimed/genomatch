import os
import subprocess
import sys
import sysconfig
from pathlib import Path

from utils import REPO_ROOT


def installed_site_packages(prefix: Path) -> Path:
    return Path(
        sysconfig.get_path(
            "purelib",
            vars={"base": str(prefix), "platbase": str(prefix)},
        )
    )


def installed_scripts_dir(prefix: Path) -> Path:
    return Path(
        sysconfig.get_path(
            "scripts",
            vars={"base": str(prefix), "platbase": str(prefix)},
        )
    )


def test_pip_install_exposes_cli_tools(tmp_path):
    prefix = tmp_path / "install"
    install_result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "install",
            "--no-build-isolation",
            "--no-deps",
            "--prefix",
            str(prefix),
            str(REPO_ROOT),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert install_result.returncode == 0, install_result.stderr

    scripts_dir = installed_scripts_dir(prefix)
    site_packages = installed_site_packages(prefix)
    prepare_variants = scripts_dir / "prepare_variants.py"
    project_payload = scripts_dir / "project_payload.py"
    assert prepare_variants.exists()
    assert project_payload.exists()

    env = os.environ.copy()
    env["PYTHONPATH"] = str(site_packages) + os.pathsep + env.get("PYTHONPATH", "")

    help_result = subprocess.run(
        [str(prepare_variants), "--help"],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    assert help_result.returncode == 0, help_result.stderr
    assert "Prepare raw input variants" in help_result.stdout

    import_result = subprocess.run(
        [sys.executable, "-c", "import genomatch; print(genomatch.__version__)"],
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )
    assert import_result.returncode == 0, import_result.stderr
    assert import_result.stdout.strip()
