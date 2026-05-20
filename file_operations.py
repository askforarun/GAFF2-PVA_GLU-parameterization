"""Local cleanup helpers for standalone parametrization workflows."""

from __future__ import annotations

import glob
import os
import shutil


# Signac metadata files are preserved when these helpers are used in a Signac
# workspace, but this module does not otherwise depend on Signac.
_SIGNAC_PRESERVE = frozenset({"signac_statepoint.json", "signac_job_document.json"})


def cleanup_beginning(target_dir: str | None = None) -> None:
    """Remove non-Python files from a working directory before regeneration.

    Subdirectories are left untouched. Python source files and Signac metadata
    are preserved so scripts can safely call this from their own directories.
    """
    if target_dir is None:
        target_dir = os.getcwd()
    target_dir = os.path.abspath(target_dir)

    files_removed = 0
    for entry in os.scandir(target_dir):
        if not entry.is_file(follow_symlinks=False):
            continue
        if entry.name.lower().endswith(".py"):
            continue
        if entry.name in _SIGNAC_PRESERVE:
            continue
        try:
            os.remove(entry.path)
            files_removed += 1
            print(f"Removed: {entry.name}")
        except OSError as exc:
            print(f"Could not remove {entry.name}: {exc}")

    print(
        "Beginning cleanup completed. "
        f"Removed {files_removed} files (Python files preserved)."
    )


def cleanup_end(target_dir: str | None = None) -> None:
    """Remove temporary files and sandbox directories after a run."""
    if target_dir is None:
        target_dir = os.getcwd()
    target_dir = os.path.abspath(target_dir)

    patterns_to_remove = [
        "leap.*",
        "*.out",
        "log*",
        "*.err",
        "ANTECHAMBER*",
        "sqm.*",
        "ATOMTYPE.INF",
        "*.AC",
        "*.AC0",
    ]

    files_removed = 0
    for pattern in patterns_to_remove:
        full_pattern = os.path.join(target_dir, pattern)
        for file_path in glob.glob(full_pattern):
            try:
                if os.path.isfile(file_path):
                    os.remove(file_path)
                    files_removed += 1
                    print(f"Removed temp file: {os.path.basename(file_path)}")
            except OSError as exc:
                print(f"Could not remove {file_path}: {exc}")

    sandbox_pattern = os.path.join(target_dir, "sandbox_*")
    for sandbox_dir in glob.glob(sandbox_pattern):
        if not os.path.isdir(sandbox_dir):
            continue
        try:
            shutil.rmtree(sandbox_dir)
            files_removed += 1
            print(f"Removed sandbox directory: {os.path.basename(sandbox_dir)}")
        except OSError as exc:
            print(f"Could not remove sandbox directory {sandbox_dir}: {exc}")

    print(
        "End cleanup completed. "
        f"Removed {files_removed} temporary files and directories."
    )
