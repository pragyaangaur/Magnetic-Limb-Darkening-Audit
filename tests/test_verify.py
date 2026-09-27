"""Run the fourteen published-number cross-checks in analysis/verify.py under pytest.

The script is the real record and stays runnable on its own. This wrapper only lets a
plain pytest run fail when a check fails, and names the check that did.
"""
import os
import sys

import pytest

ANALYSIS = os.path.join(os.path.dirname(__file__), os.pardir, "analysis")
sys.path.insert(0, os.path.abspath(ANALYSIS))

import verify  # noqa: E402


@pytest.fixture(scope="module")
def checks():
    verify.CHECKS.clear()
    verify.main()
    return list(verify.CHECKS)


def test_all_fourteen_checks_ran(checks):
    assert len(checks) == 14


def test_every_check_passes(checks):
    failed = [name for name, *_, ok in checks if not ok]
    assert not failed, f"failed checks: {failed}"
