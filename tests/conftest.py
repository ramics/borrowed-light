"""AC6: the full synthetic suite must run in under five minutes on a laptop.
Enforced here so a slow change fails CI instead of quietly eroding it."""

import time

import pytest

BUDGET_SECONDS = 300
_start = None


def pytest_sessionstart(session):
    global _start
    _start = time.monotonic()


def pytest_sessionfinish(session, exitstatus):
    elapsed = time.monotonic() - _start
    full_run = not session.config.args or session.config.args == ["tests"]
    if full_run and elapsed > BUDGET_SECONDS and exitstatus == 0:
        session.exitstatus = pytest.ExitCode.TESTS_FAILED
        print(f"\nAC6 FAILED: suite took {elapsed:.0f}s, budget {BUDGET_SECONDS}s")
