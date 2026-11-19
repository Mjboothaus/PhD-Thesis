"""Lightweight memory helpers used by the numerical core."""

from __future__ import annotations

import os
from collections.abc import Callable
from functools import wraps
from typing import ParamSpec, TypeVar

import psutil

P = ParamSpec("P")
R = TypeVar("R")


def get_memory_usage() -> float:
    """Return current process RSS in mebibytes."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024


def monitor_memory[**P, R](func: Callable[P, R]) -> Callable[P, R]:
    """No-op-friendly decorator placeholder for memory tracing.

    The historical Streamlit app logged deltas via loguru. Core library code keeps
    the decorator hook so call sites stay stable; logging stays optional and quiet
    by default so unit tests are not noisy.
    """

    @wraps(func)
    def wrapper(*args: P.args, **kwargs: P.kwargs) -> R:
        return func(*args, **kwargs)

    return wrapper
