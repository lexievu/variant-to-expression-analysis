"""Shared pytest fixtures and configuration."""

import logging
import pytest


@pytest.fixture(autouse=True)
def _disable_log_capture_for_logging_tests(request):
    """For tests marked @pytest.mark.no_capture, remove all root logger
    handlers (including pytest's LogCaptureHandler) before the test runs
    and restore them afterwards.  This lets tests inspect handler state
    set up by our own setup_logging() without interference.
    """
    marker = request.node.get_closest_marker("no_capture")
    if marker is None:
        yield
        return

    root = logging.getLogger()
    saved = root.handlers[:]
    saved_level = root.level
    root.handlers.clear()

    yield

    # Teardown — remove anything the test added and restore originals
    for h in root.handlers[:]:
        root.removeHandler(h)
        h.close()
    root.setLevel(saved_level)
    for h in saved:
        root.addHandler(h)
