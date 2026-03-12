"""
Unit tests for pyacs.message (message, verbose_message, error, warning, debug_message).
"""

import pytest

import pyacs.message as MSG


# -----------------------------------------------------------------------------
# Module exports
# -----------------------------------------------------------------------------


def test_message_module_exports():
    """Module exposes message, verbose_message, error, warning, debug_message."""
    assert hasattr(MSG, "message")
    assert hasattr(MSG, "verbose_message")
    assert hasattr(MSG, "error")
    assert hasattr(MSG, "warning")
    assert hasattr(MSG, "debug_message")
    assert callable(MSG.message)
    assert callable(MSG.verbose_message)
    assert callable(MSG.error)
    assert callable(MSG.warning)
    assert callable(MSG.debug_message)


# -----------------------------------------------------------------------------
# message
# -----------------------------------------------------------------------------


def test_message_level_zero():
    """message(ustr, level=0) runs without raising (prints [PYACS] prefix)."""
    MSG.message("test level 0", level=0)


def test_message_level_one():
    """message(ustr, level=1) runs without raising (banner, no colors)."""
    MSG.message("test level 1", level=1)


def test_message_level_two_requires_colors():
    """message(ustr, level=2) runs without raising when colors available."""
    pytest.importorskip("colors")
    MSG.message("test level 2", level=2)


# -----------------------------------------------------------------------------
# verbose_message
# -----------------------------------------------------------------------------


def test_verbose_message():
    """verbose_message(str) runs without raising."""
    MSG.verbose_message("verbose test")


# -----------------------------------------------------------------------------
# debug_message
# -----------------------------------------------------------------------------


def test_debug_message():
    """debug_message(str) runs without raising."""
    MSG.debug_message("debug test")


# -----------------------------------------------------------------------------
# warning
# -----------------------------------------------------------------------------


def test_warning_requires_colors():
    """warning(str) runs without raising when colors available."""
    pytest.importorskip("colors")
    MSG.warning("warning test")


# -----------------------------------------------------------------------------
# error
# -----------------------------------------------------------------------------


def test_error_exit_false():
    """error(str, exit=False) runs without raising (prints only)."""
    pytest.importorskip("colors")
    MSG.error("error test", exit=False)


def test_error_exit_true_raises_system_exit():
    """error(str, exit=True) raises SystemExit."""
    pytest.importorskip("colors")
    with pytest.raises(SystemExit):
        MSG.error("error exit test", exit=True)


def test_error_accepts_string():
    """error accepts various string inputs."""
    pytest.importorskip("colors")
    MSG.error("plain", exit=False)
    MSG.error("format %s" % 1, exit=False)
