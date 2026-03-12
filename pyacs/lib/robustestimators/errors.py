"""Exceptions for robustestimators module."""


class Error(Exception):
    """Base class for exceptions in module robustestimators."""

    pass


class UnboundedFunctionError(Exception):
    """Exception raised for unbounded objective function."""

    pass
