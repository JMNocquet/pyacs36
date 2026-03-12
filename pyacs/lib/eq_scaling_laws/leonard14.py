"""
leonard14_full.py

Complete implementation scaffold for Leonard (2014) earthquake scaling
relations (Table 3), covering all six deformation regimes:

 - interplate_ss
 - interplate_ds
 - crustal_ss
 - crustal_ds
 - scr_ss
 - scr_ds

All equations normalized to the form:

    log10(X) = a + b * log10(Y)

Units follow Leonard (2014):

    L, W in kilometers (km)
    A in square kilometers (km²)
    D, Dav in meters (m)
    M0 in Newton-meters (N·m)
    LSR / SRL dimensionless

Domain rules:
    - Upper branches use >= threshold.
    - Branch selection uses the domain variable printed in Leonard’s table
      (L, A, W, D, or LSR).

Alias system:
    D_avg  -> Dav
    D_mean -> Dav
    D_max  -> D
    SRL    -> LSR
"""

import math

# -----------------------------------------------------------
# VARIABLE ALIASES
# -----------------------------------------------------------

ALIASES = {
    "D_avg": "Dav",
    "D_mean": "Dav",
    "D_max": "D",
    "SRL": "LSR",
}

def _resolve(var):
    """Resolve aliases to the canonical Leonard variable names."""
    return ALIASES.get(var, var)

# -----------------------------------------------------------
# MAIN DICTIONARY
# This will be populated in Blocks 2–6
# -----------------------------------------------------------

LEONARD14 = {
    "interplate_ss": {},
    "interplate_ds": {},
    "crustal_ss": {},
    "crustal_ds": {},
    "scr_ss": {},
    "scr_ds": {},
}

# ============================================================
# INTERPLATE SS (Strike-Slip at subduction interface)
# ============================================================

LEONARD14["interplate_ss"] = {

    # --------------------------------------------------------
    # log(A) = a + b log(L)
    # Domains: L < 3400, L >= 3400
    # --------------------------------------------------------
    ("A", "L"): [
        {"a": 0.0,     "b": 2.0,     "domain": ("L", 0, 3400)},
        {"a": 1.176,   "b": 1.667,   "domain": ("L", 3400, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(L)
    # same L thresholds
    # --------------------------------------------------------
    ("W", "L"): [
        {"a": 0.0,     "b": 1.0,     "domain": ("L", 0, 3400)},
        {"a": 0.588,   "b": 0.833,   "domain": ("L", 3400, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(A)
    # Domains: A < 11.5e6, A >= 11.5e6
    # --------------------------------------------------------
    ("W", "A"): [
        {"a": 0.0,     "b": 0.5,     "domain": ("A", 0, 1.15e7)},
        {"a": 0.393,   "b": 0.4,     "domain": ("A", 1.15e7, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(A)
    # For interplate SS: one branch only
    # --------------------------------------------------------
    ("Dav", "A"): [
        {"a": -4.121,  "b": 0.5,     "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(L)
    # same L thresholds: <3400, >=3400
    # --------------------------------------------------------
    ("Dav", "L"): [
        {"a": -4.121,  "b": 1.0,     "domain": ("L", 0, 3400)},
        {"a": -3.528,  "b": 0.833,   "domain": ("L", 3400, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(W)
    # one branch? Screenshot shows one value.
    # --------------------------------------------------------
    ("Dav", "W"): [
        {"a": -4.829,  "b": 1.25,    "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(A)
    # One branch
    # --------------------------------------------------------
    ("M0", "A"): [
        {"a": 5.769,   "b": 1.5,     "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(L)
    # same L thresholds
    # --------------------------------------------------------
    ("M0", "L"): [
        {"a": 5.769,   "b": 3.0,     "domain": ("L", 0, 3400)},
        {"a": 7.564,   "b": 2.5,     "domain": ("L", 3400, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(W)
    # one branch
    # --------------------------------------------------------
    ("M0", "W"): [
        {"a": 2.935,   "b": 3.75,    "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(D)
    # one branch
    # --------------------------------------------------------
    ("M0", "D"): [
        {"a": 14.307,  "b": 3.0,     "domain": ("D", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(LSR)
    # one branch
    # --------------------------------------------------------
    ("M0", "LSR"): [
        {"a": 8.86,    "b": 2.27,    "domain": ("LSR", 0, None)},
    ],
}

# ============================================================
# INTERPLATE DS (Dip-Slip at subduction interface)
# ============================================================

LEONARD14["interplate_ds"] = {

    # --------------------------------------------------------
    # log(A) = a + b log(L)
    # Domains: L < 5360, L >= 5360
    # --------------------------------------------------------
    ("A", "L"): [
        {"a": 0.0,      "b": 2.0,      "domain": ("L", 0, 5360)},
        {"a": 1.243,    "b": 1.667,    "domain": ("L", 5360, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(L)
    # same L thresholds
    # --------------------------------------------------------
    ("W", "L"): [
        {"a": 0.0,      "b": 1.0,      "domain": ("L", 0, 5360)},
        {"a": 1.243,    "b": 0.667,    "domain": ("L", 5360, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(A)
    # Domains: A < 28.7e6, A >= 28.7e6
    # --------------------------------------------------------
    ("W", "A"): [
        {"a": 0.0,      "b": 0.5,      "domain": ("A", 0, 2.87e7)},
        {"a": 0.746,    "b": 0.4,      "domain": ("A", 2.87e7, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(A)
    # one branch
    # --------------------------------------------------------
    ("Dav", "A"): [
        {"a": -4.420,   "b": 0.5,      "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(L)
    # Domains: L < 5360, L >= 5360
    # --------------------------------------------------------
    ("Dav", "L"): [
        {"a": -4.420,   "b": 1.0,      "domain": ("L", 0, 5360)},
        {"a": -3.799,   "b": 0.833,    "domain": ("L", 5360, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(W)
    # One branch: slope 1.25, intercept -5.352
    # --------------------------------------------------------
    ("Dav", "W"): [
        {"a": -5.352,   "b": 1.25,     "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(A)
    # one branch
    # --------------------------------------------------------
    ("M0", "A"): [
        {"a": 6.098,    "b": 1.5,      "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(L)
    # Domains: L < 5360, L >= 5360
    # --------------------------------------------------------
    ("M0", "L"): [
        {"a": 6.098,    "b": 3.0,      "domain": ("L", 0, 5360)},
        {"a": 7.963,    "b": 2.5,      "domain": ("L", 5360, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(W)
    # One branch: slope 3.75, intercept 3.301
    # --------------------------------------------------------
    ("M0", "W"): [
        {"a": 3.301,    "b": 3.75,     "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(D)
    # One branch
    # --------------------------------------------------------
    ("M0", "D"): [
        {"a": 14.939,   "b": 3.0,      "domain": ("D", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(LSR)
    # One branch: slope 2.27, intercept 9.25
    # --------------------------------------------------------
    ("M0", "LSR"): [
        {"a": 9.25,     "b": 2.27,     "domain": ("LSR", 0, None)},
    ],
}

# ============================================================
# CRUSTAL SS (Crustal Strike-Slip Faulting)
# ============================================================

LEONARD14["crustal_ss"] = {

    # --------------------------------------------------------
    # log(A) = a + b log(L)
    # Domains: L < 400, L >= 400
    # --------------------------------------------------------
    ("A", "L"): [
        {"a": 0.0,       "b": 2.0,       "domain": ("L", 0, 400)},
        {"a": 1.046,     "b": 1.667,     "domain": ("L", 400, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(L)
    # same L thresholds
    # --------------------------------------------------------
    ("W", "L"): [
        {"a": 0.0,       "b": 1.0,       "domain": ("L", 0, 400)},
        {"a": 0.523,     "b": 0.833,     "domain": ("L", 400, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(A)
    # Domains: A < 11.5e6, A >= 11.5e6
    # (same threshold as interplate SS)
    # --------------------------------------------------------
    ("W", "A"): [
        {"a": 0.0,       "b": 0.5,       "domain": ("A", 0, 1.15e7)},
        {"a": 0.393,     "b": 0.4,       "domain": ("A", 1.15e7, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(A)
    # one branch
    # --------------------------------------------------------
    ("Dav", "A"): [
        {"a": -3.983,    "b": 0.5,       "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(L)
    # same L thresholds
    # --------------------------------------------------------
    ("Dav", "L"): [
        {"a": -3.983,    "b": 1.0,       "domain": ("L", 0, 400)},
        {"a": -3.460,    "b": 0.833,     "domain": ("L", 400, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(W)
    # one branch
    # --------------------------------------------------------
    ("Dav", "W"): [
        {"a": -4.687,    "b": 1.25,      "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(A)
    # one branch
    # --------------------------------------------------------
    ("M0", "A"): [
        {"a": 6.016,     "b": 1.5,       "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(L)
    # Domains: L < 400, L >= 400
    # --------------------------------------------------------
    ("M0", "L"): [
        {"a": 6.016,     "b": 3.0,       "domain": ("L", 0, 400)},
        {"a": 7.559,     "b": 2.5,       "domain": ("L", 400, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(W)
    # one branch
    # --------------------------------------------------------
    ("M0", "W"): [
        {"a": 3.060,     "b": 3.75,      "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(D)
    # one branch
    # --------------------------------------------------------
    ("M0", "D"): [
        {"a": 14.972,    "b": 3.0,       "domain": ("D", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(LSR)
    # one branch
    # --------------------------------------------------------
    ("M0", "LSR"): [
        {"a": 8.876,     "b": 2.27,      "domain": ("LSR", 0, None)},
    ],
}

# ============================================================
# CRUSTAL DS (Crustal Dip-Slip: Reverse + Normal)
# ============================================================

LEONARD14["crustal_ds"] = {

    # --------------------------------------------------------
    # log(A) = a + b log(L)
    # Domains: L < 200, L >= 200
    # --------------------------------------------------------
    ("A", "L"): [
        {"a": 0.0,       "b": 2.0,       "domain": ("L", 0, 200)},
        {"a": 0.699,     "b": 1.667,     "domain": ("L", 200, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(L)
    # same L thresholds (0–200, >=200)
    # --------------------------------------------------------
    ("W", "L"): [
        {"a": 0.0,       "b": 1.0,       "domain": ("L", 0, 200)},
        {"a": 0.349,     "b": 0.833,     "domain": ("L", 200, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(A)
    # Domains identical to interplate SS & crustal SS: A < 11.5e6, >= 11.5e6
    # --------------------------------------------------------
    ("W", "A"): [
        {"a": 0.0,       "b": 0.5,       "domain": ("A", 0, 1.15e7)},
        {"a": 0.393,     "b": 0.4,       "domain": ("A", 1.15e7, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(A)
    # one branch
    # --------------------------------------------------------
    ("Dav", "A"): [
        {"a": -4.008,    "b": 0.5,       "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(L)
    # Domains L < 200, L >= 200
    # --------------------------------------------------------
    ("Dav", "L"): [
        {"a": -4.008,    "b": 1.0,       "domain": ("L", 0, 200)},
        {"a": -3.659,    "b": 0.833,     "domain": ("L", 200, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(W)
    # one branch
    # --------------------------------------------------------
    ("Dav", "W"): [
        {"a": -4.712,    "b": 1.25,      "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(A)
    # one branch
    # --------------------------------------------------------
    ("M0", "A"): [
        {"a": 6.021,     "b": 1.5,       "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(L)
    # Domains L < 200, >= 200
    # --------------------------------------------------------
    ("M0", "L"): [
        {"a": 6.021,     "b": 3.0,       "domain": ("L", 0, 200)},
        {"a": 7.522,     "b": 2.5,       "domain": ("L", 200, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(W)
    # one branch
    # --------------------------------------------------------
    ("M0", "W"): [
        {"a": 2.957,     "b": 3.75,      "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(D)
    # one branch
    # --------------------------------------------------------
    ("M0", "D"): [
        {"a": 14.940,    "b": 3.0,       "domain": ("D", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(LSR)
    # one branch
    # --------------------------------------------------------
    ("M0", "LSR"): [
        {"a": 8.823,     "b": 2.27,      "domain": ("LSR", 0, None)},
    ],
}

# ============================================================
# SCR_SS  (Stable Continental Region — Strike-Slip)
# ============================================================

LEONARD14["scr_ss"] = {

    # --------------------------------------------------------
    # log(A) = a + b log(L)
    # Domains: L < 900, L >= 900
    # --------------------------------------------------------
    ("A", "L"): [
        {"a": 0.0,       "b": 2.0,       "domain": ("L", 0, 900)},
        {"a": 1.176,     "b": 1.667,     "domain": ("L", 900, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(L)
    # Domains: L < 900, L >= 900
    # --------------------------------------------------------
    ("W", "L"): [
        {"a": 0.0,       "b": 1.0,       "domain": ("L", 0, 900)},
        {"a": 0.588,     "b": 0.833,     "domain": ("L", 900, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(A)
    # Threshold: A < 11.5e6 km²
    # --------------------------------------------------------
    ("W", "A"): [
        {"a": 0.0,       "b": 0.5,       "domain": ("A", 0, 1.15e7)},
        {"a": 0.393,     "b": 0.4,       "domain": ("A", 1.15e7, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(A)
    # one branch
    # --------------------------------------------------------
    ("Dav", "A"): [
        {"a": -4.121,    "b": 0.5,       "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(L)
    # Domains: <900, >=900
    # --------------------------------------------------------
    ("Dav", "L"): [
        {"a": -4.121,    "b": 1.0,       "domain": ("L", 0, 900)},
        {"a": -3.528,    "b": 0.833,     "domain": ("L", 900, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(W)
    # one branch
    # --------------------------------------------------------
    ("Dav", "W"): [
        {"a": -4.829,    "b": 1.25,      "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(A)
    # one branch
    # --------------------------------------------------------
    ("M0", "A"): [
        {"a": 5.769,     "b": 1.5,       "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(L)
    # Domains: <900, >=900
    # --------------------------------------------------------
    ("M0", "L"): [
        {"a": 5.769,     "b": 3.0,       "domain": ("L", 0, 900)},
        {"a": 7.564,     "b": 2.5,       "domain": ("L", 900, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(W)
    # one branch
    # --------------------------------------------------------
    ("M0", "W"): [
        {"a": 2.935,     "b": 3.75,      "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(D)
    # one branch
    # --------------------------------------------------------
    ("M0", "D"): [
        {"a": 14.307,    "b": 3.0,       "domain": ("D", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(LSR)
    # one branch
    # --------------------------------------------------------
    ("M0", "LSR"): [
        {"a": 8.86,      "b": 2.27,      "domain": ("LSR", 0, None)},
    ],
}

# ============================================================
# SCR_DS  (Stable Continental Region — Dip-Slip)
# ============================================================

LEONARD14["scr_ds"] = {

    # --------------------------------------------------------
    # log(A) = a + b log(L)
    # Domains: L < 1000, L >= 1000
    # --------------------------------------------------------
    ("A", "L"): [
        {"a": 0.0,       "b": 2.0,       "domain": ("L", 0, 1000)},
        {"a": 1.176,     "b": 1.667,     "domain": ("L", 1000, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(L)
    # same thresholds
    # --------------------------------------------------------
    ("W", "L"): [
        {"a": 0.0,       "b": 1.0,       "domain": ("L", 0, 1000)},
        {"a": 0.588,     "b": 0.833,     "domain": ("L", 1000, None)},
    ],

    # --------------------------------------------------------
    # log(W) = a + b log(A)
    # Domains: A < 11.5e6, >=11.5e6
    # --------------------------------------------------------
    ("W", "A"): [
        {"a": 0.0,       "b": 0.5,       "domain": ("A", 0, 1.15e7)},
        {"a": 0.393,     "b": 0.4,       "domain": ("A", 1.15e7, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(A)
    # one branch
    # --------------------------------------------------------
    ("Dav", "A"): [
        {"a": -4.121,    "b": 0.5,       "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(L)
    # Domains: <1000, >=1000
    # --------------------------------------------------------
    ("Dav", "L"): [
        {"a": -4.121,    "b": 1.0,       "domain": ("L", 0, 1000)},
        {"a": -3.528,    "b": 0.833,     "domain": ("L", 1000, None)},
    ],

    # --------------------------------------------------------
    # log(Dav) = a + b log(W)
    # one branch
    # --------------------------------------------------------
    ("Dav", "W"): [
        {"a": -4.829,    "b": 1.25,      "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(A)
    # one branch
    # --------------------------------------------------------
    ("M0", "A"): [
        {"a": 5.769,     "b": 1.5,       "domain": ("A", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(L)
    # Domains: <1000, >=1000
    # --------------------------------------------------------
    ("M0", "L"): [
        {"a": 5.769,     "b": 3.0,       "domain": ("L", 0, 1000)},
        {"a": 7.564,     "b": 2.5,       "domain": ("L", 1000, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(W)
    # one branch
    # --------------------------------------------------------
    ("M0", "W"): [
        {"a": 2.935,     "b": 3.75,      "domain": ("W", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(D)
    # one branch
    # --------------------------------------------------------
    ("M0", "D"): [
        {"a": 14.307,    "b": 3.0,       "domain": ("D", 0, None)},
    ],

    # --------------------------------------------------------
    # log(M0) = a + b log(LSR)
    # one branch
    # --------------------------------------------------------
    ("M0", "LSR"): [
        {"a": 8.86,      "b": 2.27,      "domain": ("LSR", 0, None)},
    ],
}


# -----------------------------------------------------------
# HELPER: SELECT REGRESSION BRANCH
# -----------------------------------------------------------

def _select_branch(branches, value):
    """
    Given a list of regression branches of the form:
        {
            "a": float,
            "b": float,
            "domain": (varname, lo, hi)
        }
    return (a, b) for the branch whose domain matches the value.
    """
    for br in branches:
        dom_var, lo, hi = br["domain"]
        lo_ok = (lo is None) or (value >= lo)
        hi_ok = (hi is None) or (value < hi)
        if lo_ok and hi_ok:
            return br["a"], br["b"]
    raise ValueError(f"No matching branch for value {value}")

# -----------------------------------------------------------
# PUBLIC API: PREDICT log10(X)
# -----------------------------------------------------------

def predict(X, Y, Y_value, regime):
    """
    Predict log10(X) using Leonard (2014) regression:

        log10(X) = a + b * log10(Y_value)

    Parameters
    ----------
    X : str
        Dependent variable (A, W, Dav, D, M0).
    Y : str
        Independent variable (L, A, W, D, LSR).
    Y_value : float
        The numeric value of Y in Leonard's units.
    regime : str
        One of the six regimes in LEONARD14.

    Returns
    -------
    float
        log10(X)

    Raises
    ------
    NotImplementedError
        If the regime has not been populated.
    ValueError
        If no branch matches the given value.
    """

    X = _resolve(X)
    Y = _resolve(Y)

    if regime not in LEONARD14 or not LEONARD14[regime]:
        raise NotImplementedError(f"Regime '{regime}' has no data loaded.")

    # Search for matching (X, Y)
    for (XX, YY), branches in LEONARD14[regime].items():
        if XX == X and YY == Y:

            # Identify domain variable
            dom_var = branches[0]["domain"][0]
            dom_var = _resolve(dom_var)

            # Use Y_value as domain value if domain var == Y
            # (Leonard's table defines which variable controls the split)
            if dom_var == Y:
                domain_value = Y_value
            else:
                # If domain variable is not Y, user must provide it
                raise ValueError(
                    f"Domain for ({X},{Y}) regression is based on {dom_var}, "
                    f"but only value for {Y} was provided."
                )

            a, b = _select_branch(branches, domain_value)
            return a + b * math.log10(Y_value)

    raise ValueError(f"No regression for ({X},{Y}) in regime {regime}.")

# -----------------------------------------------------------
# OPTIONAL: RETURN X IN LINEAR SPACE
# -----------------------------------------------------------

def predict_linear(X, Y, Y_value, regime):
    """
    Return X (not log10(X)).

    X = 10^(log10(X))
    """
    return 10 ** predict(X, Y, Y_value, regime)

def predict_inverse(Y, X, X_value, regime):
    """
    Invert a relationship: if we have X, predict Y.
    Given log10(X) = a + b*log10(Y), solve for Y:
    log10(Y) = (log10(X_value) - a) / b
    Y = 10^((log10(X_value) - a) / b)
    
    Parameters
    ----------
    Y : str
        Variable to predict (e.g., 'L', 'W', 'A')
    X : str
        Variable we have (e.g., 'M0')
    X_value : float
        Value of X
    regime : str
        Mechanism regime
    """
    X = _resolve(X)
    Y = _resolve(Y)
    
    if regime not in LEONARD14 or not LEONARD14[regime]:
        raise NotImplementedError(f"Regime '{regime}' has no data loaded.")
    
    # Search for matching (X, Y) relationship where X is predicted from Y
    # We want to invert this to get Y from X
    for (XX, YY), branches in LEONARD14[regime].items():
        if XX == X and YY == Y:
            # For inversion, we need to handle domain branches
            # Try each branch and see which one gives a valid result
            best_Y = None
            
            for branch in branches:
                a = branch["a"]
                b = branch["b"]
                domain_var, domain_min, domain_max = branch["domain"]
                domain_var = _resolve(domain_var)
                
                # Invert: log10(Y) = (log10(X_value) - a) / b
                if abs(b) < 1e-10:
                    continue  # Skip if b is essentially zero
                
                log10_Y = (math.log10(X_value) - a) / b
                Y_candidate = 10 ** log10_Y
                
                # Check if Y_candidate is in the domain for this branch
                if domain_var == Y:
                    if domain_min is not None and Y_candidate < domain_min:
                        continue
                    if domain_max is not None and Y_candidate >= domain_max:
                        continue
                    best_Y = Y_candidate
                    break
                else:
                    # Domain is based on another variable, use this as candidate
                    # We'll need to verify later, but for now take first valid
                    if best_Y is None:
                        best_Y = Y_candidate
            
            if best_Y is not None:
                return best_Y
    
    raise ValueError(f"Cannot invert relationship ({X},{Y}) for regime {regime}.")

# -----------------------------------------------------------
# COMPUTE SLIP, LENGTH, WIDTH FROM Mw
# -----------------------------------------------------------

def predict_leonard14(Mw, mechanism='all'):
    """
    Compute slip, length, width from Mw using Leonard (2014) regression.

    Parameters
    ----------
    Mw : float
        Moment magnitude.
    mechanism : str
        One of the six regimes in LEONARD14.
    """

    from pyacs.lib.units import magnitude_to_moment
    M0 = magnitude_to_moment(Mw)

    import pyacs.lib.eq_scaling_laws.wc94 as wc94 

    
    # Compute length and width from M0 by inverting relationships
    # Try to compute A from M0, then L and W from A
    try:
        A = predict_inverse('A', 'M0', M0, mechanism)
        length = predict_inverse('L', 'A', A, mechanism)
        width = predict_inverse('W', 'A', A, mechanism)
        slip = predict_linear('Dav', 'A', A, mechanism)
    except (ValueError, NotImplementedError):
        # Fallback: compute directly from M0 if available
        try:
            length = predict_inverse('L', 'M0', M0, mechanism)
        except (ValueError, NotImplementedError):
            raise ValueError(f"Cannot compute length for mechanism {mechanism}")
        
        try:
            width = predict_inverse('W', 'M0', M0, mechanism)
        except (ValueError, NotImplementedError):
            raise ValueError(f"Cannot compute width for mechanism {mechanism}")
        try:
            slip = predict_linear('Dav', 'A', A, mechanism)
        except (ValueError, NotImplementedError):
            raise ValueError(f"Cannot compute slip for mechanism {mechanism}")
    return slip, length*1.E-3, width * 1.E-3 # km