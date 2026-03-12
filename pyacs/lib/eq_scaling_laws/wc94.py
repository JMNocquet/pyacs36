
import numpy as np

# Wells & Coppersmith 1994 full tables (a,b,sigma)
WC94 = {
    'length': {
        'all':    {'a': -3.22, 'b': 0.69, 'sigma': 0.22},
        'strike': {'a': -3.55, 'b': 0.74, 'sigma': 0.23},
        'reverse':{'a': -2.42, 'b': 0.58, 'sigma': 0.25},
        'normal': {'a': -2.57, 'b': 0.62, 'sigma': 0.29},
    },
    'width': {
        'all':    {'a': -1.01, 'b': 0.32, 'sigma': 0.22},
        'strike': {'a': -1.14, 'b': 0.35, 'sigma': 0.26},
        'reverse':{'a': -0.76, 'b': 0.27, 'sigma': 0.22},
        'normal': {'a': -1.61, 'b': 0.41, 'sigma': 0.29},
    },
    'avg_slip': {
        'strike': {'a': -6.32, 'b': 0.90, 'sigma': 0.30},
        'reverse': {'a': -0.74, 'b': 0.08, 'sigma': 0.30},
        'normal': {'a': -4.45, 'b': 0.63, 'sigma': 0.30},
        'all': {'a': -4.80, 'b': 0.69, 'sigma': 0.30},
    },
    'max_slip': {
            'strike': {'a': -7.03, 'b': 1.03, 'sigma': 0.30},
            'reverse': {'a': -1.84, 'b': 0.29, 'sigma': 0.30},
            'normal': {'a': -5.90, 'b': 0.89, 'sigma': 0.30},
            'all': {'a': -5.46, 'b': 0.82, 'sigma': 0.30},

    }
}

def _compute(Mw, key, mech='all'):
    c = WC94[key][mech]
    logX = c['a'] + c['b'] * Mw
    X = 10**logX
    s = c['sigma']
    return X, 10**(logX - s), 10**(logX + s)

def length(Mw, mechanism='all'): return _compute(Mw,'length',mechanism)
def width(Mw, mechanism='all'): return _compute(Mw,'width',mechanism)
def avg_slip(Mw, mechanism='all'): return _compute(Mw,'avg_slip',mechanism)
def max_slip(Mw, mechanism='all'): return _compute(Mw,'max_slip',mechanism)

