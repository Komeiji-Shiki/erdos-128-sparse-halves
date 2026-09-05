"""Optional symbolic audit of the manuscript, using SymPy.

These are polynomial identities, checked by expansion, not parameter
sampling. The signs, graph constructions and external theorems still
need mathematical review. The main numerical replay uses only stdlib.
"""

import json
from pathlib import Path

import sympy as sp


def audit():
    a, b, z, cap, s, top, mu, d, psi, k, rho, cut = sp.symbols(
        "a b z A s T mu d psi k rho I")
    r = sp.Rational(1, 2)-b
    outside = 1-cap-b
    threshold = sp.Rational(1, 8)
    raw = (threshold*(b-z)+r*(1-r)/4+(r/2+threshold)*s
           -sp.Rational(3, 4)*s*s-threshold*outside/2)
    low = (12*cap-32*b*b+12*b-24*z+9)/192
    high = (cap-8*b*b+5*b-2*z)/16
    stationary = (r+sp.Rational(1, 4))/3
    c = 1-a
    identities = {
        "cut_baseline_after_clearing_denominators":
            4*c*c*(cut+a*(rho/2-cut))+(a*a-1)*mu
            -(4*c*c*(a*rho/2+c*cut)-(1+a)*c*mu),
        "reciprocal_tangent_equality":
            (19*d-18*mu)*mu*d*d+mu*d*d*(d+18*mu)-20*mu*d**3,
        "product_moment_factorization":
            d*psi*(k-psi+d)-k*d*d-d*(psi-d)*(k-psi),
        "curvature_low_after_multiplication_by_T":
            top*d*(top*d/2-mu*d)-top**2*(top*d/2)
            -(top*mu*(1-top)-2*mu**2)*d+top**2*mu-2*top*mu**2
            -(top-d)*(top/2-mu)*(2*mu-top*d),
        "curvature_high_after_multiplication_by_T":
            top*(d*mu-mu*d*d)-(top*mu*(1-top)-2*mu**2)*d-2*top*mu**2
            -mu*(top-d)*(top*d-2*mu),
        "extension_completion_of_square":
            raw-low+sp.Rational(3, 4)*(s-stationary)**2,
        "extension_high_endpoint": raw.subs(s, r)-high,
        "extension_low_derivative": sp.diff(low, b)+(16*b-3)/48,
        "extension_high_derivative": sp.diff(high, b)+(16*b-5)/16,
        "independence_J_minus_K":
            sp.Rational(9, 128)-z/8-(sp.Rational(3, 64)-z*z/6)
            -(z-sp.Rational(3, 8))**2/6,
        "extension_piece_continuity":
            (low-high).subs(b, sp.Rational(3, 8)),
        "large_neighbourhood_quadratic_maximum":
            threshold*(b-z)+b*(sp.Rational(1, 2)-b)/2
            -(sp.Rational(9, 128)-z/8-(b-sp.Rational(3, 8))**2/2),
    }
    for name, expression in identities.items():
        if sp.expand(expression) != 0:
            raise ValueError(f"Failed symbolic identity: {name}")
    return {"verified": True, "sympy_version": sp.__version__,
            "polynomial_identities": list(identities), "count": len(identities),
            "scope": "exact symbolic identities, not an end-to-end mathematical proof"}


if __name__ == "__main__":
    result = audit()
    root = Path(__file__).resolve().parents[1]
    (root/"results/algebra-audit.json").write_text(
        json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
