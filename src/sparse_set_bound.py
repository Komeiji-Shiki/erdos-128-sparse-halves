"""A robust sparse-set extension bound, plus numerical exploration.

The proof of phi belongs in the manuscript. Floating-point optimization here
only proposes constants for the separate rational certificate checker.
"""

import json
from fractions import Fraction as F
from pathlib import Path

def independent_bound(z):
    if z >= F(3,8):
        return z*(F(1,2)-z)/2
    return F(9,128)-z/8


def phi(z, alpha=F(3,8)):
    if z < F(1,3) or z > F(1,2):
        raise ValueError("Sparse-set mass must lie between one third and one half")
    if not F(1,3) <= alpha <= F(3,8):
        raise ValueError("This form uses an independence cap at most three eighths")
    second = ((12*alpha-32*z*z-12*z+9)/192 if z <= F(3,8)
              else (alpha-8*z*z+3*z)/16)
    return max(independent_bound(z), second)


def bound(a, mu, rho=3197/10000, cut=2/47, alpha=3/8):
    from perturbation_bound import savings
    c = 1-a
    p = a*rho/2+c*cut
    baseline = min(mu/4, p-(1+a)*mu/(4*c))
    perturbed = baseline-min(savings(a,mu))/(4*c*c)
    if a <= .25:
        z = 1-1/(2*c)
        extension = cut-mu/(4*c*c)+float(phi(z, alpha))
        return min(perturbed, extension)
    return perturbed


def main():
    from scipy.optimize import differential_evolution
    fit = differential_evolution(lambda x:-bound(x[0], x[1]*8*(1-x[0])**2/47),
                                 [(0,.5),(0,1)], seed=128, tol=1e-11, polish=True)
    a, fraction = map(float,fit.x)
    mu = fraction*8*(1-a)**2/47
    result = {"status":"numerical exploration only", "a":a,"mu":mu,
              "mass_cap_fraction":fraction,"proposed_upper":float(-fit.fun),
              "small_part_mass":1-1/(2*(1-a)),
              "small_part_edge_upper":2/47-mu/(4*(1-a)**2)}
    Path("results/sparse-set-numerical.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
