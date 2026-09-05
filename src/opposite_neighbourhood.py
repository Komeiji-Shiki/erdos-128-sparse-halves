"""Cross-cut neighbourhood bounds with separate exploratory utilities.

Only interval_saving is used by the rational proof checker. The numerical
optimizer and sampled LP are discovery tools, not certificates.
"""


def saving(a, mu, rho, cut, alpha=9/25):
    c = 1-a
    h = 2*c-1
    top = 2*alpha*c
    cap = min(top,h)
    cross = 4*c*c*(rho/2-cut)
    if h <= 0 or cap-cross+2*mu <= 0:
        return 0
    mean_s = (4*top*cap*mu*mu/(h*(cap-cross+2*mu))
              +mu*(1-top)*cross/h-top*mu)
    if mu <= top*top/2:
        mean_s += 2*mu*mu*(1-cross/(h*top))
    return min(4*a*a,a*c)*max(0,mean_s)


def interval_saving(a_interval, mu_interval, rho_lower, cut, alpha):
    """Rational lower bound, enabled only under explicit sufficient conditions."""
    al, ah = a_interval
    ml, mh = mu_interval
    if 4*ah > 1 or rho_lower < 2*cut or ml <= 0:
        return 0
    cl, ch = 1-ah, 1-al
    tl, th = 2*alpha*cl, 2*alpha*ch
    hl, hh = 2*cl-1, 2*ch-1
    kl, kh = min(tl,hl), min(th,hh)
    cross_density = rho_lower/2-cut
    zl, zh = 4*cl*cl*cross_density, 4*ch*ch*cross_density
    coefficient = ml*(1-th)-2*mh*mh/tl
    if mh > tl*tl/2 or coefficient < 0 or zh > kl:
        return 0
    lower_s = (4*tl*kl*ml*ml/(hh*(kh-zl+2*mh))
               +coefficient*zl/hh-mh*th+2*ml*ml)
    radius = min(4*al*al,al*(1-al))
    return radius*max(0,lower_s)/(4*ch*ch)


def sampled_affine_bound(a, mu, rho, cut, alpha=9/25):
    """Exploratory LP over sampled moment constraints, never a certificate."""
    import numpy as np
    from scipy.optimize import linprog
    c = 1-a
    h = 2*c-1
    top = 2*alpha*c
    cap = min(top,h)
    cross = 4*c*c*(rho/2-cut)
    mean_d = cross/h
    mean_e = 4*cap*mu*mu/(h*(cap-cross+2*mu))
    constraints, values = [], []
    for d in np.linspace(0,top,2001):
        lower = max(0,mu-(1-d)**2/4)
        upper = min(mu,top*d/2)
        if lower > upper:
            continue
        for e in (lower,upper):
            constraints.append([e,d,1])
            values.append(d*e-mu*d*d)
    fit = linprog([-mean_e,-mean_d,-1],A_ub=constraints,b_ub=values,
                  bounds=[(0,None),(None,None),(None,None)],method="highs")
    return {"plane":fit.x.tolist(),"sampled_lower":-fit.fun,
            "mean_d":mean_d,"mean_e_lower":mean_e,"degree_cap":top}


if __name__ == "__main__":
    import sys
    from scipy.optimize import differential_evolution
    from sparse_set_bound import phi
    from explore_tangent import tangent_saving
    def upper(x, rho, lower):
        a, fraction = x
        c = 1-a
        mu = fraction*8*c*c/47
        base = min(mu/4, a*rho/2+c*2/47-(1+a)*mu/(4*c))
        gain = max(tangent_saving(a,mu,.95)[0], saving(a,mu,lower,2/47))
        result = base-gain/(4*c*c)
        if a <= .25:
            result = min(result,2/47-mu/(4*c*c)+float(phi(1-1/(2*c),9/25)))
        return result
    for rho in [.25,.275,.29,.3,.305,.31,.315,.3197]:
        lower = rho
        fit = differential_evolution(lambda x:-upper(x,rho,lower),[(.001,.49),(.001,1)],
                                     seed=128,tol=1e-9,popsize=8)
        print(rho,lower,fit.x.tolist(),-fit.fun)
