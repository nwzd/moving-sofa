"""Gerver's sofa from just two numbers, phi and theta.

All other constants follow from phi and theta by explicit formulas (smoothness of the
corner path at tau = phi and tau = theta), and phi, theta themselves solve two equations
(the contact condition). This script checks both against the 40-digit solver in
gerver_path.py.
"""
import mpmath as mp

mp.mp.dps = 40
pi = mp.pi


def constants(phi, theta):
    """The eight path constants as explicit functions of phi and theta."""
    sf, cf, st, ct = mp.sin(phi), mp.cos(phi), mp.sin(theta), mp.cos(theta)
    c = ((mp.mpf(1) / 2 + (2 - phi - pi / 4 * (theta - phi) + (theta**2 - phi**2) / 4) * sf
          - (phi + mp.mpf(3) / 2 - pi / 2) * cf) / (cf - (1 + (theta - phi) / 2) * sf))
    a = (c + phi + mp.mpf(3) / 2 - pi / 2 - cf / 2) / (2 * sf)
    b1 = (pi / 2 - 2 - c) / 2
    b2 = c - theta + mp.mpf(1) / 2 - b1 * theta + theta**2 / 4
    h = (c + theta - pi / 2) / 2
    m = a * cf - sf / 4 - mp.mpf(1) / 2
    K2 = -ct / 2 - h * st
    K1 = K2 + m * cf + sf / 4
    L2 = mp.mpf(1) / 4 + cf / 4 - m * sf
    L3 = L2 + st / 2 - h * ct
    return dict(a=a, K1=K1, b1=b1, b2=b2, K2=K2, L2=L2, c=c, L3=L3)


def contact(phi, theta):
    """Zero exactly when the inner corner first touches the sofa where the inner wall last does."""
    k = constants(phi, theta)
    h = (k['c'] + theta - pi / 2) / 2
    X1 = k['a'] * mp.cos(2 * phi) - mp.sin(2 * phi) / 4 - mp.cos(phi) + mp.sin(phi) / 2 + k['K1']
    Y1 = k['a'] * mp.sin(2 * phi) + mp.cos(2 * phi) / 4 - mp.sin(phi) - mp.cos(phi) / 2 + mp.mpf(1) / 4
    return [X1 - (2 * h * mp.sin(theta) + mp.cos(theta)),
            Y1 - (2 * h * mp.cos(theta) - mp.sin(theta) + k['L3'])]


if __name__ == '__main__':
    phi, theta = mp.findroot(lambda p, t: contact(p, t), (mp.mpf('0.04'), mp.mpf('0.7')))
    print('phi   =', mp.nstr(phi, 25))
    print('theta =', mp.nstr(theta, 25))
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        import gerver_path as gp
    ref = dict(a=gp.C[gp.a1], K1=1 - gp.C[gp.a1] - gp.C[gp.k31], b1=gp.C[gp.b1], b2=gp.C[gp.b2],
               K2=gp.C[gp.k21] - gp.C[gp.k31], L2=gp.C[gp.k22], c=gp.C[gp.c1], L3=gp.C[gp.k32])
    print('phi, theta vs solver:', mp.nstr(abs(phi - gp.C[gp.phi]), 3), mp.nstr(abs(theta - gp.C[gp.th]), 3))
    for name, val in constants(phi, theta).items():
        print(f'{name:3s} = {mp.nstr(val, 20):>24s}   |diff vs solver| = {mp.nstr(abs(val - ref[name]), 3)}')
