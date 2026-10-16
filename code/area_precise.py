"""Sofa area from a corner path A(alpha), without the angle-sampling bias of verify.sofa_area:
for every slice x, the wall envelopes max_a lower(x, a) and min_a upper(x, a) are found on a grid
and then refined by golden-section search, and the slice heights are integrated by Gauss-Legendre."""
import sys
import numpy as np
from verify import corner, f_closed, g_closed

GR = (np.sqrt(5) - 1) / 2


def walls(x, al, cornerfun):
    Ax, Ay = cornerfun(al)
    t = al / 2
    ct, st = np.cos(t), np.sin(t)
    dx = x - Ax
    q01 = ct * dx - st * Ay
    q02 = -st * dx - ct * Ay
    with np.errstate(divide='ignore', invalid='ignore'):
        lower = np.minimum(-q01 / st, -q02 / ct)
        upper = np.minimum((1 - q01) / st, (1 - q02) / ct)
    return lower, upper


def refine(x, a0, cornerfun, which, sign, halfwidth, iters=60):
    """Extremum of a wall over [a0-h, a0+h] by golden-section search (vectorized over x):
    sign=+1 returns max of the wall, sign=-1 returns min."""
    lo = np.clip(a0 - halfwidth, 1e-12, np.pi - 1e-12)
    hi = np.clip(a0 + halfwidth, 1e-12, np.pi - 1e-12)
    F = lambda a: sign * walls(x, a, cornerfun)[which]
    for _ in range(iters):
        c, d = hi - GR * (hi - lo), lo + GR * (hi - lo)
        m = F(c) > F(d)
        hi, lo = np.where(m, d, hi), np.where(m, lo, c)
    return sign * F((lo + hi) / 2)


def heights(xs, cornerfun, n_grid=24000):
    grid = (np.arange(n_grid) + 0.5) * np.pi / n_grid
    h = np.pi / n_grid
    out = np.empty_like(xs)
    for i in range(0, len(xs), 256):
        x = xs[i:i + 256, None]
        L, U = walls(x, grid[None, :], cornerfun)
        aL = grid[np.argmax(L, axis=1)][:, None]
        aU = grid[np.argmin(U, axis=1)][:, None]
        Lmax = np.maximum(refine(x, aL, cornerfun, 0, +1, 1.5 * h)[:, 0], L.max(axis=1))
        Umin = np.minimum(refine(x, aU, cornerfun, 1, -1, 1.5 * h)[:, 0], U.min(axis=1))
        out[i:i + 256] = np.clip(np.minimum(Umin, 1) - np.maximum(Lmax, 0), 0, None)
    return out


def area(cornerfun, n_sub, order=8):
    x_end = 1 + cornerfun(np.array([0.0]))[0][0]
    xg, wg = np.polynomial.legendre.leggauss(order)
    edges = np.linspace(0, x_end, n_sub + 1)
    mids, half = (edges[1:] + edges[:-1]) / 2, (edges[1:] - edges[:-1]) / 2
    xs = (mids[:, None] + half[:, None] * xg[None, :]).ravel()
    ws = (half[:, None] * wg[None, :]).ravel()
    return 2 * np.sum(ws * heights(xs, cornerfun))      # sofa is symmetric in x


if __name__ == '__main__':
    gerver = lambda a: (f_closed(a) * np.cos(a), g_closed(a) * np.sin(a))
    ell = lambda r, tt: (lambda a: (r * np.cos(a), tt * np.sin(a)))
    for n in (500, 1000, 2000):
        print(f'n_sub={n:5d}  Gerver(f,g)={area(gerver, n):.12f}  '
              f'Hammersley={area(ell(2/np.pi, 2/np.pi), n):.12f}  '
              f'best ellipse={area(ell(0.605513519698965, 0.6678342468712839), n):.12f}', flush=True)
    print('exact: Gerver 2.219531668871973  Hammersley', 2 / np.pi + np.pi / 2, ' ellipse 2.217856997942074')
