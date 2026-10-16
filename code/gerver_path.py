"""Gerver's (optimal) sofa as r = f(alpha), t = g(alpha) in the parametrization of
https://math.stackexchange.com/q/1787466 :  inner corner A(alpha) = (f cos(alpha), g sin(alpha)),
hallway rotated by alpha/2, alpha in [0, pi].

The rotation path comes from D. Romik, "Differential equations and exact solutions in the
moving sofa problem" (Exp. Math. 2018), Section 6: five ODE phases glued C^1, solved here
for the 22 constants with mpmath.
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 40

t = sp.symbols('t', real=True)
(phi, th, a1, a2, b1, b2, c1, c2, d1, d2, e1, e2,
 k11, k12, k21, k22, k31, k32, k41, k42, k51, k52) = U = sp.symbols(
    'phi th a1 a2 b1 b2 c1 c2 d1 d2 e1 e2 k11 k12 k21 k22 k31 k32 k41 k42 k51 k52', real=True)


def R(v):
    return sp.Matrix([sp.cos(t) * v[0] - sp.sin(t) * v[1], sp.sin(t) * v[0] + sp.cos(t) * v[1]])


Q = sp.Rational(1, 4)
H = sp.Rational(1, 2)
X1 = R([a1 * sp.cos(t) + a2 * sp.sin(t) - 1, -a2 * sp.cos(t) + a1 * sp.sin(t) - H]) + sp.Matrix([k11, k12])
X2 = R([-Q * t**2 + b1 * t + b2, H * t - b1 - 1]) + sp.Matrix([k21, k22])
X3 = R([c1 - t, c2 + t]) + sp.Matrix([k31, k32])
X4 = R([-H * t + d1 - 1, -Q * t**2 + d1 * t + d2]) + sp.Matrix([k41, k42])
X5 = R([e1 * sp.cos(t) + e2 * sp.sin(t) - H, -e2 * sp.cos(t) + e1 * sp.sin(t) - 1]) + sp.Matrix([k51, k52])
D = lambda X: X.diff(t)
e_t = sp.Matrix([sp.cos(t), sp.sin(t)])
f_t = sp.Matrix([-sp.sin(t), sp.cos(t)])
A = lambda X: X + D(X).dot(e_t) * f_t + e_t        # outer-wall contact point
B = lambda X: X + D(X).dot(e_t) * f_t              # inner-wall contact point

P2 = sp.pi / 2
eqs = [e1 - a1, e2 + a2, d1 - (sp.pi / 4 - b1), d2 - (b2 + sp.pi / 4 * (2 * b1 - sp.pi / 4)), c2 - (c1 - P2)]
eqs += list(X1.subs(t, 0)) + [A(X1)[1].subs(t, 0)]    # x(0)=0, A(0)=(1,0); A_1(0)=1 is then automatic
for Xa, Xb, s, deriv in [(X1, X2, phi, True), (X2, X3, th, True), (X3, X4, P2 - th, False), (X4, X5, P2 - phi, False)]:
    eqs += list((Xa - Xb).subs(t, s))
    if deriv:
        eqs += list(D(Xa - Xb).subs(t, s))
eqs += list(X1.subs(t, phi) - B(X3).subs(t, P2 - th))
assert len(eqs) == 22

F = sp.lambdify(U, eqs, 'mpmath')
guess = [0.039177364790083641, 0.681301509382724894, 1.210322422072688751, -0.25,
         -0.527624598026784624, 0.920258385160637622, 0.626045522848465867, -0.944750803946430751,
         1.313022761424232933, -0.525382670414554437, 1.210322422072688751, 0.25,
         -0.210322422072688751, 0.25, -0.919179292771593322, 0.472406619750805465,
         -0.613763229430251668, 0.889626479003221860, -0.308347166088910014, 0.472406619750805465,
         -1.017204036787814585, 0.25]
sol = mp.findroot(lambda *u: F(*u), [mp.mpf(g) for g in guess], tol=mp.mpf(10)**-35)
C = {s: v for s, v in zip(U, sol)}
SUB = {s: sp.Float(str(v), 40) for s, v in C.items()}

# Asker's coordinates: alpha = 2t, shift x so the path is symmetric about x = 0.
alpha = sp.symbols('alpha', real=True)
shift = sp.N(X3[0].subs(SUB).subs(t, sp.pi / 4), 40)   # equals k31 because c2 = c1 - pi/2
pieces_t = [(X1, 0, C[phi]), (X2, C[phi], C[th]), (X3, C[th], mp.pi / 2 - C[th]),
            (X4, mp.pi / 2 - C[th], mp.pi / 2 - C[phi]), (X5, mp.pi / 2 - C[phi], mp.pi / 2)]
pieces = [(sp.lambdify(t, (X[0].subs(SUB) - shift), 'mpmath'), sp.lambdify(t, X[1].subs(SUB), 'mpmath'),
           2 * lo, 2 * hi) for X, lo, hi in pieces_t]


def A_of_alpha(al):
    """Inner corner A(alpha) = (X, Y) of the hallway, alpha in [0, pi]."""
    al = mp.mpf(al)
    for fx, fy, lo, hi in pieces:
        if al <= hi:
            return fx(al / 2), fy(al / 2)
    raise ValueError(al)


def f(al):
    al = mp.mpf(al)
    if abs(al - mp.pi / 2) < mp.mpf(10)**-12:
        return -mp.diff(lambda s: A_of_alpha(s)[0], mp.pi / 2)   # L'Hopital
    return A_of_alpha(al)[0] / mp.cos(al)


def g(al):
    al = mp.mpf(al)
    if al < mp.mpf(10)**-12:
        return C[a1] - mp.mpf(1) / 2                                # = Y'(0)
    if mp.pi - al < mp.mpf(10)**-12:
        return C[a1] - mp.mpf(1) / 2
    return A_of_alpha(al)[1] / mp.sin(al)


if __name__ == '__main__':
    print('residual', max(abs(v) for v in F(*sol)))
    for s in U:
        print(f'{str(s):4s} = {mp.nstr(C[s], 25)}')
    print('shift (x of path midpoint) =', mp.nstr(mp.mpf(str(shift)), 25))
    print('alpha breakpoints:', [mp.nstr(p[2], 18) for p in pieces], mp.nstr(pieces[-1][3], 18))
    for al in [0, 0.05, 0.0783547, 0.1, 0.5, 1, 1.3626, 1.5, mp.pi / 2]:
        print(f'alpha={mp.nstr(al, 8):>10s}  f={mp.nstr(f(al), 15)}  g={mp.nstr(g(al), 15)}')
