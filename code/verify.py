"""Check the closed-form f(alpha), g(alpha) against the solver, then compute the sofa area
directly from them with the asker's construction (rotate the hallway by alpha/2 about the
inner corner A = (f cos a, g sin a), keep what stays inside every hallway)."""
import numpy as np

# Constants (Gerver's sofa, from Romik's ODE solution; see gerver_path.py)
PHI = 0.039177364790083641863      # rotation angle where the inner corner starts touching
THETA = 0.68130150938272489447     # rotation angle where the 'pure corner' phase starts
A1 = 1.2103224220726887514
B1 = -0.52762459802678462416
B2 = 0.92025838516063762289
C1 = 0.62604552284846586755
K1 = 1 - A1 + 0.61376322943025166855          # 1 - a1 - k31
K2X = -0.91917929277159332227 + 0.61376322943025166855   # k21 - k31
K2Y = 0.47240661975080546518                  # k22
K3Y = 0.88962647900322186073                  # k32


def corner(al):
    """Inner corner A(alpha) = (X, Y) for alpha in [0, pi] (symmetric: X(pi-a)=-X(a), Y(pi-a)=Y(a))."""
    al = np.asarray(al, float)
    m = al > np.pi / 2
    a = np.where(m, np.pi - al, al)
    s = a / 2
    # phase 1: 0 <= a <= 2*phi
    X1 = A1 * np.cos(a) - np.sin(a) / 4 - np.cos(s) + np.sin(s) / 2 + K1
    Y1 = A1 * np.sin(a) + (1 + np.cos(a)) / 4 - np.sin(s) - np.cos(s) / 2
    # phase 2: 2*phi <= a <= 2*theta
    P = B2 + B1 * s - s * s / 4
    Qv = 1 + B1 - s / 2
    X2 = P * np.cos(s) + Qv * np.sin(s) + K2X
    Y2 = P * np.sin(s) - Qv * np.cos(s) + K2Y
    # phase 3: 2*theta <= a <= pi/2
    p, w = C1 - np.pi / 4, np.pi / 4 - s
    X3 = p * (np.cos(s) - np.sin(s)) + w * (np.cos(s) + np.sin(s))
    Y3 = p * (np.cos(s) + np.sin(s)) - w * (np.cos(s) - np.sin(s)) + K3Y
    X = np.where(a <= 2 * PHI, X1, np.where(a <= 2 * THETA, X2, X3))
    Y = np.where(a <= 2 * PHI, Y1, np.where(a <= 2 * THETA, Y2, Y3))
    return np.where(m, -X, X), Y


def f_closed(al):
    al = np.asarray(al, float)
    a = np.minimum(al, np.pi - al)
    s = a / 2
    f1 = A1 + (-np.sin(a) / 4 - np.cos(s) + np.sin(s) / 2 + K1) / np.cos(a)
    f2 = corner(a)[0] / np.cos(a)
    with np.errstate(divide='ignore', invalid='ignore'):
        f3 = ((C1 - np.pi / 4) / np.sin(np.pi / 4 + s) +
              np.where(np.abs(np.pi / 4 - s) < 1e-12, 1.0, (np.pi / 4 - s) / np.sin(np.pi / 4 - s))) / np.sqrt(2)
    return np.where(a <= 2 * PHI, f1, np.where(a <= 2 * THETA, f2, f3))


def g_closed(al):
    al = np.asarray(al, float)
    a = np.minimum(al, np.pi - al)
    s = a / 2
    g1 = A1 - 1 / (2 * np.cos(s)) - np.tan(a / 4) / 4
    with np.errstate(divide='ignore', invalid='ignore'):
        g23 = corner(a)[1] / np.sin(a)
    return np.where(a <= 2 * PHI, g1, g23)


def sofa_area(cornerfun, n_alpha=20000, n_x=20000):
    al = (np.arange(n_alpha) + 0.5) * np.pi / n_alpha       # interior angles
    Ax, Ay = cornerfun(al)
    t = al / 2
    ct, st = np.cos(t), np.sin(t)
    x0 = cornerfun(np.array([0.0]))[0][0]                    # = f(0)
    xs = -(1 + x0) + (np.arange(n_x) + 0.5) * (2 * (1 + x0)) / n_x
    total = 0.0
    for chunk in np.array_split(xs, max(1, n_x // 500)):
        dx = chunk[:, None] - Ax[None, :]
        q01 = ct * dx - st * Ay          # q1 = q01 + y sin t   (hallway coords, rotated back)
        q02 = -st * dx - ct * Ay         # q2 = q02 + y cos t
        upper = np.minimum((1 - q01) / st, (1 - q02) / ct).min(axis=1)        # outer walls
        lower = np.minimum(-q01 / st, -q02 / ct).max(axis=1)                  # inner walls
        total += np.clip(np.minimum(upper, 1) - np.maximum(lower, 0), 0, None).sum()
    return total * 2 * (1 + x0) / n_x


if __name__ == '__main__':
    import mpmath as mp
    import gerver_path as gp
    grid = np.linspace(0, np.pi, 2001)
    ref_f = np.array([float(gp.f(a)) for a in grid])
    ref_g = np.array([float(gp.g(a)) for a in grid])
    print('max |f_closed - f_solver| =', np.max(np.abs(f_closed(grid) - ref_f)))
    print('max |g_closed - g_solver| =', np.max(np.abs(g_closed(grid) - ref_g)))
    ell = lambda r, tt: (lambda a: (r * np.cos(a), tt * np.sin(a)))
    print('ellipse r=.605513519698965 t=.6678342468712839:',
          sofa_area(ell(0.605513519698965, 0.6678342468712839), 8000, 8000), '(thread: 2.217856997942)')
    print('Hammersley r=t=2/pi:', sofa_area(ell(2 / np.pi, 2 / np.pi), 8000, 8000),
          '(exact', 2 / np.pi + np.pi / 2, ')')
    for n in (4000, 8000, 16000):
        print(f'Gerver via f,g  n={n}:', sofa_area(lambda a: (f_closed(a) * np.cos(a), g_closed(a) * np.sin(a)), n, n))
    print('Gerver constant: 2.21953166887197...')
