# The optimal moving sofa as r = f(α), t = g(α)

**Read the article: https://nwzd.github.io/moving-sofa/**

The Math.SE question [“Modelling the moving sofa”](https://math.stackexchange.com/questions/1787466/modelling-the-moving-sofa) builds a sofa by turning an L-shaped hallway of width 1 around it, with the hallway's inner corner at

```math
A(\alpha) = \bigl(r\cos\alpha,\; t\sin\alpha\bigr), \qquad 0 \le \alpha \le \pi,
```

after the hallway has turned through α/2. Constant r and t give ellipses; the best one has area 2.217857. The question asked whether r = f(α), t = g(α) could describe the optimal sofa.

It can. This repository gives f and g in closed form for Gerver's sofa (area 2.21953166887…, proven optimal by [Baek, 2024](https://arxiv.org/abs/2411.19826)), checks them, and animates the motion. The shape is Gerver's sofa itself, not a new one; what is new is its description in the question's coordinates.

## The result

f and g are mirror-symmetric about α = π/2 and have three pieces on [0, π/2], split at α = 2φ ≈ 4.49° and α = 2θ ≈ 78.07°. With s = α/2:

| Piece | f(α) | g(α) |
| --- | --- | --- |
| 0 ≤ α ≤ 2φ | a + (K₁ − cos s + ½ sin s − ¼ sin α) / cos α | a − 1/(2 cos s) − ¼ tan(α/4) |
| 2φ ≤ α ≤ 2θ | (P cos s + Q sin s + K₂) / cos α | (P sin s − Q cos s + L₂) / sin α |
| 2θ ≤ α ≤ π/2 | [(c − π/4)/sin(π/4 + s) + (π/4 − s)/sin(π/4 − s)] / √2 | [√2(c − π/4) sin(π/4 + s) − √2(π/4 − s) sin(π/4 − s) + L₃] / sin α |

where P = b₂ + b₁s − s²/4 and Q = 1 + b₁ − s/2, and

| φ | θ | a | K₁ | b₁ | b₂ | K₂ | L₂ | c | L₃ |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 0.0391773648 | 0.6813015094 | 1.2103224221 | 0.4034408074 | −0.5276245980 | 0.9202583852 | −0.3054160633 | 0.4724066198 | 0.6260455228 | 0.8896264790 |

Full-precision constants, the derivation from [Romik's exact solution (2018)](https://arxiv.org/abs/1606.08111), and the checks are in the article.

## Checks

| Check | Result | Reference |
| --- | --- | --- |
| Closed forms vs. the 40-digit solver (2,001 angles) | max difference 8 × 10⁻¹⁶ | 0 |
| Area recomputed from f and g | 2.2195316686 | 2.2195316689 (Gerver) |
| Same area method, circle r = t = 2/π | 2.20741616 | 2.20741610 (exact) |
| Same area method, best ellipse | 2.217856994 | 2.217856998 |

## Code

```
pip install -r code/requirements.txt
cd code
python gerver_path.py    # solve Romik's 22 equations for the constants (40 digits), print f and g
python verify.py         # closed forms vs. solver; also a quick fixed-angle area estimate, which runs high
python area_precise.py   # sofa area from f and g, about 3 minutes
```

- `code/gerver_path.py` solves for the 22 constants of Gerver's rotation path and maps it into (r, t) form.
- `code/verify.py` holds the closed-form f and g and compares them with the solver.
- `code/area_precise.py` computes the sofa area directly from f and g: for each vertical slice it optimizes over the hallway angle, then integrates by Gauss–Legendre quadrature.
- `docs/` is the article (GitHub Pages). Its figures and animation are computed in the browser from the same closed forms (`docs/sofa.js`).

## License

Code: MIT (see `LICENSE`). Article text and figures in `docs/`: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
