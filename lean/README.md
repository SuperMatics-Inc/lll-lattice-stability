# Lean verification of the GST lattice certificate argument

Toolchain: Lean 4.35.0, Mathlib v4.35.0, `autoImplicit=false`.

The library checks algebra, uniform near-origin scaling, exact rational
rounding and the geometry of two adaptive subdivisions. The
spectral-determinant theorems assume the analytic enclosures listed below.

## Build and axiom dependencies

Run `lake build`. The default target imports `GSTStability.AxiomCheck`, which
lists each main theorem's axiom dependencies. The expected list contains only
`propext`, `Classical.choice` and `Quot.sound`. The proofs use no custom axiom,
`sorry` or `native_decide`. The rational trace checks use `decide +kernel`.

For a portable checkout use this dependency in `lakefile.toml`:

```toml
[[require]]
name = "mathlib"
scope = "leanprover-community"
rev = "v4.35.0"
```

Fetch the precompiled Mathlib files with `lake exe cache get` before `lake build`. GitHub Actions builds this dependency configuration.

## Proof map

- `Algebra.lean`: the real-coordinate expansion from the symmetric moments
  proves `moment_quadratic_cancellation`; finite Cauchy products prove
  `cancellation_degree_le_three` and `quartic_identity` (Proposition 3.1).
  These identities hold for arbitrary real parameters and arbitrary directions.
  `omitted_factor_count` proves the degree-20/eleven-factor combinatorial step.
- `Remainder.lean`: nonnegative polynomial majorants are monotone; the
  exponential envelope is nonnegative and monotone; the degree-eleven
  exponential remainder and the product error factor respectively as rho^11
  and rho^4 times nonnegative monotone expressions. `near_lower_bound` transfers
  an exact margin at r to a uniform lower enclosure; `near_positive` proves
  positivity for **every** `0 < rho <= r`, not only on an annulus.
- `Coverage.lean`: four closed children cover each parent, the eight initial
  squares cover the closed half-torus, and structural induction proves
  `adaptive_exhaustion` for any terminating adaptive subdivision. A skipped
  cell is proved to lie entirely within the near square.
- `Trace.lean`: an exact-rational checker proves that every omitted near leaf
  of a recorded tree meets the omission rule. `far_lower_bounds_positive`
  derives strict positivity from lower enclosures and a positive endpoint;
  `checked_trace_sound` converts a checked tree into a geometric cover.
- `Recorded.lean`: actual subdivision shapes generated from
  `data/trace-minus.json` and `data/trace-plus.json` at the repository root. Each tree has the same
  dyadic centers and half-widths as its trace; the geometry of **every** near
  leaf is kernel-checked, with no arbitrary coverage hypothesis. The final
  `minus_punctured_torus_positive` and `plus_punctured_torus_positive` instantiate
  the torus assembly at radii `1/100` and `1/40`, deriving positivity from the
  printed rational near/far endpoints. Exact kernel-checked (near, far) leaf
  counts are `(20, 756)` and `(12, 707)`.
- `Main.lean`: a Taylor lower-bound rule, exact rational positivity and
  downward-rounding checks for all six printed P4/near/far bounds, reflection
  across the origin, exact integer centering, and positivity on R^2 minus Z^2.
  For a periodic symbol the latter is precisely positivity on the punctured
  torus, including fundamental-domain boundaries.

## Explicit analytic hypotheses and unformalized material

1. The Gaussian theta sums have the stated stationary value, real/even jets,
   moment symmetries, and termwise derivatives. The algebraic consequences
   are checked; normal convergence, Poisson summation, and the strip-basis
   overlap are not formalized.
2. The interval enclosures establish the lower bound for the homogeneous P4
   on the unit boundary, all higher coefficient majorants, and the Taylor and
   Gaussian-shell remainder bounds. `near_lower_bound` takes the actual Taylor
   inequality, monotone envelope bounds, and an endpoint margin as arguments.
   The final recorded theorems require the resulting continuous enclosure
   `near_endpoint * radius^4 <= Sigma`, not near positivity as a hypothesis.
   The near endpoints are exactly `33019570/100000000` and `384161/100000000`.
3. Every far leaf requires the continuous-cell lower enclosure in
   `QuadTree.farLowerBounds`, with common rational endpoint
   `104595/1000000000000000` or `424467/100000000000000`.
   The final assembly proves that these endpoints imply strict positivity.
4. The function’s evenness and integer periodicity are explicit hypotheses.
   Lean checks arithmetic on the printed enclosures; the interval calculation
   supplies the enclosures themselves.
5. The Lean conclusion is determinant positivity, conditional on those analytic
   inputs. It does not formalize fibre eigenvalues, the zero Jordan block,
   modular non-equivalence, persistence in tau, the hexagonal shape theorem,
   weighted/dispersive estimates, or nonlinear stability.

The interval verifier and the accompanying report supply the analytic arguments outside Lean.
