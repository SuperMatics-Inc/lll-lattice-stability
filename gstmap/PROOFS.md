# PROOFS: spectral stability map of the LLL Abrikosov lattice (GST Hamiltonian setting)

Setting (Germain--Schwinte--Thomann, arXiv:2404.06085; notation of the associated article):
`tau = t1 + i t2 in H`, `c = sqrt(t2)/pi`,

    G_tau(x,y) = c exp( -pi t2 (x^2+y^2) + 2 pi i t1 x y ),      lambda_tau = sum_{m,n} G_tau(m,n),
    K(s,xi)  = sum_{m,n} G(m+s, n) e^{2 pi i xi n}             (real),
    B(s,xi)  = sum_{m,n} G(m+s, n-s) e^{2 pi i xi (n-m)}       (complex),
    A = 2K - lambda,      Sigma_tau(s,xi) = A^2 - |B|^2.

The fibre spectrum over Bloch momentum `(s,xi)` is `{ +-sqrt(Sigma) }` (paper, Thm. on the fibre operator; at `Sigma=0`
a nilpotent block). **Spectral stability of tau** means `Sigma_tau >= 0` on the whole Bloch torus.
All sums converge absolutely (Gaussian), so every identity below is an identity of real-analytic functions.

Conventions: PROVED denotes the analytic arguments below, CERTIFIED an interval computation with the stated semantics,
LEAN kernel-checked algebra/coverage with explicit analytic hypotheses, and NUMERICAL a non-rigorous calculation.

## 0. Well-definedness on the torus; energy form

`K(s+1,xi)=K(s,xi)=K(s,xi+1)` (shift `m`, resp. `e^{2 pi i n}=1`), and `B(s+1,xi)=e^{4 pi i xi}B(s,xi)`,
`B(s,xi+1)=B(s,xi)`. Hence `|B|`, `A`, `Sigma` are 1-periodic in both variables: `Sigma` lives on the torus `R^2/Z^2`.

**Energy equivalence (PROVED).** The second variation of the Hamiltonian energy on Bloch perturbations is the Hermitian form
`W(s,xi) = [[A, B],[conj B, A]]` (paper, Sec. 2). `W >= 0` iff `A >= |B|`. At `(0,0)`, `A = lambda > 0`.
If `Sigma > 0` on the torus minus `{0}`, then `A != 0` there, the punctured torus is connected, so `A>0` throughout and
`W >= 0`; its two eigenvalues are `A +/- |B|`, strictly positive away from zero. This is a fibrewise quadratic-energy statement, not a nonlinear local-minimum theorem. Conversely `W>=0` gives `Sigma>=0`.

## 1. Modular covariance (PROVED)

**Proposition 1.** For all `tau in H`, `(s,xi) in R^2`:

    (T) Sigma_{tau+1}(s,xi)   = Sigma_tau(s, xi+s),
    (R) Sigma_{-conj tau}(s,xi) = Sigma_tau(s, -xi),
    (S) Sigma_{-1/tau}(s,xi)  = Sigma_tau(xi, -s).

Each right-hand side is the composition of `Sigma_tau` with a unimodular linear bijection of the Bloch plane, so
`{Sigma_tau >= 0}` holds for one member of an orbit of `<T,S,R>` (the extended modular group `PGL_2(Z)`) iff for all.

*Proof of (T).* `G_{tau+1}(x,y) = G_tau(x,y) e^{2 pi i x y}`. In `K_{tau+1}(s,xi)` the factor is `e^{2 pi i (m+s) n}
= e^{2 pi i s n}` (since `e^{2 pi i m n}=1`), so `K_{tau+1}(s,xi) = K_tau(s, xi+s)`, `lambda_{tau+1} = lambda_tau`. In `B`:
`(m+s)(n-s) = mn + s(n-m) - s^2`, so `B_{tau+1}(s,xi) = e^{-2 pi i s^2} B_tau(s, xi+s)`. Thus `A_{tau+1}(s,xi)=A_tau(s,xi+s)`
and `|B|` likewise. QED.

*Proof of (R).* For real `x,y`, `G_{-conj tau}(x,y) = conj G_tau(x,y)` and `G_tau(x,-y) = conj G_tau(x,y)`. Therefore `K_tau`
is real and `K_{-conj tau}(s,xi) = conj K_tau(s,-xi) = K_tau(s,-xi)`, `B_{-conj tau}(s,xi) = conj B_tau(s,-xi)`,
`lambda_{-conj tau} = lambda_tau`. QED.

*Proof of (S).* Write `G_tau(x,y) = c e^{-pi v^T Q v}`, `v=(x,y)`, `Q = [[t2, -i t1],[-i t1, t2]]`
(`Re Q > 0`). Then `det Q = t2^2 + t1^2 = |tau|^2` and `Q^{-1} = |tau|^{-2}[[t2, i t1],[i t1, t2]] = Q_{tau'}` with
`tau' = -1/tau = (-t1 + i t2)/|tau|^2`, i.e. `t2' = t2/|tau|^2`, `t1' = -t1/|tau|^2`. The Fourier transform
`f^(p) = int f(v) e^{-2 pi i p.v} dv` of `e^{-pi v^T Q v}` is `(det Q)^{-1/2} e^{-pi p^T Q^{-1} p}` (analytic continuation in Q).
Since `c/|tau| = sqrt(t2)/(pi |tau|) = sqrt(t2')/pi = c'`:

    G_tau^ = G_{tau'}.

Poisson summation in `Z^2` for `F(x,y) = G_tau(x+s,y) e^{2 pi i xi y}`: `F^(k,l) = e^{2 pi i k s} G_{tau'}(k, l-xi)`, hence
`K_tau(s,xi) = sum_{k,l} e^{2 pi i k s} G_{tau'}(k,l-xi)`. Using `G(x,y)=G(y,x)=G(-x,-y)` and relabelling
`(m,n)=(-l,-k)`: `K_tau(s,xi) = K_{tau'}(xi, -s)`. For `lambda`: `lambda_tau = sum G_tau^ = lambda_{tau'}`.
For `F(x,y) = G_tau(x+s,y-s) e^{2 pi i xi (y-x)}`: with `u=x+s, v=y-s`,
`F^(k,l) = e^{4 pi i xi s} e^{2 pi i (k-l) s} G_{tau'}(k+xi, l-xi)`, so
`B_tau(s,xi) = e^{4 pi i xi s} B_{tau'}(xi, -s)`. Therefore `A_tau(s,xi) = A_{tau'}(xi,-s)`, `|B_tau(s,xi)| = |B_{tau'}(xi,-s)|`,
`Sigma_tau(s,xi) = Sigma_{tau'}(xi,-s)`. Applying this with `tau'` in place of `tau` (`(tau')' = tau`) gives (S) as stated. QED.

**Corollary 1 (reduction to F+).** Stability of `tau` depends only on the `PGL_2(Z)`-orbit. A fundamental domain of `PGL_2(Z)`
on `H` is `F+ = { 0 <= Re tau <= 1/2, |tau| >= 1 }`. Formalised in `lean/GSTMap/Modular.lean`
(`stable_T`, `stable_S`, `stable_R`, `stable_smul`, `reduce_to_FPlus`) from Mathlib's `SL2Z_generators` and
`exists_smul_mem_fd`, with (T),(S),(R) as the only analytic hypotheses (`Covariant`).

## 2. Local structure at the Goldstone point and the long-wave mechanism (quartic form)

`Sigma` is even, `Sigma(0)=0` and the Taylor coefficients of degree <= 3 vanish identically:
`Sigma = P4_tau(theta) rho^4 + P6_tau(theta) rho^6 + O(rho^8)`, `(s,xi)=rho(cos theta, sin theta)`, with
`P4 = 2 lambda (A4 - Re B4) - (Im B2)^2` (`A2 = Re B2`; subscripts are homogeneous Taylor parts of A and of the suitably gauged B).

* If `min_theta P4_tau(theta) < 0` then `Sigma<0` for small `rho` in the minimising direction: **tau is unstable** by a
  long-wave (`omega ~ |q|^2`, Tkachenko-type) mode.
* If `min_theta P4_tau > 0` then `Sigma>0` near `0`.
* Define the **quartic-degeneracy set** `Gamma4` by `min_theta P4_tau(theta)=0`. This definition does not identify it with the full spectral boundary. Section 10 encloses one analytic branch and its positive sextic coefficient.


## 3. Interval-arithmetic certificates (CERTIFIED)

**Unstable certificate (box U).** A box `Omega = [c1-w1,c1+w1] x [c2-w2,c2+w2] subset H` is *unstable-certified* if a rational Bloch
point `q0=(s0,xi0)` exists with `sup_{tau in Omega} Sigma_tau(q0) < 0`. The engine (`unst.py`, `gsteng.py: NearData.sigma_upper_at`)
evaluates `K,B,lambda` at `q0` by Arb balls with `tau` as a ball, truncating the lattice at `|m|,|n| <= M=6` with the explicit Gaussian
tail `T(M)` added to the radii; the output upper bound is rigorous for all `tau in Omega` simultaneously.

**Stable certificate (box S).** Claim `Sigma_tau(q) > 0` for all `q in torus minus {0}`, all `tau in Omega`.
(a) NEAR `|q| <= rmax`: Taylor expansion of `Sigma_tau` about `0` of order `N` with all coefficients as Arb balls in `tau in Omega`
(`jets_full`); each homogeneous part is bounded on angular cells by `lo_j`, using a centred value plus interval derivatives in `tau`.
The terms above degree `N` are enclosed by the positive exponential/Gaussian coefficient majorant `taylor_tail_eps`.
For `rin<=rho<=rmax` the resulting Taylor lower bound is positive; for `0<rho<=rin` the quartic dominates.
The inner-disc bound does **not** treat `E(rin)` as a constant error down to `rho=0`: if `R` starts in degree `N+1`,
its coefficientwise absolute majorant gives `|R(rho)| <= (rho/rin)^(N+1) E_R(rin)`.
The omitted-lattice part has the same property: bound its positive majorant at `rin` and then discard its terms of degree at most `N`.
Products with the bounded low-degree `A,B` series preserve this scaling, so `E_Sigma(rho)/rho^4 <= E_Sigma(rin)/rin^4` for `N+1>=4`.
This is the justification for the inner-disc subtraction in `near_cover`, not an inference from a pointwise error bound alone.
(b) FAR `|q|>=rmax`: adaptive subdivision of the torus, on each cell `Sigma>= Sigma(center) - L h` with `L` from rigorous second-derivative
bounds `KAP6`; all cells must have a positive lower bound.
Both steps are in Arb balls with 192-bit precision. `boxcert.try_stable` returns `(ok, info)` only if every cell closes.

The mathematical content (a)+(b) is certified *per box*, and the boxes are glued by `lean/GSTMap/Tree.lean` (coverage of F+ cap {Im<=3}).

## 4. Cusp lemma (CERTIFIED, computer-assisted): every tau with Im tau >= 3 is unstable

**Lemma C.** For all `t1 in R` and `t2 >= 3`, `Sigma_tau(s*,0) < 0` where `s* = sqrt(a/(pi t2))`, `a = ln(10/7)` (so `e^{-pi t2 s*^2} = x := 7/10`).

*Proof.* Put `E0 = sum_m e^{-pi t2 (m+s*)^2}` and `u = 1 - s*`. The `m=0` term is `x`; for `m != 0`, `|m+s*| >= u + (|m|-1)` (for `m<=-1`) resp. `>= 1 >= u + ...`, so
`0 <= E0 - x <= delta := 2 e^{-pi t2 u^2}/(1-e^{-2 pi t2 u})`  (use `(u+j)^2 >= u^2 + 2uj`).
With `T = sum_{n != 0} e^{-pi t2 n^2} <= T_n := 2 e^{-pi t2}/(1-e^{-3 pi t2})`:

* `K/c = sum_n e^{-pi t2 n^2} sum_m e^{-pi t2 (m+s*)^2} e^{2 pi i t1 (m+s*) n}`; the `n=0` row is `E0`, the other rows have modulus `<= e^{-pi t2 n^2} E0`; `K` is real, so `K/c in E0 [1-T_n, 1+T_n]`.
* `lambda/c = K(0,0)/c` (same argument at `s=0`, `n=0`-row `1+T`): `lambda/c in [(1+T)(1-T), (1+T)^2] subset [1-T_n^2, (1+T_n)^2]`.
* Hence `|A/c - (2x-1)| <= eta_A := max( 2 delta + 2 (x+delta) T_n + T_n^2, 2 (x+1) T_n + T_n^2 )`.
* `B(s*,0) = sum G(m+s*, n-s*)`: the `(0,0)` term has modulus `c x^2`; the sum of all term moduli is `c E0(s*) E0(-s*) = c E0^2`.
  So `|B|/c >= x^2 - (E0^2 - x^2) >= b := x^2 - 2 x delta - delta^2`.

Thus `Sigma/c^2 <= (2x - 1 + eta_A)^2 - b^2` for all `t1`. The right-hand side is monotone non-increasing in `t2` (each of `delta, T_n, eta_A` decreases, `b` increases; `u` increases with `t2`, and `pi t2 u^2 = (sqrt(pi t2) - sqrt(a))^2` increases).
At `t2 = 3`, Arb (`cusp_arb.py`, 200 bits) certifies `eta_A < 0.009069`, `b > 0.48379` and a strictly negative bound. This is one evaluation of majorants uniform over all real `t1`, without a horizontal cell subdivision. QED.

Lean: the arithmetic step (`|a-2/5| <= 1/100`, `b >= 12/25` imply `a^2 - b^2 < 0`) is `lean/GSTMap/Cusp.lean`; the numerics `delta, T_n` are the only
inputs (hypotheses, Arb-checked). Cusp bound used in the tree: every cell with `t2 >= 3` is unstable; the tree is built on `t2 <= 51/16`.

## 5. What `Sigma_tau(q) < 0` means (instability mechanism)

On the Bloch fibre at momentum `q` the linearisation of `i d_t u = Pi(|u|^2 u)` has eigenfrequencies `omega = +-sqrt(Sigma_tau(q))`
(`omega^2 = Sigma`, the paper's convention). `Sigma_tau(q) < 0` therefore gives `omega = +-i sqrt(-Sigma_tau(q))`: one exponentially growing and
one exponentially decaying mode, i.e. spectral instability at Bloch momentum `q`. `Sigma >= 0` on the whole torus gives only real `omega` (spectral
stability). The mechanism for the shapes in the unstable region is the loss of positive definiteness of the quartic long-wave form `P4_tau`
(Section 2), i.e. the energy form `W` of Section 0 becoming indefinite on long-wave Bloch modes. Of the 25088 final U leaves, 22107 use witnesses near `|q|=1/512` and 2981 near `0.02`. These certify a mathematical spectral instability, not an observable growth rate. No theorem excludes a separate finite-frequency instability inside the remaining strip.

## 6. Statement of the partial stability diagram

* **Unstable region (CERTIFIED, computer-assisted).** `{Im tau >= 3}` (every Re tau; Lemma C), plus the union of the unstable leaves of the tree.
* **Stable region (CERTIFIED).** The union of the stable leaves: the closed box `Sigma>0` on torus minus 0 for every tau in the leaf.
* **Marginal set (CERTIFIED brackets).** On each vertical segment between an S box below and a U box above, closedness of the stable set gives a last stable point; this need not form a unique graph.
* **Quartic curve `Gamma4` (section 10).** The implicit equations `P4(1,r)=0`, `d_r P4(1,r)=0` describe an analytic branch. Identifying it with the full spectral boundary still requires control of every nonzero Bloch momentum.
* **Undecided set.** The explicit X leaves of the final tree. Its coordinate area is a certificate-dependent upper bound, not a physical transition width.

## 7. Final exact tree (`tree_smallq.jsonl`)

The domain `[0,1/2] x [13/16,51/16]` has 76 root squares and 40570 leaves: K=66, S=4178, U=25088, X=11238. Exact rational validation (`gen_data.py`, pure `Fraction`) gives zero partition errors and total area `19/16`. The respective coordinate areas are K=`16235067/268435456`, S=`1478129/67108864`, U=`296608283/268435456`, X=`5619/134217728`.

K is strictly outside F+. S means `Sigma_tau(q)>0` for every tau in the closed leaf and every nonzero q on the Bloch torus. U means an explicit rational q has negative Sigma throughout the closed leaf. X is undecided. S is interval-certified, single implementation, sampling-corroborated. U has independent Arb and mpmath.iv certificates.

Sampled vertical brackets are x=0.3: `[0.95751953125,0.95782470703125]`, x=0.4: `[1.026611328125,1.02685546875]`, x=0.5: `[1.0489501953125,1.04913330078125]`. These are brackets, not a uniform strip-width or uniqueness theorem.

The certified S and U unions resemble a petal and its exterior, but the X leaves are not assigned either sign. The existence of brackets does not prove uniqueness, smoothness, or equality of the true spectral boundary with `Gamma4`.

## 8. Independence and limits

The independent floating-point S-dual implementation samples all 4178 S leaves. Sampling is NUMERICAL and does not replace a continuous cover. A second-arithmetic attempt covered only 12 leaves; the remaining 4166 have no second-arithmetic coverage. S retains a single interval implementation. U has a separate cancellation-preserving interval verifier, described below. The analytic formulas and tail proofs remain part of the mathematical trust base.

## 9. Analytic cancellation, P4/P6, and a second rigorous U verifier

Write `tau=x+iy`, `q=t(s,xi)`, `w_mn=G_tau(m,n)`, and

    LK=-2 pi y m s + 2 pi i n(x s+xi),       QK=-pi y s^2,
    LB=(n-m)[2 pi y s+2 pi i(x s+xi)],       QB=-2 pi(y+i x)s^2.

Then `K(tq)=sum w exp(LK t+QK t^2)` and `B(tq)=sum w exp(LB t+QB t^2)`.
On compact subsets of `y>0`, every differentiated summand is bounded by a polynomial times
`exp(-a(m^2+n^2)+b(|m|+|n|))`, `a>0`. All differentiations, reindexings and Taylor expansions below
are therefore justified by locally uniform absolute convergence. The coefficients of `exp(Lt+Qt^2)`
satisfy `c0=1,c1=L,cj=(L c(j-1)+2Q c(j-2))/j`.

**Exact quadratic cancellation (PROVED).** Let `M2=sum w m^2=sum w n^2` and `M11=sum w mn`.
Pairing `n` with `-n` shows that `lambda` and `M2` are real and `M11` is purely imaginary.
Put `u=2 pi y s`, `v=2 pi(x s+xi)`. Direct expansion gives

    A2-B2=(u^2-v^2)M11-2 i u v M2+2 pi i x s^2 lambda.

The right side is purely imaginary; since `A2` is real, `A2=Re B2`.
The substitution `(m,n)->(-m,-n)` makes K and B even in t. Also `A0=B0=lambda`.
Consequently the constant and quadratic coefficients of Sigma vanish exactly. Multiplying the
convergent series gives the homogeneous identities

    P4=2 lambda(A4-Re B4)-(Im B2)^2,
    P6=2 lambda(A6-Re B6)+2 A2(A4-Re B4)-2 Im B2 Im B4.

These are identities for the infinite lattice sums, not fitted finite differences. The remaining
terms are `O(|q|^8)`, uniformly in angle and in compact tau sets. Lean `Quartic.lean` checks the
polynomial algebra once `A2=Re B2` is supplied; it does not formalize differentiation of the lattice sum.
If `P4(q)<0`, taking `t>0` small gives `Sigma(tq)<0`. Hence spectral stability necessarily implies
`P4>=0` in every direction. Positive-definite P4 implies positivity only for sufficiently small q.

**Independent interval semantics.** `cancel_iv.py` has no Arb or original-engine import. It expands
the scalar t series FIRST and inserts the exact cancellations above before enclosing tau; this is why
its tau gradient retains the fourth power of the small witness momentum. Only the quartic uses a
centred tau mean-value enclosure. Higher even coefficients through degree 2N are enclosed directly.
It reuses the supplied rational witness, not the original search or Taylor implementation.

Here is a proof of its tail bounds. For `|t|<=1` complex, either K or B summand is bounded in modulus by

    sqrt(y)/pi exp(q) exp[-a(m^2+n^2)+beta(|m|+|n|)],
    a=pi*y_min, beta=2*pi*(y_max*|s|+|x*s+xi|)_max,
    q=2*pi*(y_max+|x|_max)*s^2.

For the one-dimensional sequence, ratios after the first term and after M are
`r0=exp(-3a+beta)` and `rM=exp(-(2M+3)a+beta)`. If both are below one, the full 1D sum is at most
`H=1+2exp(-a+beta)/(1-r0)` and the omitted 1D tail is at most
`E=2exp(-a(M+1)^2+beta(M+1))/(1-rM)`. The omitted 2D lattice sum is bounded by `2EH`.
Multiplying by `sqrt(y)_max/pi*exp(q)` gives T. Cauchy's coefficient bound on the unit t circle
bounds each omitted Taylor coefficient by T. The omitted-lattice Taylor remainder is at most `(N+2)T`.

For a finite term, the exponential-series remainder is bounded by `rem(l,q,N)`:
sum `q^k/k! * exp(l)*l^(N-2k+1)/(N-2k+1)!` for `0<=k<=floor(N/2)`, plus
`exp(l+q)*q^(floor(N/2)+1)/(floor(N/2)+1)!`; here `l>=|L|`, `q>=|Q|`.
The inequality `sum_{j>=h} l^j/j! <= exp(l) l^h/h!` proves it coefficientwise.
For tau derivatives, bound the logarithmic derivative by
`2*pi*(m^2+n^2)+4*pi*s^2+1/(2*y_min)`, consume half the Gaussian decay, and use
`sup_{z>=0} z exp(-pi*y_min*z/2)=2/(pi*y_min*e)`.
Thus the half-decay tail times `4/(y_min*e)+4*pi*s^2+1/(2*y_min)` bounds each differentiated coefficient.
All constants are evaluated with `mpmath.iv`, 35 decimal digits, not ordinary floating point.

If `|A-A_N|<=epsA`, `|B-B_N|<=epsB`, then the error in Sigma is at most
`2|A_N|epsA+epsA^2+2|B_N|epsB+epsB^2`. Its strictly negative upper bound proves each U claim.
Tau subdivision is exact rational subdivision; every subbox must pass. All 25088 U leaves have strictly negative upper bounds. The negative-control script accepts a true U box and rejects a stable box, a curve-straddling box, a zero witness and an inadequate truncation.

## 10. Parametric enclosure of the quartic branch (not yet the full spectral boundary)

Set `f(x,y,r)=P4_(x+iy)(1,r)`. `gamma_implicit.py` encloses the two equations `F=(f,f_r)=0`
on the exact x interval `[37/128,1/2]`. Each patch has rational x centre and half-width hx,
rational seeds `(y0,r0)`, and radii `hy=4hx`, `hr=8hx`. Coefficients and both tau derivatives
are enclosed by the mpmath.iv expansion in section 9, shared with `cancel_iv.py`; this is a single interval backend, not an independent second Gamma4 certificate. Five exact interpolation
nodes `r=-1,-1/2,0,1/2,1` and a rational inverse Vandermonde matrix recover all coefficients.
The small evaluation scale `s=1/64` keeps Cauchy lattice-tail bounds sharp; division by `s^4` is exact.
Centred coefficient enclosures retain the shared tau dependence.

Let J enclose `d_(y,r)F` throughout a patch, R be any rational approximate inverse, z0 the rational
seed, and Z the centred `(y,r)` box. For every x in the patch the code verifies

    K=-R F(x,z0)+(I-RJ)(Z-z0) is strictly inside Z-z0,
    max_i sum_j |(I-RJ)_ij| radius_j/radius_i < 1.

The associated Newton map is a contraction from Z into itself, so it has exactly one zero for each
x. These are ordinary interval inclusions, independent of the float seed's accuracy. The code also
checks `f_y<0`, `f_rr>0`. At a root the Jacobian determinant is `f_y*f_rr`, so the implicit-function
theorem gives a real-analytic branch. `verify_gamma_chain.py` proves existence of a root in the
intersection of each neighbouring pair of closed (y,r) boxes at their exact common x endpoint.
Uniqueness in each patch then glues them; merely overlapping seed boxes would NOT suffice.
The chain has 396 joins, two unit-arc bracket endpoint roots and three additional fixed-x roots. The negative-control script accepts a genuine double root and adjacent join, and rejects a shifted centre, wrong null direction, omitted patch and shifted common endpoint.

All 397 patches certify `f_y < -1.6427`, `f_rr > 12.7189` and `P6(1,r) > 0.74577`, including patches below x=0.3. The sextic bound uses the unnormalised direction `(1,r)`. These decimals conservatively summarize interval signs; the verifier uses the interval values.

To identify the branch as a global quartic minimum, divide at its double root r*:

    f(r)=(r-r*)^2 [c4*r^2+b*r+d],
    b=c3+2r* c4, d=c2+2r* c3+3r*^2 c4,
    D=4 c4 d-b^2.

Every patch proves `c4>0`. Patches lying wholly at `x>=0.3` prove `D>0`.
Patches meeting `x<0.3` instead prove along the branch `D'>0` and `(x^2+y^2)'>0`, using
`y'=-f_x/f_y`, `r'=-(f_rx+f_ry*y')/f_rr` and interval differentiation of the displayed D.
At `|tau|=1`, combining modular S with reflection R gives `P4(s,xi)=P4(xi,s)`.
The enclosed r* is nonzero and not +/-1; its reciprocal is therefore a second distinct double root.
It follows exactly that `D=0` there. The endpoint certificates and gluing give a unique arc crossing
in `x in (297/1000,2971/10000)`. There cannot be another crossing for `x>=0.3`, since `D>0`
contradicts the exact arc identity `D=0`. Thus on the F+ portion of this branch P4 is nonnegative,
with one null projective direction in the interior and two at the unit-arc endpoint.
Inside the arc the same analytic solution branch need not be a global-minimum-zero branch.

Outward-rounded rational decimal enclosures from the separate fixed-x Krawczyk checks are

| x | lower y | upper y |
|---|---|---|
| 0.3 | 0.95772312050217 | 0.95772312050218 |
| 0.4 | 1.02681009799052 | 1.02681009799054 |
| 0.5 | 1.04912619784646 | 1.04912619784647 |

These locate the exact quartic roots, not necessarily the full spectral boundary. At x=0.297
the root has `|tau|^2-1 < -3.6165e-5`, while at x=0.2971 it has `|tau|^2-1 > 2.0713e-4`.

P6 is positive in the null direction throughout the certified branch; modular symmetry supplies the
second null direction at the endpoint. By continuity, P6 stays positive near each null angle;
away from them P4 has a positive lower bound. The uniform `O(rho^8)` remainder therefore proves
`Sigma>0` for every sufficiently small nonzero q at each branch point. Immediately above the branch,
`f_y<0` gives a negative quartic direction and spectral instability. On the branch the null-direction
frequency scales as `|q|^3`, not the generic `|q|^2`.

**Remaining logical gap:** none of this rules out `Sigma<0` at a finite nonzero momentum on the
branch or between the S union and the branch. Therefore `Gamma4` is a certified long-wave transition
curve, NOT yet a proved full spectral-stability boundary.

