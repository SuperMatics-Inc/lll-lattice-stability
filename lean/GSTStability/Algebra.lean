import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Tactic

noncomputable section

/-! Algebraic part of Proposition 3.1. These are homogeneous jets along an
arbitrary direction, not assertions that theta series have these jets. The
analytic moment identities supply the common real quadratic coefficient. -/
namespace GST

/-- K₂ obtained from the symmetric moments M20=M02=mu and M11=i nu. -/
def momentK₂ (pi t₁ t₂ lam mu nu s xi : ℝ) : ℝ :=
  (-pi * t₂ * lam + 2 * pi ^ 2 * (t₂ ^ 2 - t₁ ^ 2) * mu +
    4 * pi ^ 2 * t₁ * t₂ * nu) * s ^ 2 -
  4 * pi ^ 2 * (t₁ * mu - t₂ * nu) * s * xi - 2 * pi ^ 2 * mu * xi ^ 2

/-- Real part of -2 pi i conjugate(tau) lam s² + Psi₂ z², with
Psi₂=-4 pi²(mu-i nu) and z=xi+(t₁-i t₂)s. Complex multiplication is
written explicitly in real coordinates. -/
def momentReB₂ (pi t₁ t₂ lam mu nu s xi : ℝ) : ℝ :=
  -2 * pi * t₂ * lam * s ^ 2 -
  4 * pi ^ 2 * mu * ((xi + t₁ * s) ^ 2 - (t₂ * s) ^ 2) +
  8 * pi ^ 2 * nu * t₂ * s * (xi + t₁ * s)

/-- The algebraic moment reduction proves A₂=Re B₂ in every direction.
The symmetries and termwise differentiation identifying the moments are
analytic prerequisites, not part of this ring identity. -/
theorem moment_quadratic_cancellation (pi t₁ t₂ lam mu nu s xi : ℝ) :
    2 * momentK₂ pi t₁ t₂ lam mu nu s xi =
      momentReB₂ pi t₁ t₂ lam mu nu s xi := by
  unfold momentK₂ momentReB₂
  ring

/-- Even jets, retaining precisely the coefficients relevant through degree four. -/
def evenJet (a₀ a₂ a₄ : ℝ) (n : ℕ) : ℝ :=
  if n = 0 then a₀ else if n = 2 then a₂ else if n = 4 then a₄ else 0

def convolution (a b : ℕ → ℝ) (n : ℕ) : ℝ :=
  ∑ i ∈ Finset.range (n + 1), a i * b (n - i)

/-- Coefficients of A² - (Re B)² - (Im B)², assuming A₂ = Re B₂. -/
def sigmaJet (lam a₂ a₄ b₄ v₂ v₄ : ℝ) (n : ℕ) : ℝ :=
  convolution (evenJet lam a₂ a₄) (evenJet lam a₂ a₄) n -
  convolution (evenJet lam a₂ b₄) (evenJet lam a₂ b₄) n -
  convolution (evenJet 0 v₂ v₄) (evenJet 0 v₂ v₄) n

/-- Exact cancellation, not the assertion that interval balls contain zero. -/
theorem cancellation_degree_le_three (lam a₂ a₄ b₄ v₂ v₄ : ℝ)
    (n : ℕ) (hn : n ≤ 3) : sigmaJet lam a₂ a₄ b₄ v₂ v₄ n = 0 := by
  interval_cases n <;>
    norm_num [sigmaJet, convolution, evenJet, Finset.sum_range_succ]

/-- The degree-four coefficient is the P4 identity of Proposition 3.1. -/
theorem quartic_identity (lam a₂ a₄ b₄ v₂ v₄ : ℝ) :
    sigmaJet lam a₂ a₄ b₄ v₂ v₄ 4 = 2 * lam * (a₄ - b₄) - v₂ ^ 2 := by
  norm_num [sigmaJet, convolution, evenJet, Finset.sum_range_succ]
  ring

/-- The elementary cancellation responsible for the negative imaginary square. -/
theorem quadratic_norm_cancellation (a₂ v₂ : ℝ) :
    a₂ ^ 2 - (a₂ ^ 2 + v₂ ^ 2) = -(v₂ ^ 2) := by ring

/-- Degree >20 requires at least eleven linear/quadratic exponent factors. -/
theorem omitted_factor_count (linear quadratic : ℕ)
    (h : 20 < linear + 2 * quadratic) : 11 ≤ linear + quadratic := by omega

end GST
