import GSTStability.Algebra
import Mathlib.Analysis.SpecialFunctions.Exp

noncomputable section

namespace GST

/-- A finite nonnegative-coefficient polynomial, used after dividing by rho^4. -/
def majorant (c : ℕ → ℝ) (N : ℕ) (rho : ℝ) : ℝ :=
  ∑ i ∈ Finset.range N, c i * rho ^ i

theorem majorant_nonneg (c : ℕ → ℝ) (N : ℕ) (rho : ℝ)
    (hc : ∀ i, 0 ≤ c i) (hrho : 0 ≤ rho) : 0 ≤ majorant c N rho := by
  exact Finset.sum_nonneg fun i _ => mul_nonneg (hc i) (pow_nonneg hrho i)

theorem majorant_mono (c : ℕ → ℝ) (N : ℕ) (rho r : ℝ)
    (hc : ∀ i, 0 ≤ c i) (hrho : 0 ≤ rho) (hle : rho ≤ r) :
    majorant c N rho ≤ majorant c N r := by
  apply Finset.sum_le_sum
  intro i _
  exact mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hrho hle i) (hc i)

/-- The product-error estimate E=2(Abar+Bbar)eps+2eps², factored by rho^4.
Here eps=rho^11*e; e is a nonnegative monotone envelope. -/
def normalizedError (A B e rho : ℝ) : ℝ :=
  2 * (A + B) * rho ^ 7 * e + 2 * rho ^ 18 * e ^ 2

theorem product_error_factorization (A B e rho : ℝ) :
    2 * (A + B) * (rho ^ 11 * e) + 2 * (rho ^ 11 * e) ^ 2 =
      rho ^ 4 * normalizedError A B e rho := by
  unfold normalizedError
  ring

theorem normalizedError_mono (A B e Ar Br er rho r : ℝ)
    (hA : 0 ≤ A) (hB : 0 ≤ B) (he : 0 ≤ e) (hrho : 0 ≤ rho)
    (hAr : A ≤ Ar) (hBr : B ≤ Br) (her : e ≤ er) (hr : rho ≤ r) :
    normalizedError A B e rho ≤ normalizedError Ar Br er r := by
  have hAr0 : 0 ≤ Ar := le_trans hA hAr
  have hBr0 : 0 ≤ Br := le_trans hB hBr
  have hr0 : 0 ≤ r := le_trans hrho hr
  unfold normalizedError
  gcongr

/-- Transfer an exact lower endpoint at r to a uniform rho^4 lower enclosure.
The proofs for the two recorded subdivisions use this bound. -/
theorem near_lower_bound (sigma p lower rho r : ℝ) (c : ℕ → ℝ) (N : ℕ)
    (A B e : ℝ → ℝ)
    (hrho : 0 ≤ rho) (hr : rho ≤ r)
    (hc : ∀ i, 0 ≤ c i)
    (hA : 0 ≤ A rho) (hB : 0 ≤ B rho) (he : 0 ≤ e rho)
    (hAm : A rho ≤ A r) (hBm : B rho ≤ B r) (hem : e rho ≤ e r)
    (hTaylor : rho ^ 4 * (p - majorant c N rho -
      normalizedError (A rho) (B rho) (e rho) rho) ≤ sigma)
    (hmargin : lower ≤ p - majorant c N r -
      normalizedError (A r) (B r) (e r) r) : lower * rho ^ 4 ≤ sigma := by
  have hp := majorant_mono c N rho r hc hrho hr
  have he' := normalizedError_mono (A rho) (B rho) (e rho)
    (A r) (B r) (e r) rho r hA hB he hrho hAm hBm hem hr
  have hb : lower ≤ p - majorant c N rho -
      normalizedError (A rho) (B rho) (e rho) rho := by linarith
  calc
    lower * rho ^ 4 = rho ^ 4 * lower := mul_comm _ _
    _ ≤ rho ^ 4 * (p - majorant c N rho -
        normalizedError (A rho) (B rho) (e rho) rho) :=
      mul_le_mul_of_nonneg_left hb (pow_nonneg hrho 4)
    _ ≤ sigma := hTaylor

/-- Uniform version of (4.2): no absolute error is divided by arbitrarily small rho.
The pointwise theta/Taylor enclosure is the explicit hypothesis hTaylor. -/
theorem near_positive (sigma p rho r : ℝ) (c : ℕ → ℝ) (N : ℕ)
    (A B e : ℝ → ℝ)
    (hrho : 0 < rho) (hr : rho ≤ r)
    (hc : ∀ i, 0 ≤ c i)
    (hA : 0 ≤ A rho) (hB : 0 ≤ B rho) (he : 0 ≤ e rho)
    (hAm : A rho ≤ A r) (hBm : B rho ≤ B r) (hem : e rho ≤ e r)
    (hTaylor : rho ^ 4 * (p - majorant c N rho -
      normalizedError (A rho) (B rho) (e rho) rho) ≤ sigma)
    (hmargin : 0 < p - majorant c N r -
      normalizedError (A r) (B r) (e r) r) : 0 < sigma := by
  have hb := near_lower_bound sigma p
    (p - majorant c N r - normalizedError (A r) (B r) (e r) r)
    rho r c N A B e hrho.le hr hc hA hB he hAm hBm hem hTaylor (le_refl _)
  exact lt_of_lt_of_le (mul_pos hmargin (pow_pos hrho 4)) hb

/-- Explicit exponential envelope for one omitted lattice summand after factoring
rho^11. Gaussian summability and shell-tail bounds remain analytic inputs. -/
def exponentialEnvelope (alpha beta rho : ℝ) : ℝ :=
  (alpha + beta * rho) ^ 11 * Real.exp (alpha * rho + beta * rho ^ 2) / 39916800

theorem exponentialEnvelope_nonneg (alpha beta rho : ℝ)
    (ha : 0 ≤ alpha) (hb : 0 ≤ beta) (hr : 0 ≤ rho) :
    0 ≤ exponentialEnvelope alpha beta rho := by
  unfold exponentialEnvelope
  positivity

theorem exponentialEnvelope_mono (alpha beta rho r : ℝ)
    (ha : 0 ≤ alpha) (hb : 0 ≤ beta) (hrho : 0 ≤ rho) (hr : rho ≤ r) :
    exponentialEnvelope alpha beta rho ≤ exponentialEnvelope alpha beta r := by
  have hr0 : 0 ≤ r := le_trans hrho hr
  unfold exponentialEnvelope
  gcongr

/-- The standard y^11 exp(y)/11! remainder majorant factors with rho^11,
not a constant error. This is the key reason division by rho^4 is harmless. -/
theorem exponential_remainder_factorization (alpha beta rho : ℝ) :
    (alpha * rho + beta * rho ^ 2) ^ 11 *
      Real.exp (alpha * rho + beta * rho ^ 2) / 39916800 =
    rho ^ 11 * exponentialEnvelope alpha beta rho := by
  have h : alpha * rho + beta * rho ^ 2 = rho * (alpha + beta * rho) := by ring
  rw [h, mul_pow]
  unfold exponentialEnvelope
  rw [← h]
  ring

end GST
