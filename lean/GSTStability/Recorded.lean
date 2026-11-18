import GSTStability.Trace
import GSTStability.Main

/-! Generated exact subdivision shapes from `data/trace-minus.json` and
`data/trace-plus.json`. Only near/far labels and dyadic coordinates are used.
No floating-point theta bound is imported into the kernel. -/
namespace GST
set_option maxRecDepth 4096
set_option maxHeartbeats 0

def minusTrace (i j : ℕ) : QuadTree :=
  match i, j with
  | 0, 0 => (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split .far .far .far .far) (.split .far .far .far .far))
  | 0, 1 => (.split (.split (.split .far .far .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) .far (.split .far .far .far .far)) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far .near .far .near) (.split .far .far .far .far) (.split .far .near .far .near)) .near (.split .far (.split .far .near .far .far) .far .far) (.split (.split .near .near .far .far) (.split .near .near .far .far) (.split .far .far .far .far) (.split .far .far .far .far))) (.split .far (.split .far .far .far .far) .far .far) (.split (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far .far) (.split .far .far .far .far) (.split .far .far .far .far))) (.split .far (.split .far .far .far .far) .far .far) (.split (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far .far) (.split .far .far .far .far) (.split .far .far .far .far))) (.split .far (.split .far .far .far .far) .far .far) (.split (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far .far) (.split .far .far .far .far) (.split .far .far .far .far))) (.split (.split .far .far .far .far) (.split .far .far .far .far) .far (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far (.split .far .far .far .far)) (.split .far .far .far .far) (.split .far .far .far .far)))
  | 0, 2 => (.split (.split (.split (.split (.split .near (.split (.split .near .far .near .far) (.split .far .far .far .far) (.split .near .far .near .far) (.split .far .far .far .far)) (.split (.split .near .near (.split .far .far .far .far) .far) (.split .near .near .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .near .far .far .far) .far .far .far)) (.split (.split (.split .far .far .far .far) .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) .far (.split .far .far .far .far) .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) .far .far .far)) (.split (.split (.split .far .far .far .far) .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) .far)) (.split (.split (.split .far .far .far .far) .far (.split .far .far .far .far) .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far))) (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split .far .far .far .far) (.split (.split .far .far .far .far) .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)))
  | 0, 3 => (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) .far .far .far))
  | 1, 0 => (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far (.split .far .far .far .far) .far) (.split .far .far .far .far))
  | 1, 1 => (.split (.split .far .far .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split .far .far .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)))
  | 1, 2 => (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) .far .far .far))
  | 1, 3 => (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split .far .far (.split .far .far .far .far) .far) (.split .far (.split .far .far .far .far) .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) .far))
  | _, _ => .far

theorem minus_trace_geometry (i j : ℕ) (hi : i < 2) (hj : j < 4) :
    (minusTrace i j).geometry (1 / 100) (rationalRoot i j) = true := by
  interval_cases i <;> interval_cases j <;> decide +kernel

theorem minus_recorded_leaf_counts : rootLeafCounts minusTrace = (20, 756) := by decide

theorem minus_recorded_trees (sigma : Point → ℝ)
    (far : ∀ i < 2, ∀ j < 4, (minusTrace i j).farEnclosures sigma (rationalRoot i j)) :
    ∀ i < 2, ∀ j < 4, CertifiedTree (1 / 100) sigma (root i j) := by
  intro i hi j hj
  have h := checked_trace_sound (minusTrace i j) (1 / 100) sigma
    (rationalRoot i j) (minus_trace_geometry i j hi hj) (far i hi j hj)
  simpa only [rational_root_real, Rat.cast_div, Rat.cast_one, Rat.cast_ofNat] using h

/-- The actual minus certificate subdivision exhausts the punctured torus,
conditional only on near and per-far-leaf analytic bounds and symbol symmetries. -/
theorem minus_punctured_torus_positive (sigma : Point → ℝ)
    (near : ∀ p, 0 < radius p → radius p ≤ 1 / 100 →
      (33019570 / 100000000 : ℝ) * radius p ^ 4 ≤ sigma p)
    (far : ∀ i < 2, ∀ j < 4, (minusTrace i j).farLowerBounds
      (104595 / 1000000000000000) sigma (rationalRoot i j))
    (even : ∀ p : Point, sigma (-p.1, -p.2) = sigma p)
    (periodic : ∀ (p : Point) (m n : ℤ), sigma (p.1 + m, p.2 + n) = sigma p)
    (p : Point) (hp : ∀ m n : ℤ, p ≠ ((m : ℝ), (n : ℝ))) : 0 < sigma p := by
  have hn : ∀ q, 0 < radius q → radius q ≤ 1 / 100 → 0 < sigma q := by
    intro q hq hqr
    exact lt_of_lt_of_le (mul_pos (by norm_num) (pow_pos hq 4)) (near q hq hqr)
  have hf : ∀ i < 2, ∀ j < 4, (minusTrace i j).farEnclosures sigma (rationalRoot i j) := by
    intro i hi j hj
    exact far_lower_bounds_positive _ _ (by norm_num) _ _ (far i hi j hj)
  exact punctured_torus_positive sigma (1 / 100) hn
    (minus_recorded_trees sigma hf) even periodic p hp

def plusTrace (i j : ℕ) : QuadTree :=
  match i, j with
  | 0, 0 => (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split .far (.split .far .far .far .far) .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)))
  | 0, 1 => (.split (.split (.split .far .far .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far .far)) (.split (.split (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far (.split .far .far .far .far)) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) .near (.split .far .far .far .far) .near) .near (.split .far .near .far (.split .far .far .far .far)) (.split .near .near (.split .far .far .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) .far))) (.split (.split .far .far .far .far) (.split .far .far .far .far) .far (.split .far .far .far .far)) (.split (.split .far (.split .far .far .far .far) .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split .far .far .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) .far))) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far (.split .far .far .far .far) .far (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split .far .far .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) .far .far))) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far (.split .far .far .far .far) .far (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split .far .far .far .far) (.split .far (.split .far .far .far .far) .far .far)))
  | 0, 2 => (.split (.split (.split (.split .near (.split .near (.split .far .far .far .far) .near (.split .far .far .far .far)) (.split .near .near (.split .far .far .far .far) (.split .far .far .far .far)) (.split .near .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) .far .far .far)) (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split (.split .far .far .far .far) .far (.split .far .far .far .far) .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) .far))
  | 0, 3 => (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) .far .far .far) (.split .far .far .far .far))
  | 1, 0 => (.split (.split .far .far .far (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) .far (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) .far .far))
  | 1, 1 => (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) .far (.split .far .far .far .far)) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)))
  | 1, 2 => (.split (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split (.split .far .far .far .far) .far .far .far) (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far)) (.split .far .far .far .far))
  | 1, 3 => (.split (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far .far) (.split .far .far .far (.split .far .far .far .far)))
  | _, _ => .far

theorem plus_trace_geometry (i j : ℕ) (hi : i < 2) (hj : j < 4) :
    (plusTrace i j).geometry (1 / 40) (rationalRoot i j) = true := by
  interval_cases i <;> interval_cases j <;> decide +kernel

theorem plus_recorded_leaf_counts : rootLeafCounts plusTrace = (12, 707) := by decide

theorem plus_recorded_trees (sigma : Point → ℝ)
    (far : ∀ i < 2, ∀ j < 4, (plusTrace i j).farEnclosures sigma (rationalRoot i j)) :
    ∀ i < 2, ∀ j < 4, CertifiedTree (1 / 40) sigma (root i j) := by
  intro i hi j hj
  have h := checked_trace_sound (plusTrace i j) (1 / 40) sigma
    (rationalRoot i j) (plus_trace_geometry i j hi hj) (far i hi j hj)
  simpa only [rational_root_real, Rat.cast_div, Rat.cast_one, Rat.cast_ofNat] using h

/-- The actual plus certificate subdivision exhausts the punctured torus,
conditional only on near and per-far-leaf analytic bounds and symbol symmetries. -/
theorem plus_punctured_torus_positive (sigma : Point → ℝ)
    (near : ∀ p, 0 < radius p → radius p ≤ 1 / 40 →
      (384161 / 100000000 : ℝ) * radius p ^ 4 ≤ sigma p)
    (far : ∀ i < 2, ∀ j < 4, (plusTrace i j).farLowerBounds
      (424467 / 100000000000000) sigma (rationalRoot i j))
    (even : ∀ p : Point, sigma (-p.1, -p.2) = sigma p)
    (periodic : ∀ (p : Point) (m n : ℤ), sigma (p.1 + m, p.2 + n) = sigma p)
    (p : Point) (hp : ∀ m n : ℤ, p ≠ ((m : ℝ), (n : ℝ))) : 0 < sigma p := by
  have hn : ∀ q, 0 < radius q → radius q ≤ 1 / 40 → 0 < sigma q := by
    intro q hq hqr
    exact lt_of_lt_of_le (mul_pos (by norm_num) (pow_pos hq 4)) (near q hq hqr)
  have hf : ∀ i < 2, ∀ j < 4, (plusTrace i j).farEnclosures sigma (rationalRoot i j) := by
    intro i hi j hj
    exact far_lower_bounds_positive _ _ (by norm_num) _ _ (far i hi j hj)
  exact punctured_torus_positive sigma (1 / 40) hn
    (plus_recorded_trees sigma hf) even periodic p hp

end GST
