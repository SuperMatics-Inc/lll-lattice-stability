import GSTStability.Coverage

noncomputable section

namespace GST

/-- Finite Taylor lower-bound rule used for every accepted far cell.
The coefficients and analytic remainder enclosure are hypotheses, not axioms. -/
theorem far_cell_lower_bound (sigma c₀ error : ℝ) (terms bounds : ℕ → ℝ) (N : ℕ)
    (hterms : ∀ i < N, |terms i| ≤ bounds i)
    (herror : |sigma - (c₀ + ∑ i ∈ Finset.range N, terms i)| ≤ error) :
    c₀ - (∑ i ∈ Finset.range N, bounds i) - error ≤ sigma := by
  have hs : -(∑ i ∈ Finset.range N, bounds i) ≤ ∑ i ∈ Finset.range N, terms i := by
    rw [← Finset.sum_neg_distrib]
    apply Finset.sum_le_sum
    intro i hi
    have h := hterms i (Finset.mem_range.mp hi)
    have h' := neg_abs_le (terms i)
    linarith
  have he := (abs_le.mp herror).1
  linarith

/-- Every decimal used here is an exact rational, proved by kernel arithmetic. -/
theorem exact_rational_thresholds :
    (0 : ℝ) < 33019570 / 100000000 ∧
    (0 : ℝ) < 384161 / 100000000 ∧
    (0 : ℝ) < 104595 / 1000000000000000 ∧
    (0 : ℝ) < 424467 / 100000000000000 ∧
    (0 : ℝ) < 35094829 / 100000000 ∧
    (0 : ℝ) < 13501414 / 100000000 ∧
    (0 : ℝ) < 1 / 100 ∧ (0 : ℝ) < 1 / 40 := by norm_num

/-- Exact downward rounding from the *printed enclosing balls* in the private
historical provenance logs `cert_sigma_{m04_095,p04_100}.log`, lines 12,16,24/25.
Regenerate these enclosures with `certificates/cert_sigma.py` using the two
parameter lists in `scripts/reproduce.py`; see the root README reproduction commands
and `lean/README.md` for the conditional formal scope.
This checks arithmetic on the reported enclosures, not their theta provenance. -/
theorem printed_ball_threshold_checks :
    (35094830 : ℝ) / 10^8 - 363 / 10^11 > 35094829 / 10^8 ∧
    (13501415 : ℝ) / 10^8 - 394 / 10^11 > 13501414 / 10^8 ∧
    (33019571 : ℝ) / 10^8 - 414 / 10^11 > 33019570 / 10^8 ∧
    (38416147 : ℝ) / 10^10 - 404 / 10^13 > 384161 / 10^8 ∧
    (104596 : ℝ) / 10^15 - 151 / 10^18 > 104595 / 10^15 ∧
    (424468 : ℝ) / 10^14 - 4 / 10^15 > 424467 / 10^14 := by norm_num

/-- Assembly on the half-torus from a uniform near certificate and eight
terminating adaptive trees. No point in the punctured near square is lost. -/
theorem half_torus_positive (sigma : Point → ℝ) (r : ℝ)
    (near : ∀ p, 0 < radius p → radius p ≤ r → 0 < sigma p)
    (trees : ∀ i < 2, ∀ j < 4, CertifiedTree r sigma (root i j))
    (p : Point) (hp : 0 < radius p)
    (hx0 : 0 ≤ p.1) (hx1 : p.1 ≤ 1 / 2)
    (hy0 : -(1 / 2) ≤ p.2) (hy1 : p.2 ≤ 1 / 2) : 0 < sigma p := by
  rcases half_torus_covered r sigma trees p hx0 hx1 hy0 hy1 with h | h
  · exact near p hp h
  · exact h

/-- Evenness supplies the other half of the closed fundamental square. -/
theorem fundamental_square_positive (sigma : Point → ℝ) (r : ℝ)
    (near : ∀ p, 0 < radius p → radius p ≤ r → 0 < sigma p)
    (trees : ∀ i < 2, ∀ j < 4, CertifiedTree r sigma (root i j))
    (even : ∀ p : Point, sigma (-p.1, -p.2) = sigma p)
    (p : Point) (hp : 0 < radius p)
    (hx : |p.1| ≤ 1 / 2) (hy : |p.2| ≤ 1 / 2) : 0 < sigma p := by
  by_cases hx0 : 0 ≤ p.1
  · exact half_torus_positive sigma r near trees p hp hx0
      (abs_le.mp hx).2 (abs_le.mp hy).1 (abs_le.mp hy).2
  · have hrad : radius (-p.1, -p.2) = radius p := by simp [radius]
    have hpos := half_torus_positive sigma r near trees (-p.1, -p.2)
      (by simpa only [hrad] using hp) (by linarith)
      (by have := (abs_le.mp hx).1; linarith)
      (by have := (abs_le.mp hy).2; linarith)
      (by have := (abs_le.mp hy).1; linarith)
    simpa only [even] using hpos

/-- Centering by an exact integer; no rounding of torus representatives. -/
def centered (x : ℝ) : ℝ := x - (⌊x + 1 / 2⌋ : ℤ)

theorem centered_bounds (x : ℝ) : |centered x| ≤ 1 / 2 := by
  have h₁ := Int.floor_le (x + 1 / 2)
  have h₂ := Int.lt_floor_add_one (x + 1 / 2)
  unfold centered
  apply abs_le.mpr
  constructor <;> linarith

theorem radius_positive_of_ne_zero (p : Point) (hp : p ≠ (0, 0)) : 0 < radius p := by
  have h0 : 0 ≤ radius p := le_trans (abs_nonneg p.1) (le_max_left _ _)
  by_contra h
  have hz : radius p = 0 := le_antisymm (le_of_not_gt h) h0
  have hax : |p.1| ≤ 0 := calc
    |p.1| ≤ radius p := le_max_left _ _
    _ = 0 := hz
  have hay : |p.2| ≤ 0 := calc
    |p.2| ≤ radius p := le_max_right _ _
    _ = 0 := hz
  have hx : p.1 = 0 := abs_eq_zero.mp (le_antisymm hax (abs_nonneg p.1))
  have hy : p.2 = 0 := abs_eq_zero.mp (le_antisymm hay (abs_nonneg p.2))
  exact hp (Prod.ext hx hy)

/-- Theorem 1.1, determinant-positivity part, conditional on analytic enclosures.
A periodic function on R² is positive on T²\{0} precisely when it is positive
at every p not in Z². Periodicity, evenness, the near bounds and tree-leaf
enclosures are all explicit arguments. The theta function is NOT postulated. -/
theorem punctured_torus_positive (sigma : Point → ℝ) (r : ℝ)
    (near : ∀ p, 0 < radius p → radius p ≤ r → 0 < sigma p)
    (trees : ∀ i < 2, ∀ j < 4, CertifiedTree r sigma (root i j))
    (even : ∀ p : Point, sigma (-p.1, -p.2) = sigma p)
    (periodic : ∀ (p : Point) (m n : ℤ), sigma (p.1 + m, p.2 + n) = sigma p)
    (p : Point) (hp : ∀ m n : ℤ, p ≠ ((m : ℝ), (n : ℝ))) : 0 < sigma p := by
  let m : ℤ := ⌊p.1 + 1 / 2⌋
  let n : ℤ := ⌊p.2 + 1 / 2⌋
  let q : Point := (centered p.1, centered p.2)
  have hq : q ≠ (0, 0) := by
    intro hz
    have hx := congrArg Prod.fst hz
    have hy := congrArg Prod.snd hz
    apply hp m n
    apply Prod.ext
    · change p.1 = (m : ℝ)
      change p.1 - (m : ℝ) = 0 at hx
      linarith
    · change p.2 = (n : ℝ)
      change p.2 - (n : ℝ) = 0 at hy
      linarith
  have hpos := fundamental_square_positive sigma r near trees even q
    (radius_positive_of_ne_zero q hq) (centered_bounds p.1) (centered_bounds p.2)
  have heq := periodic q m n
  have hx : q.1 + (m : ℝ) = p.1 := by dsimp [q, centered, m]; ring
  have hy : q.2 + (n : ℝ) = p.2 := by dsimp [q, centered, n]; ring
  rw [hx, hy] at heq
  rw [Prod.mk.eta] at heq
  rw [heq]
  exact hpos

end GST
