import Mathlib.Data.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-! Algebraic core of PROOFS.md section 9. Analytic Taylor/moment identities
and coefficient bounds remain explicit external inputs, not axioms. -/
namespace GSTMap

theorem quartic_coefficient (l a2 a4 b2r b2i b4r : ℝ) (h2 : a2 = b2r) :
    (2*l*a4+a2^2)-(2*l*b4r+(b2r^2+b2i^2)) =
      2*l*(a4-b4r)-b2i^2 := by
  rw [h2]
  ring

theorem sextic_coefficient (l a2 a4 a6 b2r b2i b4r b4i b6r : ℝ)
    (h2 : a2 = b2r) :
    (2*l*a6+2*a2*a4)-(2*l*b6r+2*(b2r*b4r+b2i*b4i)) =
      2*l*(a6-b6r)+2*a2*(a4-b4r)-2*b2i*b4i := by
  rw [h2]
  ring

theorem quartic_remainder_sign {p e r sigma : ℝ} (hp : p < 0)
    (he : e ≤ -p/2) (hr : 0 < r) (hrem : sigma ≤ r^4*(p+e)) :
    sigma < 0 := by
  have h : p+e < 0 := by linarith
  exact lt_of_le_of_lt hrem (mul_neg_of_pos_of_neg (pow_pos hr 4) h)

end GSTMap
