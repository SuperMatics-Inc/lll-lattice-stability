import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Topology.Order.Compact
import Mathlib.Topology.Order.Monotone
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Push
import GSTMap.Modular

/-!
# Cusp lemma (arithmetic step) and the marginal-point lemma

`cusp_not_stable`: the analytic input of Lemma C (PROOFS.md, section 4) is that for every `τ` with `Im τ ≥ 3`
the symbol at the Bloch point `(s*, 0)` has the form `c² (a² - b²)` with `c = √(Im τ)/π`, `|a - 2/5| ≤ 1/100` and
`b ≥ 12/25`.  (Arb gives, at `t₂ = 3`: `|a - 2/5| ≤ 0.009069`, `b ≥ 0.48379`, and both bounds improve with `t₂`.)
Here we only do the arithmetic that turns this into instability.

`boundary_exists`: the marginal-point lemma used for the neutral curve `Γ`: along a segment, if `f ≥ 0` at the
left end and `f < 0` at the right end for a continuous `f` (`f(τ) = min_q Σ_τ(q)` is continuous in `τ`, a
hypothesis), there is a last stable point after which everything is unstable.
-/

open UpperHalfPlane

namespace GSTMap

/-- The analytic input of the cusp lemma (Arb-verified; see `gstmap/PROOFS.md` section 4
and `gstmap/lean/README.md` for the external-input scope). -/
def CuspData (Sig : ℍ → ℝ → ℝ → ℝ) : Prop :=
  ∀ τ : ℍ, 3 ≤ τ.im → ∃ (a b s ξ : ℝ), |a - 2 / 5| ≤ 1 / 100 ∧ 12 / 25 ≤ b ∧
    Sig τ s ξ = (Real.sqrt τ.im / Real.pi) ^ 2 * (a ^ 2 - b ^ 2)

theorem cusp_arith {a b : ℝ} (ha : |a - 2 / 5| ≤ 1 / 100) (hb : 12 / 25 ≤ b) : a ^ 2 - b ^ 2 < 0 := by
  obtain ⟨h1, h2⟩ := abs_le.mp ha
  nlinarith

theorem cusp_not_stable {Sig : ℍ → ℝ → ℝ → ℝ} (h : CuspData Sig) :
    ∀ τ : ℍ, 3 ≤ τ.im → ¬ Stable Sig τ := by
  intro τ hτ hst
  obtain ⟨a, b, s, ξ, ha, hb, hS⟩ := h τ hτ
  have hneg := cusp_arith ha hb
  have hpos : 0 < (Real.sqrt τ.im / Real.pi) ^ 2 := by
    have : 0 < τ.im := τ.im_pos
    positivity
  have := hst s ξ
  rw [hS] at this
  nlinarith

/-- Marginal points: a continuous `f` which is `≥ 0` at `a` and `< 0` at `b` has a last nonnegative point `t`,
after which `f < 0` up to `b`. -/
theorem boundary_exists {f : ℝ → ℝ} (hf : Continuous f) {a b : ℝ} (hab : a < b) (ha : 0 ≤ f a) (hb : f b < 0) :
    ∃ t ∈ Set.Ico a b, 0 ≤ f t ∧ ∀ u ∈ Set.Ioc t b, f u < 0 := by
  set S : Set ℝ := {t | t ∈ Set.Icc a b ∧ 0 ≤ f t} with hSdef
  have hSclosed : IsClosed S := by
    have : S = Set.Icc a b ∩ f ⁻¹' Set.Ici 0 := by ext t; simp [hSdef]
    rw [this]
    exact isClosed_Icc.inter (isClosed_Ici.preimage hf)
  have hne : S.Nonempty := ⟨a, ⟨le_refl a, hab.le⟩, ha⟩
  have hbdd : BddAbove S := ⟨b, fun t ht => ht.1.2⟩
  have hmem : sSup S ∈ S := hSclosed.csSup_mem hne hbdd
  refine ⟨sSup S, ⟨hmem.1.1, ?_⟩, hmem.2, ?_⟩
  · rcases lt_or_eq_of_le hmem.1.2 with h | h
    · exact h
    · exfalso; have := hmem.2; rw [h] at this; linarith
  · intro u hu
    by_contra hcon
    push Not at hcon
    have huS : u ∈ S := ⟨⟨by linarith [hmem.1.1, hu.1], hu.2⟩, hcon⟩
    have := le_csSup hbdd huS
    linarith [hu.1]

end GSTMap
