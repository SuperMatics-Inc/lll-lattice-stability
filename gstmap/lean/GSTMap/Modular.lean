import Mathlib.NumberTheory.Modular
import Mathlib.LinearAlgebra.Matrix.FixedDetMatrices
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Ring

/-!
# Modular reduction of the GST stability problem

`Sig τ s ξ` is the determinant symbol `Σ_τ(s,ξ)`.  The three covariances `T`, `S`, `R` are the analytic
input (proved in `PROOFS.md` by 2D Poisson summation, and checked in float64 / Arb); they are hypotheses
here.  Everything in this file is then a theorem of Mathlib's modular group API:
spectral stability of `τ` depends only on the orbit of `τ` under `⟨T,S,R⟩`, so it suffices to
classify `F⁺ = {0 ≤ Re τ ≤ 1/2, |τ| ≥ 1}`.
-/

open UpperHalfPlane ModularGroup MatrixGroups Subgroup

namespace GSTMap

/-- The orientation reversing reflection `τ ↦ -conj τ` of the upper half plane. -/
def refl (τ : ℍ) : ℍ := ⟨-(starRingEnd ℂ) (τ : ℂ), by simpa using τ.im_pos⟩

@[simp] lemma refl_re (τ : ℍ) : (refl τ).re = -τ.re := by
  show (-(starRingEnd ℂ) (τ : ℂ)).re = -(τ : ℂ).re
  simp

@[simp] lemma refl_im (τ : ℍ) : (refl τ).im = τ.im := by
  show (-(starRingEnd ℂ) (τ : ℂ)).im = (τ : ℂ).im
  simp

/-- `T • τ = τ + 1` and `S • τ = -1/τ` (sanity: the abstract `SL(2,ℤ)` generators are the usual maps). -/
lemma T_smul_coe (τ : ℍ) : ((ModularGroup.T • τ : ℍ) : ℂ) = (τ : ℂ) + 1 := by
  rw [modular_T_smul]; simp [add_comm]

lemma S_smul_coe (τ : ℍ) : ((ModularGroup.S • τ : ℍ) : ℂ) = -((τ : ℂ))⁻¹ := by
  rw [modular_S_smul]; simp

/-- Covariance of the determinant symbol (analytic input, hypotheses of every theorem below). -/
structure Covariant (Sig : ℍ → ℝ → ℝ → ℝ) : Prop where
  T : ∀ τ s ξ, Sig (ModularGroup.T • τ) s ξ = Sig τ s (ξ + s)
  S : ∀ τ s ξ, Sig (ModularGroup.S • τ) s ξ = Sig τ ξ (-s)
  R : ∀ τ s ξ, Sig (refl τ) s ξ = Sig τ s (-ξ)

/-- Spectral stability: the symbol is nonnegative on the whole Bloch torus (real `(s,ξ)` suffice, by periodicity). -/
def Stable (Sig : ℍ → ℝ → ℝ → ℝ) (τ : ℍ) : Prop := ∀ s ξ, 0 ≤ Sig τ s ξ

theorem stable_T {Sig : ℍ → ℝ → ℝ → ℝ} (h : Covariant Sig) (τ : ℍ) :
    Stable Sig (ModularGroup.T • τ) ↔ Stable Sig τ := by
  constructor
  · intro hs s ξ
    have := hs s (ξ - s)
    rwa [h.T, sub_add_cancel] at this
  · intro hs s ξ
    rw [h.T]; exact hs s _

theorem stable_S {Sig : ℍ → ℝ → ℝ → ℝ} (h : Covariant Sig) (τ : ℍ) :
    Stable Sig (ModularGroup.S • τ) ↔ Stable Sig τ := by
  constructor
  · intro hs s ξ
    have := hs (-ξ) s
    rwa [h.S, neg_neg] at this
  · intro hs s ξ
    rw [h.S]; exact hs _ _

theorem stable_R {Sig : ℍ → ℝ → ℝ → ℝ} (h : Covariant Sig) (τ : ℍ) :
    Stable Sig (refl τ) ↔ Stable Sig τ := by
  constructor
  · intro hs s ξ
    have := hs s (-ξ)
    rwa [h.R, neg_neg] at this
  · intro hs s ξ
    rw [h.R]; exact hs _ _

theorem stable_smul {Sig : ℍ → ℝ → ℝ → ℝ} (h : Covariant Sig) (g : SL(2, ℤ)) (τ : ℍ) :
    Stable Sig (g • τ) ↔ Stable Sig τ := by
  have hg : g ∈ Subgroup.closure {ModularGroup.S, ModularGroup.T} := by
    simp [SpecialLinearGroup.SL2Z_generators]
  revert τ
  refine Subgroup.closure_induction (p := fun g _ => ∀ τ : ℍ, Stable Sig (g • τ) ↔ Stable Sig τ)
    ?_ ?_ ?_ ?_ hg
  · intro x hx τ
    rcases hx with rfl | hx
    · exact stable_S h τ
    · rw [Set.mem_singleton_iff] at hx
      subst hx
      exact stable_T h τ
  · intro τ; simp
  · intro x y _ _ ihx ihy τ
    rw [mul_smul, ihx, ihy]
  · intro x _ ih τ
    have := ih (x⁻¹ • τ)
    rw [smul_inv_smul] at this
    exact this.symm

/-- Points of `F⁺ = {0 ≤ Re τ ≤ 1/2, |τ| ≥ 1}`. -/
def FPlus (τ : ℍ) : Prop := 0 ≤ τ.re ∧ τ.re ≤ 1 / 2 ∧ 1 ≤ τ.re ^ 2 + τ.im ^ 2

/-- **Modular reduction.**  Every `τ ∈ ℍ` has the same stability as some point of `F⁺`. -/
theorem reduce_to_FPlus {Sig : ℍ → ℝ → ℝ → ℝ} (h : Covariant Sig) (τ : ℍ) :
    ∃ τ' : ℍ, FPlus τ' ∧ (Stable Sig τ ↔ Stable Sig τ') := by
  obtain ⟨g, hg⟩ := ModularGroup.exists_smul_mem_fd τ
  have hD : 1 ≤ Complex.normSq ((g • τ : ℍ) : ℂ) ∧ |(g • τ).re| ≤ 1 / 2 := hg
  have hn : Complex.normSq ((g • τ : ℍ) : ℂ) = (g • τ).re ^ 2 + (g • τ).im ^ 2 := by
    rw [Complex.normSq_apply]
    show _ = (((g • τ : ℍ) : ℂ).re) ^ 2 + (((g • τ : ℍ) : ℂ).im) ^ 2
    ring
  rw [hn] at hD
  have hr := abs_le.mp hD.2
  by_cases hre : 0 ≤ (g • τ).re
  · exact ⟨g • τ, ⟨hre, hr.2, hD.1⟩, (stable_smul h g τ).symm⟩
  · refine ⟨refl (g • τ), ⟨?_, ?_, ?_⟩, ?_⟩
    · rw [refl_re]; linarith
    · rw [refl_re]; linarith
    · rw [refl_re, refl_im]; nlinarith [hD.1]
    · rw [stable_R h, stable_smul h]

end GSTMap
