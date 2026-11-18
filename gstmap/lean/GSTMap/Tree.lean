import Mathlib.Tactic.NormNum
import GSTMap.Modular

/-!
# Aggregation of the finite certificate

The adaptive quadtree computed by `tree.py` is a finite combinatorial object: every leaf is a closed
rational square `[x-h,x+h] × [y-h,y+h]` in the `τ = t₁ + i t₂` plane with one of four labels

* `skip`  – the farthest corner lies strictly inside the unit disc (checked here, exact rationals);
* `stab`  – Arb certificate: `Σ_τ > 0` on the Bloch torus minus `0` for every `τ` in the leaf;
* `unst`  – Arb certificate: a Bloch point with `Σ_τ < 0` for every `τ` in the leaf;
* `undec` – no certificate (half-width below the resolution limit).

The arithmetic content of `stab`/`unst` is hypothesis `hS`/`hU`; the cusp `Im τ ≥ 3` is hypothesis `hCusp`.
What is proved here: the leaves cover `F⁺ ∩ {Im τ ≤ 3}`, together with modular reduction (`Modular.lean`).
-/

open UpperHalfPlane

namespace GSTMap

structure RCell where
  x : ℚ
  y : ℚ
  h : ℚ
  deriving DecidableEq, Repr

def RCell.Contains (c : RCell) (p : ℝ × ℝ) : Prop :=
  (c.x : ℝ) - c.h ≤ p.1 ∧ p.1 ≤ c.x + c.h ∧ (c.y : ℝ) - c.h ≤ p.2 ∧ p.2 ≤ c.y + c.h

def RCell.child (c : RCell) (sx sy : Bool) : RCell :=
  ⟨c.x + (if sx then c.h / 2 else -(c.h / 2)), c.y + (if sy then c.h / 2 else -(c.h / 2)), c.h / 2⟩

/-- Four closed children cover the parent (boundaries included). -/
theorem RCell.children_cover (c : RCell) (p : ℝ × ℝ) (hp : c.Contains p) :
    ∃ sx sy, (c.child sx sy).Contains p := by
  obtain ⟨hx0, hx1, hy0, hy1⟩ := hp
  by_cases hx : p.1 ≤ c.x <;> by_cases hy : p.2 ≤ c.y
  · refine ⟨false, false, ?_⟩
    simp only [RCell.child, RCell.Contains, Bool.false_eq_true, ↓reduceIte]
    push_cast
    refine ⟨?_, ?_, ?_, ?_⟩ <;> linarith
  · refine ⟨false, true, ?_⟩
    simp only [RCell.child, RCell.Contains, Bool.false_eq_true, ↓reduceIte]
    push_cast
    refine ⟨?_, ?_, ?_, ?_⟩ <;> linarith
  · refine ⟨true, false, ?_⟩
    simp only [RCell.child, RCell.Contains, Bool.false_eq_true, ↓reduceIte]
    push_cast
    refine ⟨?_, ?_, ?_, ?_⟩ <;> linarith
  · refine ⟨true, true, ?_⟩
    simp only [RCell.child, RCell.Contains, ↓reduceIte]
    push_cast
    refine ⟨?_, ?_, ?_, ?_⟩ <;> linarith

inductive Kind
  | skip | stab | unst | undec
  deriving DecidableEq, Repr

inductive QT
  | leaf (k : Kind)
  | split (ll lu ul uu : QT)
  deriving Repr

/-- Leaves (kind, closed cell) of a quadtree rooted at `c`. -/
def QT.leaves : QT → RCell → List (Kind × RCell)
  | .leaf k, c => [(k, c)]
  | .split ll lu ul uu, c =>
      ll.leaves (c.child false false) ++ lu.leaves (c.child false true) ++
        ul.leaves (c.child true false) ++ uu.leaves (c.child true true)

/-- Structural counting avoids materializing rational leaf cells during kernel reduction. -/
def QT.count (kind : Kind) : QT → ℕ
  | .leaf k => if k = kind then 1 else 0
  | .split a b c d => a.count kind + b.count kind + c.count kind + d.count kind

theorem QT.count_eq (t : QT) (kind : Kind) (c : RCell) :
    ((t.leaves c).filter (fun kc => kc.1 == kind)).length = t.count kind := by
  induction t generalizing c with
  | leaf k => cases k <;> cases kind <;> simp [QT.leaves, QT.count]
  | split a b d e ia ib id ie =>
    simp [QT.leaves, QT.count, List.filter_append, ia, ib, id, ie, Nat.add_assoc]

/-- Only skip leaves require rational arithmetic; all other leaves return `true`. -/
def QT.skipOK : QT → RCell → Bool
  | .leaf k, c =>
      if k = Kind.skip then decide ((c.x+c.h)^2+(c.y+c.h)^2 < 1) else true
  | .split a b d e, c =>
      a.skipOK (c.child false false) && b.skipOK (c.child false true) &&
        d.skipOK (c.child true false) && e.skipOK (c.child true true)

theorem QT.skipOK_correct (t : QT) (c : RCell) (h : t.skipOK c = true) :
    ∀ kc ∈ t.leaves c, kc.1 = Kind.skip →
      (kc.2.x+kc.2.h)^2+(kc.2.y+kc.2.h)^2 < 1 := by
  induction t generalizing c with
  | leaf k =>
    intro kc hk hks
    have he : kc = (k,c) := by simpa [QT.leaves] using hk
    subst kc
    cases k <;> simp_all [QT.skipOK]
  | split a b d e ia ib id ie =>
    simp only [QT.skipOK, Bool.and_eq_true] at h
    rcases h with ⟨⟨⟨ha,hb⟩,hd⟩,he⟩
    intro kc hk hks
    simp only [QT.leaves, List.mem_append] at hk
    rcases hk with ((hk | hk) | hk) | hk
    · exact ia _ ha kc hk hks
    · exact ib _ hb kc hk hks
    · exact id _ hd kc hk hks
    · exact ie _ he kc hk hks

/-- **Coverage of a quadtree**: every point of the root cell lies in a leaf cell (no cell-count identity needed). -/
theorem QT.covers (t : QT) (c : RCell) (p : ℝ × ℝ) (hp : c.Contains p) :
    ∃ kc ∈ t.leaves c, kc.2.Contains p := by
  induction t generalizing c with
  | leaf k => exact ⟨(k, c), by simp [QT.leaves], hp⟩
  | split ll lu ul uu ihll ihlu ihul ihuu =>
    obtain ⟨sx, sy, hch⟩ := RCell.children_cover c p hp
    cases sx <;> cases sy
    · obtain ⟨kc, hk, hc⟩ := ihll _ hch
      exact ⟨kc, by simp [QT.leaves, hk], hc⟩
    · obtain ⟨kc, hk, hc⟩ := ihlu _ hch
      exact ⟨kc, by simp [QT.leaves, hk], hc⟩
    · obtain ⟨kc, hk, hc⟩ := ihul _ hch
      exact ⟨kc, by simp [QT.leaves, hk], hc⟩
    · obtain ⟨kc, hk, hc⟩ := ihuu _ hch
      exact ⟨kc, by simp [QT.leaves, hk], hc⟩

/-- Root squares: columns `x = 1/16 + k/8` (`k<4`), rows `y = (7+j)/8` (`j<19`), half-width `1/16`. -/
def rootCell (k j : ℕ) : RCell := ⟨1 / 16 + (k : ℚ) / 8, (7 + (j : ℚ)) / 8, 1 / 16⟩

theorem rootCells_cover (p : ℝ × ℝ) (h0 : 0 ≤ p.1) (h1 : p.1 ≤ 1 / 2)
    (h2 : 13 / 16 ≤ p.2) (h3 : p.2 ≤ 51 / 16) :
    ∃ k : ℕ, k < 4 ∧ ∃ j : ℕ, j < 19 ∧ (rootCell k j).Contains p := by
  -- column index
  have hk : ∃ k : ℕ, k < 4 ∧ (k : ℝ) / 8 ≤ p.1 ∧ p.1 ≤ ((k : ℝ) + 1) / 8 := by
    by_cases hc : (⌊8 * p.1⌋₊ : ℕ) < 4
    · refine ⟨⌊8 * p.1⌋₊, hc, ?_, ?_⟩
      · have := Nat.floor_le (by linarith : (0 : ℝ) ≤ 8 * p.1); linarith
      · have := Nat.lt_floor_add_one (8 * p.1); linarith
    · refine ⟨3, by norm_num, ?_, ?_⟩
      · have h4 : (4 : ℝ) ≤ (⌊8 * p.1⌋₊ : ℝ) := by exact_mod_cast (not_lt.mp hc)
        have := Nat.floor_le (by linarith : (0 : ℝ) ≤ 8 * p.1)
        norm_num; linarith
      · norm_num; linarith
  have hj : ∃ j : ℕ, j < 19 ∧ (13 + 2 * (j : ℝ)) / 16 ≤ p.2 ∧ p.2 ≤ (15 + 2 * (j : ℝ)) / 16 := by
    by_cases hc : (⌊(16 * p.2 - 13) / 2⌋₊ : ℕ) < 19
    · refine ⟨⌊(16 * p.2 - 13) / 2⌋₊, hc, ?_, ?_⟩
      · have := Nat.floor_le (by linarith : (0 : ℝ) ≤ (16 * p.2 - 13) / 2); linarith
      · have := Nat.lt_floor_add_one ((16 * p.2 - 13) / 2); linarith
    · refine ⟨18, by norm_num, ?_, ?_⟩
      · have h19 : (19 : ℝ) ≤ (⌊(16 * p.2 - 13) / 2⌋₊ : ℝ) := by exact_mod_cast (not_lt.mp hc)
        have := Nat.floor_le (by linarith : (0 : ℝ) ≤ (16 * p.2 - 13) / 2)
        norm_num; linarith
      · norm_num; linarith
  obtain ⟨k, hk4, hk0, hk1⟩ := hk
  obtain ⟨j, hj19, hj0, hj1⟩ := hj
  refine ⟨k, hk4, j, hj19, ?_⟩
  simp only [rootCell, RCell.Contains]
  push_cast
  refine ⟨?_, ?_, ?_, ?_⟩ <;> linarith

/-- All leaves of a forest indexed by the 76 root squares. -/
def allLeaves (forest : ℕ → ℕ → QT) : List (Kind × RCell) :=
  (List.range 4).flatMap fun k => (List.range 19).flatMap fun j => (forest k j).leaves (rootCell k j)

def forestCount (forest : ℕ → ℕ → QT) (kind : Kind) : ℕ :=
  ((List.range 4).map fun k => ((List.range 19).map fun j => (forest k j).count kind).sum).sum

theorem allLeaves_count (forest : ℕ → ℕ → QT) (kind : Kind) :
    ((allLeaves forest).filter (fun kc => kc.1 == kind)).length = forestCount forest kind := by
  simp [allLeaves, forestCount, List.filter_flatMap, List.length_flatMap, QT.count_eq]

theorem allLeaves_cover (forest : ℕ → ℕ → QT) (p : ℝ × ℝ) (h0 : 0 ≤ p.1) (h1 : p.1 ≤ 1 / 2)
    (h2 : 13 / 16 ≤ p.2) (h3 : p.2 ≤ 51 / 16) :
    ∃ kc ∈ allLeaves forest, kc.2.Contains p := by
  obtain ⟨k, hk, j, hj, hp⟩ := rootCells_cover p h0 h1 h2 h3
  obtain ⟨kc, hkc, hc⟩ := (forest k j).covers (rootCell k j) p hp
  refine ⟨kc, ?_, hc⟩
  simp only [allLeaves, List.mem_flatMap, List.mem_range]
  exact ⟨k, hk, j, hj, hkc⟩

/-- **Global classification, assuming the leaf certificates.**  For every `τ ∈ ℍ` there is a representative
`τ' ∈ F⁺` with the same stability, and `τ'` is either in a stable leaf, in an unstable leaf (or in the cusp
`Im ≥ 3`), or in an undecided leaf. -/
theorem classify {Sig : ℍ → ℝ → ℝ → ℝ} (hcov : Covariant Sig) (forest : ℕ → ℕ → QT)
    (hskip : ∀ kc ∈ allLeaves forest, kc.1 = Kind.skip →
      (kc.2.x + kc.2.h) ^ 2 + (kc.2.y + kc.2.h) ^ 2 < 1)
    (hS : ∀ kc ∈ allLeaves forest, kc.1 = Kind.stab →
      ∀ τ : ℍ, kc.2.Contains (τ.re, τ.im) → Stable Sig τ)
    (hU : ∀ kc ∈ allLeaves forest, kc.1 = Kind.unst →
      ∀ τ : ℍ, kc.2.Contains (τ.re, τ.im) → ¬ Stable Sig τ)
    (hCusp : ∀ τ : ℍ, 3 ≤ τ.im → ¬ Stable Sig τ) (τ : ℍ) :
    ∃ τ' : ℍ, FPlus τ' ∧ (Stable Sig τ ↔ Stable Sig τ') ∧
      ((Stable Sig τ' ∧ ∃ kc ∈ allLeaves forest, kc.1 = Kind.stab ∧ kc.2.Contains (τ'.re, τ'.im)) ∨
       (¬ Stable Sig τ' ∧ (3 ≤ τ'.im ∨
          ∃ kc ∈ allLeaves forest, kc.1 = Kind.unst ∧ kc.2.Contains (τ'.re, τ'.im))) ∨
       (∃ kc ∈ allLeaves forest, kc.1 = Kind.undec ∧ kc.2.Contains (τ'.re, τ'.im))) := by
  obtain ⟨τ', hF, hiff⟩ := reduce_to_FPlus hcov τ
  refine ⟨τ', hF, hiff, ?_⟩
  obtain ⟨h0, h1, hn⟩ := hF
  have him : 0 < τ'.im := τ'.im_pos
  by_cases h3 : 3 ≤ τ'.im
  · exact Or.inr (Or.inl ⟨hCusp τ' h3, Or.inl h3⟩)
  · push_neg at h3
    have h2 : (13 : ℝ) / 16 ≤ τ'.im := by
      by_contra hlt
      push_neg at hlt
      nlinarith
    obtain ⟨kc, hkc, hcont⟩ := allLeaves_cover forest (τ'.re, τ'.im) h0 h1 h2 (by simp only; linarith)
    obtain ⟨k, c⟩ := kc
    cases k with
    | skip =>
      exfalso
      have hs := hskip (Kind.skip, c) hkc rfl
      obtain ⟨_, hx1, _, hy1⟩ := hcont
      simp only at hx1 hy1
      have hs' : ((c.x : ℝ) + c.h) ^ 2 + ((c.y : ℝ) + c.h) ^ 2 < 1 := by exact_mod_cast hs
      nlinarith
    | stab => exact Or.inl ⟨hS (Kind.stab, c) hkc rfl τ' hcont, (Kind.stab, c), hkc, rfl, hcont⟩
    | unst => exact Or.inr (Or.inl ⟨hU (Kind.unst, c) hkc rfl τ' hcont, Or.inr ⟨(Kind.unst, c), hkc, rfl, hcont⟩⟩)
    | undec => exact Or.inr (Or.inr ⟨(Kind.undec, c), hkc, rfl, hcont⟩)

end GSTMap
