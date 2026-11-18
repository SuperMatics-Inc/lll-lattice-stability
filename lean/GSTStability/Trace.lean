import GSTStability.Coverage

noncomputable section

namespace GST

/-- Exact rational geometry for retaining an actual execution's subdivision tree. -/
structure RationalCell where
  x : ℚ
  y : ℚ
  h : ℚ

def RationalCell.toReal (c : RationalCell) : Cell := ⟨c.x, c.y, c.h⟩

def RationalCell.child (c : RationalCell) (sx sy : Bool) : RationalCell :=
  ⟨c.x + (if sx then c.h / 2 else -(c.h / 2)),
   c.y + (if sy then c.h / 2 else -(c.h / 2)), c.h / 2⟩

theorem rational_child_real (c : RationalCell) (sx sy : Bool) :
    (c.child sx sy).toReal = c.toReal.child sx sy := by
  cases sx <;> cases sy <;>
    simp [RationalCell.child, RationalCell.toReal, Cell.child]

inductive QuadTree where
  | near
  | far
  | split (ll lu ul uu : QuadTree)
  deriving Repr

/-- Number of near and far leaves, for comparison with the calculation logs. -/
def QuadTree.leafCounts : QuadTree → ℕ × ℕ
  | .near => (1, 0)
  | .far => (0, 1)
  | .split ll lu ul uu =>
      (ll.leafCounts.1 + lu.leafCounts.1 + ul.leafCounts.1 + uu.leafCounts.1,
       ll.leafCounts.2 + lu.leafCounts.2 + ul.leafCounts.2 + uu.leafCounts.2)

def rootLeafCounts (trees : ℕ → ℕ → QuadTree) : ℕ × ℕ :=
  ((List.ofFn fun k : Fin 8 => (trees (k.val / 4) (k.val % 4)).leafCounts)).foldl
    (fun a b => (a.1 + b.1, a.2 + b.2)) (0, 0)

/-- Only near-leaf omissions need numerical geometry checks: far leaves retain
an explicit analytic enclosure hypothesis. Child coordinates are exact rationals. -/
def QuadTree.geometry (r : ℚ) (c : RationalCell) : QuadTree → Bool
  | .near => decide (|c.x| + c.h ≤ r ∧ |c.y| + c.h ≤ r)
  | .far => true
  | .split ll lu ul uu =>
      (ll.geometry r (c.child false false) && lu.geometry r (c.child false true)) &&
      (ul.geometry r (c.child true false) && uu.geometry r (c.child true true))

/-- All and only accepted far-leaf inequalities. This is where the interval
certificates enter as hypotheses; no theta enclosure is an axiom. -/
def QuadTree.farEnclosures (sigma : Point → ℝ) (c : RationalCell) : QuadTree → Prop
  | .near => True
  | .far => ∀ p, c.toReal.Contains p → 0 < sigma p
  | .split ll lu ul uu =>
      (ll.farEnclosures sigma (c.child false false) ∧
       lu.farEnclosures sigma (c.child false true)) ∧
      (ul.farEnclosures sigma (c.child true false) ∧
       uu.farEnclosures sigma (c.child true true))

/-- A continuous determinant lower enclosure on every accepted far cell.
The lower endpoint is common to the finite certificate, not an assumption
of positivity itself. Per-cell Taylor/interval bounds imply this enclosure. -/
def QuadTree.farLowerBounds (lower : ℝ) (sigma : Point → ℝ)
    (c : RationalCell) : QuadTree → Prop
  | .near => True
  | .far => ∀ p, c.toReal.Contains p → lower ≤ sigma p
  | .split ll lu ul uu =>
      (ll.farLowerBounds lower sigma (c.child false false) ∧
       lu.farLowerBounds lower sigma (c.child false true)) ∧
      (ul.farLowerBounds lower sigma (c.child true false) ∧
       uu.farLowerBounds lower sigma (c.child true true))

theorem far_lower_bounds_positive (t : QuadTree) (lower : ℝ) (hl : 0 < lower)
    (sigma : Point → ℝ) (c : RationalCell)
    (hf : t.farLowerBounds lower sigma c) : t.farEnclosures sigma c := by
  induction t generalizing c with
  | near => trivial
  | far => exact fun p hp => lt_of_lt_of_le hl (hf p hp)
  | split ll lu ul uu ihll ihlu ihul ihuu =>
    rcases hf with ⟨⟨hll, hlu⟩, hul, huu⟩
    exact ⟨⟨ihll _ hll, ihlu _ hlu⟩, ihul _ hul, ihuu _ huu⟩

/-- Soundness of the execution-trace checker, including exact near-cell tests. -/
theorem checked_trace_sound (t : QuadTree) (r : ℚ) (sigma : Point → ℝ)
    (c : RationalCell) (hc : t.geometry r c = true)
    (hf : t.farEnclosures sigma c) : CertifiedTree (r : ℝ) sigma c.toReal := by
  induction t generalizing c with
  | near =>
    have hg : |c.x| + c.h ≤ r ∧ |c.y| + c.h ≤ r := of_decide_eq_true hc
    apply CertifiedTree.near
    intro p hp
    apply omitted_cell_inside c.toReal (r : ℝ) _ _ p hp
    · change |(c.x : ℝ)| + (c.h : ℝ) ≤ (r : ℝ)
      exact_mod_cast hg.1
    · change |(c.y : ℝ)| + (c.h : ℝ) ≤ (r : ℝ)
      exact_mod_cast hg.2
  | far => exact CertifiedTree.far c.toReal hf
  | split ll lu ul uu ihll ihlu ihul ihuu =>
    simp only [QuadTree.geometry, Bool.and_eq_true] at hc
    rcases hc with ⟨⟨hll, hlu⟩, hul, huu⟩
    rcases hf with ⟨⟨fll, flu⟩, ful, fuu⟩
    apply CertifiedTree.split
    intro sx sy
    cases sx <;> cases sy
    · simpa only [rational_child_real] using ihll (c.child false false) hll fll
    · simpa only [rational_child_real] using ihlu (c.child false true) hlu flu
    · simpa only [rational_child_real] using ihul (c.child true false) hul ful
    · simpa only [rational_child_real] using ihuu (c.child true true) huu fuu

def rationalRoot (i j : ℕ) : RationalCell :=
  ⟨1 / 8 + (i : ℚ) / 4, -3 / 8 + (j : ℚ) / 4, 1 / 8⟩

theorem rational_root_real (i j : ℕ) : (rationalRoot i j).toReal = root i j := by
  simp [rationalRoot, RationalCell.toReal, root]

end GST
