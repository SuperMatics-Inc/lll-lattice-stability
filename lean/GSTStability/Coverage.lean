import GSTStability.Remainder

noncomputable section

namespace GST

abbrev Point := ℝ × ℝ

def radius (p : Point) : ℝ := max |p.1| |p.2|

structure Cell where
  x : ℝ
  y : ℝ
  h : ℝ

def Cell.Contains (c : Cell) (p : Point) : Prop :=
  c.x - c.h ≤ p.1 ∧ p.1 ≤ c.x + c.h ∧
  c.y - c.h ≤ p.2 ∧ p.2 ≤ c.y + c.h

def Cell.child (c : Cell) (sx sy : Bool) : Cell :=
  ⟨c.x + (if sx then c.h / 2 else -(c.h / 2)),
   c.y + (if sy then c.h / 2 else -(c.h / 2)), c.h / 2⟩

/-- Four *closed* children cover the parent, including subdivision boundaries. -/
theorem children_cover (c : Cell) (p : Point) (hp : c.Contains p) :
    ∃ sx sy, (c.child sx sy).Contains p := by
  rcases hp with ⟨hx0, hx1, hy0, hy1⟩
  by_cases hx : p.1 ≤ c.x <;> by_cases hy : p.2 ≤ c.y
  · refine ⟨false, false, ?_⟩
    simp only [Cell.child, Cell.Contains, Bool.false_eq_true, ↓reduceIte]
    constructor <;> try linarith
    constructor <;> try linarith
    constructor <;> linarith
  · refine ⟨false, true, ?_⟩
    simp only [Cell.child, Cell.Contains, Bool.false_eq_true, ↓reduceIte]
    constructor <;> try linarith
    constructor <;> try linarith
    constructor <;> linarith
  · refine ⟨true, false, ?_⟩
    simp only [Cell.child, Cell.Contains, Bool.false_eq_true, ↓reduceIte]
    constructor <;> try linarith
    constructor <;> try linarith
    constructor <;> linarith
  · refine ⟨true, true, ?_⟩
    simp only [Cell.child, Cell.Contains, ↓reduceIte]
    constructor <;> try linarith
    constructor <;> try linarith
    constructor <;> linarith

/-- A terminating certificate tree. Leaves are either near-square inclusions
or validated far cells. The analytic interval enclosure is leaf evidence. -/
inductive CertifiedTree (r : ℝ) (sigma : Point → ℝ) : Cell → Prop
  | near (c : Cell) (h : ∀ p, c.Contains p → radius p ≤ r) : CertifiedTree r sigma c
  | far (c : Cell) (h : ∀ p, c.Contains p → 0 < sigma p) : CertifiedTree r sigma c
  | split (c : Cell) (h : ∀ sx sy, CertifiedTree r sigma (c.child sx sy)) :
      CertifiedTree r sigma c

/-- Structural induction proves coverage for arbitrary adaptive depth and
arbitrary choices of which cells split. No cell-count identity is substituted
for a coverage proof. -/
theorem adaptive_exhaustion (r : ℝ) (sigma : Point → ℝ) (c : Cell)
    (tree : CertifiedTree r sigma c) :
    ∀ p, c.Contains p → radius p ≤ r ∨ 0 < sigma p := by
  induction tree with
  | near c h => intro p hp; exact Or.inl (h p hp)
  | far c h => intro p hp; exact Or.inr (h p hp)
  | split c _ ih =>
    intro p hp
    obtain ⟨sx, sy, hchild⟩ := children_cover c p hp
    exact ih sx sy p hchild

/-- Root squares used verbatim by cert_sigma.py. -/
def root (i j : ℕ) : Cell := ⟨1 / 8 + (i : ℝ) / 4, -3 / 8 + (j : ℝ) / 4, 1 / 8⟩

theorem eight_roots_cover (p : Point)
    (hx0 : 0 ≤ p.1) (hx1 : p.1 ≤ 1 / 2)
    (hy0 : -(1 / 2) ≤ p.2) (hy1 : p.2 ≤ 1 / 2) :
    ∃ i < 2, ∃ j < 4, (root i j).Contains p := by
  have hxi : ∃ i < 2,
      (root i 0).x - (1 / 8 : ℝ) ≤ p.1 ∧ p.1 ≤ (root i 0).x + 1 / 8 := by
    by_cases hx : p.1 ≤ 1 / 4
    · refine ⟨0, by norm_num, ?_⟩; norm_num [root]; constructor <;> linarith
    · refine ⟨1, by norm_num, ?_⟩; norm_num [root]; constructor <;> linarith
  have hyj : ∃ j < 4,
      (root 0 j).y - (1 / 8 : ℝ) ≤ p.2 ∧ p.2 ≤ (root 0 j).y + 1 / 8 := by
    by_cases h0 : p.2 ≤ -(1 / 4)
    · refine ⟨0, by norm_num, ?_⟩; norm_num [root]; constructor <;> linarith
    · by_cases h1 : p.2 ≤ 0
      · refine ⟨1, by norm_num, ?_⟩; norm_num [root]; constructor <;> linarith
      · by_cases h2 : p.2 ≤ 1 / 4
        · refine ⟨2, by norm_num, ?_⟩; norm_num [root]; constructor <;> linarith
        · refine ⟨3, by norm_num, ?_⟩; norm_num [root]; constructor <;> linarith
  obtain ⟨i, hi, hx⟩ := hxi
  obtain ⟨j, hj, hy⟩ := hyj
  exact ⟨i, hi, j, hj, hx.1, hx.2, hy.1, hy.2⟩

theorem half_torus_covered (r : ℝ) (sigma : Point → ℝ)
    (trees : ∀ i < 2, ∀ j < 4, CertifiedTree r sigma (root i j))
    (p : Point) (hx0 : 0 ≤ p.1) (hx1 : p.1 ≤ 1 / 2)
    (hy0 : -(1 / 2) ≤ p.2) (hy1 : p.2 ≤ 1 / 2) :
    radius p ≤ r ∨ 0 < sigma p := by
  obtain ⟨i, hi, j, hj, hp⟩ := eight_roots_cover p hx0 hx1 hy0 hy1
  exact adaptive_exhaustion r sigma (root i j) (trees i hi j hj) p hp

/-- The box omission test in the program is sufficient for every point of a
skipped cell to belong to the near square. -/
theorem omitted_cell_inside (c : Cell) (r : ℝ)
    (hx : |c.x| + c.h ≤ r) (hy : |c.y| + c.h ≤ r)
    (p : Point) (hp : c.Contains p) : radius p ≤ r := by
  rcases hp with ⟨hx0, hx1, hy0, hy1⟩
  apply max_le
  · apply abs_le.mpr
    constructor <;> have hn := neg_le_abs c.x <;>
      have hp' := le_abs_self c.x <;> linarith
  · apply abs_le.mpr
    constructor <;> have hn := neg_le_abs c.y <;>
      have hp' := le_abs_self c.y <;> linarith

end GST
