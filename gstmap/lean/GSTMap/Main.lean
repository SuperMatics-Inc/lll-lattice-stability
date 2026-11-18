import GSTMap.Modular
import GSTMap.Tree
import GSTMap.Cusp
import GSTMap.Data
import GSTMap.Quartic

/-!
Theorem: conditional partial-map coverage and modular reduction

Hypotheses (all external to Lean; see `gstmap/PROOFS.md` sections 3–4
and `gstmap/lean/README.md`):
* `hcov`  : the three covariances `T`, `S`, `R` of the symbol (proved on paper, PROOFS.md section 1);
* `hS`    : each `stab` leaf of the generated forest is an Arb-certified stable box;
* `hU`    : each `unst` leaf of the generated forest is an Arb-certified unstable box;
* `hC`    : the analytic cusp data (`Im τ ≥ 3`, PROOFS.md section 4, Arb-verified at `t₂ = 3`).
Proved here: the leaves of `Data.forest` (with the exact rational data checks `skip_ok`, `count_*`) cover
`F⁺ ∩ {Im τ ≤ 3}`, and every `τ ∈ ℍ` is reduced into `F⁺`.
-/

open UpperHalfPlane

namespace GSTMap

theorem main_classification {Sig : ℍ → ℝ → ℝ → ℝ} (hcov : Covariant Sig)
    (hS : ∀ kc ∈ allLeaves forest, kc.1 = Kind.stab → ∀ τ : ℍ, kc.2.Contains (τ.re, τ.im) → Stable Sig τ)
    (hU : ∀ kc ∈ allLeaves forest, kc.1 = Kind.unst → ∀ τ : ℍ, kc.2.Contains (τ.re, τ.im) → ¬ Stable Sig τ)
    (hC : CuspData Sig) (τ : ℍ) :
    ∃ τ' : ℍ, FPlus τ' ∧ (Stable Sig τ ↔ Stable Sig τ') ∧
      ((Stable Sig τ' ∧ ∃ kc ∈ allLeaves forest, kc.1 = Kind.stab ∧ kc.2.Contains (τ'.re, τ'.im)) ∨
       (¬ Stable Sig τ' ∧ (3 ≤ τ'.im ∨
          ∃ kc ∈ allLeaves forest, kc.1 = Kind.unst ∧ kc.2.Contains (τ'.re, τ'.im))) ∨
       (∃ kc ∈ allLeaves forest, kc.1 = Kind.undec ∧ kc.2.Contains (τ'.re, τ'.im))) :=
  classify hcov forest skip_ok hS hU (cusp_not_stable hC) τ

#print axioms main_classification
#print axioms skip_ok
#print axioms cusp_not_stable
#print axioms boundary_exists
#print axioms reduce_to_FPlus
#print axioms quartic_coefficient
#print axioms sextic_coefficient
#print axioms quartic_remainder_sign

end GSTMap
