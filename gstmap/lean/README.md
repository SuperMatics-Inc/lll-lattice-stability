# GSTMap: conditional partial-map coverage

This project is **separate from** the paper's original `GSTStability` project and its `GST.*` theorems. Never merge their scopes or leaf counts.

```sh
lake build
```

Lean is pinned to 4.35.0; Mathlib to tag `v4.35.0`, resolving to `56a005855ea17e0b0958f73ec801ae140b0a1679`, with the transitive git pins in `lake-manifest.json`. The package requirement and manifest use the same input revision. The portable dependency is git-pinned. CI builds this public configuration, using focused Mathlib imports and the unchanged `-M4096` per-process limit. Earlier local builds used the identical shared Mathlib cache through a local path requirement. The theorem statements and proof bodies are unchanged.

`GSTMap.main_classification` is a combinatorial coverage/reduction theorem for an abstract symbol `Sig`. It assumes:
- `Covariant`: modular/reflection covariance;
- `hS`: strict positivity on each closed S leaf away from zero;
- `hU`: a negative Bloch witness on each closed U leaf;
- `CuspData`: the analytic cusp bounds.

The kernel checks coverage, skip geometry and counts for 76 roots and 40570 leaves (K66/S4178/U25088/X11238). It does **not** check the Gaussian symbol, theta differentiation, covariance, interval arithmetic, tails or Gamma4. `boundary_exists` is generic and assumes continuity; quartic/sextic identities are generic ring calculations with the lower-order relations as hypotheses. No PDE continuity theorem is instantiated here.

All eight axiom lists use only `propext`, `Classical.choice`, `Quot.sound`. No custom axiom, `sorry` or `native_decide` is used. Captured local-build evidence is retained privately.

The replay requires generated `Data.lean` to match the shipped source byte-for-byte in an LF checkout. `SHA256SUMS` at the repository root pins all released files.
