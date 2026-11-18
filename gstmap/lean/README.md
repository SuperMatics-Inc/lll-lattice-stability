# GSTMap: conditional partial-map coverage

`GSTMap` checks the partial shape-map argument. The two individual lattice shapes are treated in the separate `GSTStability` project.

```sh
lake build
```

The project uses Lean 4.35.0 and Mathlib v4.35.0 at revision `56a005855ea17e0b0958f73ec801ae140b0a1679`. `lake-manifest.json` records the transitive dependencies. GitHub Actions builds the project with a 4096 MB per-process limit.

`GSTMap.main_classification` is a combinatorial coverage/reduction theorem for an abstract symbol `Sig`. It assumes:
- `Covariant`: modular/reflection covariance;
- `hS`: strict positivity on each closed S leaf away from zero;
- `hU`: a negative Bloch witness on each closed U leaf;
- `CuspData`: the analytic cusp bounds.

Lean checks coverage, skip geometry and counts for 76 roots and 40570 leaves (K66/S4178/U25088/X11238). The Gaussian symbol, theta derivatives, covariance, interval enclosures and tails are analytic inputs to the map argument. The curve `Gamma4` is checked separately. `boundary_exists` assumes continuity of an abstract function. The quartic and sextic identities follow by ring algebra from the supplied lower-order identities.

The eight axiom lists contain only `propext`, `Classical.choice` and `Quot.sound`. The proofs use no custom axiom, `sorry` or `native_decide`.

The reproduction command requires the generated `Data.lean` to match the supplied LF bytes. `SHA256SUMS` lists the repository file checksums.
