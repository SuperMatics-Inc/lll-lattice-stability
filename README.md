# LLL lattice spectral stability

Davis Rens, SuperMatics Inc. Software for *Certified non-hexagonal spectral stability in the lowest-Landau-level equation*, SuperMatics Technical Report SM-TR-2026-02 (2026), [10.5281/zenodo.24326107](https://doi.org/10.5281/zenodo.24326107). The report has not been peer reviewed.

Software version 1.0.0: [10.5281/zenodo.24326109](https://doi.org/10.5281/zenodo.24326109). Source: [SuperMatics-Inc/lll-lattice-stability](https://github.com/SuperMatics-Inc/lll-lattice-stability).

## Results

Interval calculations give a positive Bloch determinant everywhere on the punctured torus for two nonhexagonal lattice shapes, `tau = -2/5 + 19i/20` and `tau = 2/5 + i`. The calculations also find two simple inflection points of the hexagonal dispersion, certify a derivative sum of at least 1.20, and bound the theta-derivative expression strictly above 26.95364.

We derive the Gaussian Bloch symbol, its covariance, and its cancellation below degree four. Lean checks the algebraic consequences of the moment identities, the near-origin scaling, and the geometry of the 776- and 719-leaf covers. It takes analytic covariance and the enclosures as hypotheses. The theta evaluations and interval arithmetic remain external to Lean.

The partial shape map certifies 4178 stable closed boxes, 25088 unstable boxes, cusp instability for `Im(tau) >= 3`, and a long-wave transition curve. The separate `GSTMap` project checks the 40570-leaf partition and the conditional reduction. Its proof assumes the analytic bounds and does not check the Gaussian symbol, interval arithmetic or the curve `Gamma4`. The undecided X leaves have coordinate area `5619/134217728`. The map locates a long-wave transition, not the full spectral boundary.

The origin has a nonzero nilpotent block. Real spectrum alone gives neither uniform boundedness nor weighted, dispersive or nonlinear stability. The examples contradict a spectral reading of the expectation of Germain, Schwinte and Thomann, rather than their theorems. The sampled landscapes, S samples and finite-torus spectra are numerical calculations, not continuum bounds.

One interval implementation checks the S leaves. Separate floating-point code samples all 4178, and a second interval calculation covered 12 of them, leaving 4166 without a second interval check. One interval implementation checks `Gamma4` and is shared with the separate U verifier. Both Arb and cancellation-preserving `mpmath.iv` have been used to check U; the public map command uses only `mpmath.iv`. The analytic formulas and tail estimates are given in `gstmap/PROOFS.md`.

## Reproduction

Install Git, Python 3.12+, [uv](https://docs.astral.sh/uv/) and, for formal checks, [elan](https://github.com/leanprover/elan).

```sh
python scripts/verify_pins.py
python scripts/verify.py --output replay-output
uv run --no-project --with-requirements requirements.txt python gstmap/replay.py --evidence --workers 2 --output map-output
uv run --no-project --with-requirements requirements.txt python figures/plot_partial_map.py --output-dir figure-output
```

The first calculation command runs the Arb and floating-point suite, the separate `mpmath.iv` certificates for the two lattice shapes and the hexagonal dispersion, and exact rational cover checks. It also runs the physical-space torus calculation at `M=256` (65536 characters). The partial-map command checks the exact partition, every U leaf, all 397 curve patches and 396 joins, the negative controls, independent S samples and the cusp. The drivers record each subprocess's exit status in JSON. `--evidence` skips the expensive primary Arb S cover; remove that flag to run it. The generated Lean partition must match the supplied LF bytes.

The `partial-map` CI job checks the exact partition and the generated Lean data, all 397 `Gamma4` patches and 396 joins, the U and `Gamma4` negative controls and the cusp:

```sh
uv run --no-project --with-requirements requirements.txt python gstmap/replay.py --ci-subset --workers 4 --output map-output
```

This job leaves the U cover, S sampling and primary Arb S cover to the reproduction commands above. The figure script draws the exact partition; thin X boxes can be narrower than a print pixel. The white region below the unit arc is outside the fundamental domain.

For each formal project:

```sh
cd lean                 # separately, gstmap/lean
lake exe cache get
lake build
```

Both projects use Lean 4.35.0 and Mathlib v4.35.0 at revision `56a005855ea17e0b0958f73ec801ae140b0a1679`. The 37 axiom checks for the original project and eight for the map list only propositional extensionality, classical choice and quotient soundness. Neither project uses custom axioms or admitted proofs. GitHub Actions builds both projects. The map build uses focused Mathlib imports and a 4096 MiB per-process limit.

`SHA256SUMS` lists checksums for an LF checkout. PDFs are treated as binary files. `requirements.txt` specifies the numerical and plotting package versions, and the drivers record the Python and package versions used. Raw execution logs and the report are not part of this repository. The supplement is licensed under Apache-2.0. It contains no experimental data or third-party manuscript.

## AI use

D.R. conceived and directed the programme, selected and framed the problems, orchestrated the research and the verification, wrote and edited the manuscript extensively, and reviewed the full manuscript and its arguments in depth. He takes responsibility for the final work.

Under Davis Rens’s direction and review, AI agents generated proofs, computer-assisted certificates, code, Lean sources and manuscript drafts and edits. Other AI agents checked the results; these checks were separate from generation, not independent human derivations. The recorded model labels were Anthropic Claude Opus 5.5, Anthropic Claude Sonnet 5.5 and OpenAI GPT-6.1 Sol. A separate map check was labelled `anthropic/claude-sonnet-5-5`. The actual models were not independently verified, and the model for an earlier finite-torus check was not recorded. Agreement among models is not proof. No AI agent is an author.
