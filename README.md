# LLL lattice spectral stability

Reproducibility supplement to Davis Rens, *Certified non-hexagonal spectral stability in the lowest-Landau-level equation*, SuperMatics Technical Report SM-TR-2026-02 (2026), report DOI [10.5281/zenodo.24326107](https://doi.org/10.5281/zenodo.24326107). The report is not peer reviewed.
Repository: https://github.com/SuperMatics-Inc/lll-lattice-stability.
Version 1.0.0; software supplement DOI: [10.5281/zenodo.24326109](https://doi.org/10.5281/zenodo.24326109). Both DOIs are reserved; the deposits are not yet published.

## Statements and limits

| Label | Statement | Lean scope |
|---|---|---|
| PROVED | Gaussian Bloch symbol, covariance and cancellation below degree four | Algebraic consequences of supplied moment identities; analytic covariance is assumed |
| CERTIFIED | Positive determinant throughout the punctured Bloch torus for tau=-2/5+19i/20 and tau=2/5+i | Conditional positivity, actual 776/719-leaf covers and near-origin scaling; analytic enclosures remain hypotheses |
| CERTIFIED | Hexagonal dispersion has two simple inflection points, derivative sum at least 1.20 and theta-derivative expression greater than 26.95364 | Theta evaluation and interval arithmetic are not formalized |
| CERTIFIED PARTIAL | 4178 stable closed shape boxes, 25088 unstable boxes, cusp instability for Im(tau)>=3 and a long-wave transition curve | Separate GSTMap checks the 40570-leaf partition and conditional reduction; not the Gaussian symbol, interval arithmetic or Gamma4 |
| NUMERICAL | Sampled landscapes, S sampling and finite-torus spectra | No continuum certificate |

The origin has a nonzero nilpotent block. Real spectrum is not uniform boundedness, weighted/dispersive stability or nonlinear stability. These examples contradict a spectral reading of GST's expectation, not their theorems. The partial map is not a classification or identification of the full spectral boundary. Its X leaves remain undecided; their coordinate area is 5619/134217728.

S has **one interval implementation**, with independent floating-point sampling of all 4178 leaves. A second-arithmetic attempt covered only **12/4178** leaves; 4166 remain uncovered by it. Gamma4 has **one interval backend**, shared with the independent U verifier. U has Arb and cancellation-preserving mpmath.iv provenance; the public map command runs only the latter. The analytic formulas and explicit tails remain part of the trust base. See `gstmap/PROOFS.md` for the mathematical arguments.

## Reproduction

Install Git, Python 3.12+, [uv](https://docs.astral.sh/uv/) and, for formal checks, [elan](https://github.com/leanprover/elan).

```sh
python scripts/verify_pins.py
python scripts/verify.py --output replay-output
uv run --no-project --with-requirements requirements.txt python gstmap/replay.py --evidence --workers 2 --output map-output
uv run --no-project --with-requirements requirements.txt python figures/plot_partial_map.py --output-dir figure-output
```

The original runner executes the Arb and floating-point suite, the independent mpmath.iv two-point/hexagonal certificates, exact rational cover checks and the physical-space M=256 torus (65536 characters). The partial-map command checks exact partition data, every U leaf, all 397 curve patches and 396 joins, negative controls, independent S sampling and the cusp. Every child exit is recorded in output JSON. `--evidence` omits the expensive primary Arb S cover; omit that flag to execute that cover, not a second S implementation. The generated Lean partition must match the shipped LF bytes exactly.

The `partial-map` CI job deliberately runs a bounded subset:
`uv run --no-project --with-requirements requirements.txt python gstmap/replay.py --ci-subset --workers 4 --output map-output`.
It checks the exact partition and LF generated-data match, all 397 Gamma4 patches and 396 joins, U and Gamma4 negative controls, and the cusp. It **does not** replay the 25088-leaf U cover, S sampling or the primary Arb S cover; these omissions are listed in its receipt. The full evidence command above still checks every U leaf and all S samples. Neither mode supplies a second S implementation or a second Gamma4 backend.

`figures/plot_partial_map.py` renders the exact partition, not a new spectral calculation. Thin X boxes may be smaller than a print pixel. White below the unit arc is outside the fundamental domain.

For each formal project:

```sh
cd lean                 # separately, gstmap/lean
lake exe cache get
lake build
```

Both use Lean 4.35.0 and Mathlib v4.35.0, resolved to `56a005855ea17e0b0958f73ec801ae140b0a1679` in both manifests. No custom axioms or admitted proofs are used. Original 37 and map eight axiom checks use only propositional extensionality, classical choice and quotient soundness. CI builds both public git-pinned projects; focused Mathlib imports keep the map project within its unchanged 4096 MiB per-process limit. It frees unused runner SDKs before fetching Lean dependencies. The theorem statements and proof bodies are unchanged.

The SHA256SUMS manifest is computed from an LF Mac checkout and verified on a fresh clone. Binary PDFs are excluded from newline translation. Raw execution logs and the manuscript are maintained privately, not distributed here. `requirements.txt` pins the exercised numerical and plotting package versions; the runners record the resolved Python/package versions.

## Attribution and trust

D.R. conceived and directed the programme, selected and framed the problems, orchestrated the research and the verification, wrote and edited the manuscript extensively, and reviewed the full manuscript and its arguments in depth. He takes responsibility for the final work. LLM-based agents generated the proofs, computer-assisted certificates, implementations, Lean sources and manuscript drafts and edits, and performed independent validation. Recorded runtime configuration labels were Anthropic Claude Opus 5.5, Anthropic Claude Sonnet 5.5 and OpenAI GPT-6.1 Sol; separate partial-map validation used runtime configuration label anthropic/claude-sonnet-5-5. Served-model identity was not independently verified. An earlier finite-torus validation agent's identity was not recorded. Independence is from the generating agent, not a non-LLM or human re-derivation. This preserves the model-attribution trust caveat: configuration labels are not proof of served identity. Model agreement is not a proof. No LLM is an author.

The supplement is Apache-2.0 licensed. No experimental data or third-party manuscript is redistributed.
