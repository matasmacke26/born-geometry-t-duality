# Generalised Born Geometry and T-Duality

Mathematica code and evaluated notebooks for generalised Born geometry, T-dual connections and curvature, and local regularity at non-extremal Killing horizons. These are companion computations for a PhD thesis.

The central construction retains the transported generalised Born structure when selecting the dual generalised Levi-Civita connection. This gives a pair of connections whose full generalised Riemann tensors are intertwined by T-duality on invariant sections. Selecting a connection independently from the dual metric and transported divergence need not give the same connection.

## Contents

The repository is intended to contain the following files. The instructions below use these filenames; remove download suffixes such as `(4)` when uploading the source and notebook.

| File | Purpose |
| --- | --- |
| `Born_TDual_GeneralisedGeometry.wl` | Readable Wolfram Language implementation. |
| `Born_TDual_GeneralisedGeometry.nb` | Self-contained notebook with an embedded implementation and saved evaluation results. |
| `README.md` | Scope, conventions, execution instructions and benchmarks. |
| `LICENSE` | MIT licence terms. |

The notebook embeds its implementation and does not load the standalone `.wl` file. Changes to one copy must therefore be reflected in the other before a release.

## What is computed

For each example, the implementation:

1. Constructs the original metric, zero initial two-form and anchor-induced Born structure.
2. Derives the dual metric and two-form from the transported generalised-metric eigenbundles.
3. Transports the tangent representative of the Born structure by
   $\hat\Psi=T_-\Psi T_+^{-1}$.
4. Constructs the original and dual Born connections from their Bismut-connection formulas.
5. Computes generalised torsion from its definition and forms
   $\mathcal D^{\mathcal B}_{\mathrm{LC}}=\mathcal D^{\mathcal B}-\frac13T_{\mathcal D^{\mathcal B}}$.
6. Evaluates the full generalised Riemann tensor, full Ricci curvature, mixed generalised Ricci tensor, generalised Ricci scalar and generalised Kretschmann scalar.
7. Transforms the dual connection into the transported regular frame, including the derivative term in the connection transformation law.
8. Recomputes curvature and torsion in that frame and checks transport identities.
9. Specialises the metric function and checks local horizon regularity.

The dual connection is selected by the **transported Born structure**, not by reconstructing the anchor-induced structure independently on the dual side.

## Requirements

- Mathematica with a Wolfram kernel for symbolic evaluation.
- The source targets Mathematica 13 or later. The supplied evaluated notebook records Mathematica 15.0 on Windows; compatibility with other versions should be checked by rerunning the validation cells.
- No external packages or datasets are loaded by the implementation.

Full symbolic evaluation can take several minutes or longer, depending on the system. The code prints progress as it constructs the geometry and evaluates curvature components. Individual horizon limits have a 30-second timeout; a timeout is reported as unresolved rather than accepted as a successful check.

## Run the evaluated notebook

1. Open `Born_TDual_GeneralisedGeometry.nb` in Mathematica.
2. Start a fresh kernel and evaluate the notebook from the beginning, including its embedded implementation cell.
3. Run the Schwarzschild and planar charged examples.
4. Confirm that both `AllChecksProved` and `AllRegularityChecksProved` return `True` for each example.

Saved outputs are a record of a previous evaluation. Re-evaluation is necessary to verify the current code in your environment.

## Run the standalone source

### Source synchronisation note

The reviewed standalone source contains an `AssociateTo` call near its end that passes two rules as separate arguments. The evaluated notebook already groups those rules correctly. Before running or releasing that version of the standalone source, replace the call with:

```wolfram
AssociateTo[checks, {
    "Original metric divergence" ->
        ZeroTensorQ[
            OriginalGeometry["LCDivergenceCoefficients"] -
            Join[metricDivOriginal, metricDivOriginal]
        ],
    "Dual Buscher divergence" ->
        ZeroTensorQ[
            TDualGeometry["LCDivergenceCoefficients"] -
            Join[metricDivDual - 2 dPhiDual,
                 metricDivDual - 2 dPhiDual]
        ]
}];
```

This README does not modify either implementation. Once the copies have been synchronised and reevaluated, update this note to record the corrected release.

### Execution

In a fresh kernel, load the source using its actual path:

```wolfram
Get["/absolute/path/to/Born_TDual_GeneralisedGeometry.wl"];

schwarzschild = BornTDual`RunBornExample["Schwarzschild"];
planarRN = BornTDual`RunBornExample["PlanarRN"];
```

Provided each calculation returns an association rather than `Failure[...]`, inspect the validation and horizon reports:

```wolfram
BornTDual`ShowValidation[schwarzschild]
BornTDual`ShowValidation[planarRN]

schwarzschild["AllChecksProved"]
planarRN["AllChecksProved"]

schwarzschildHorizon = BornTDual`HorizonReport[schwarzschild];
planarRNHorizon = BornTDual`HorizonReport[planarRN];

schwarzschildHorizon["AllRegularityChecksProved"]
planarRNHorizon["AllRegularityChecksProved"]
```

All four Boolean results should be `True`. `NOT PROVED`, `Inconclusive`, `Missing[...]` and `Failure[...]` must not be interpreted as successful verification.

Run examples sequentially. The driver uses shared package-private working variables; argument-free diagnostic helpers refer to the most recent run. Use the returned associations to retain and inspect each example separately.

## Inspect the results

Each result contains three generic geometries and their specialised versions:

| Key | Meaning |
| --- | --- |
| `"Original"` | Original geometry with generic $f(r)$. |
| `"DualPreferred"` | Dual geometry in the preferred coordinate eigenframe, away from $f=0$. |
| `"DualRegular"` | Dual geometry in the transported original eigenframe. |
| `"SpecialisedOriginal"` | Original geometry with the example's explicit $f(r)$. |
| `"SpecialisedDualPreferred"` | Specialised dual preferred-frame geometry. |
| `"SpecialisedDualRegular"` | Specialised dual regular-frame geometry. |

For example:

```wolfram
geom = schwarzschild["SpecialisedDualRegular"];

geom["RScalar"]
geom["KGen"]

BornTDual`ShowNonZeroConnection[geom, "LC"]
BornTDual`ShowNonZeroBornTorsion[geom]
BornTDual`ShowNonZeroRicci[geom]
BornTDual`ShowNonZeroGeneralisedRicci[geom]
BornTDual`ShowNonZeroIndependentRiemann[geom]
```

The component arrays include `"BornOmega"`, `"LCOmega"`, `"BornTorsionLow"`, `"RmArray"`, `"RmPair"`, `"RicFull"` and `"RcGen"`. Display routines retain entries that are not proved zero. The Riemann display removes pair-symmetry duplicates but does not reduce all algebraic Bianchi relations to a minimal independent set.

## Mathematical conventions

The coordinate order is $(v,r,\mathrm{th},\mathrm{ph})$, with $v$ the Eddington–Finkelstein Killing coordinate. In the planar example, the last two symbols represent planar coordinates.

The Courant pairing is

$$
\langle X+\xi,Y+\eta\rangle
=\frac12\bigl(\xi(Y)+\eta(X)\bigr).
$$

The ordered generalised frame is $(e_1^+,\ldots,e_4^+,e_1^-,\ldots,e_4^-)$, where $e_\mu^\pm=\partial_\mu+\iota_{\partial_\mu}b\pm g(\partial_\mu)$. In this frame, the Courant pairing and generalised metric are respectively $\text{diag}(g,-g)$ and $\text{diag}(g,g)$.

- `Omega[[A,B,C]]` denotes the coefficient of $e_C$ in $\mathcal D_{e_A}e_B$.
- The preferred flux is $H_{\mathrm{pref}}=H_{\mathrm{twist}}+db$.
- `RicFull` contracts the first and third curvature arguments with the inverse Courant pairing.
- `RcGen` retains the mixed eigenbundle blocks of `RicFull`.
- `RScalar` is one half of the inverse-generalised-metric trace of `RicFull`.
- `KGen` is the complete contraction of two curvature tensors with four inverse generalised metrics. In indefinite signature this is not a positive-definite norm.

The $28\times28$ curvature pair matrix uses increasing antisymmetric index pairs. Its Kretschmann contraction includes a factor of four to recover the unrestricted four-index sum. This bookkeeping factor is distinct from the separately checked zero-flux relation $K_{\mathrm{gen}}=4K_{\mathrm{ordinary}}$.

## Examples and benchmarks

Both models start with

$$
g=-f(r)\,dv^2+2\,dv\,dr+r^2g_\Sigma,
\qquad b=0,\qquad H_{\mathrm{twist}}=0,\qquad\Psi=I.
$$

| Example | Transverse metric | Metric function | Positive simple horizon |
| --- | --- | --- | --- |
| Schwarzschild | Unit two-sphere | $1-2M/r$ | $r_*=2M$, $M>0$ |
| Planar charged geometry | Flat two-plane | $-2M/r+Q^2/r^2$ | $r_*=Q^2/(2M)$, $M>0$, $Q\ne0$ |

The generalised Ricci scalar vanishes for both specialisations. The Schwarzschild mixed generalised Ricci tensor vanishes; the planar charged mixed tensor is non-zero for $Q\ne0$.

The generalised Kretschmann benchmarks are

$$
K_{\mathrm{Schwarzschild}}=\frac{192M^2}{r^6},
\qquad K_{\mathrm{Schwarzschild}}\big|_{r=2M}=\frac{3}{M^4},
$$

$$
K_{\mathrm{planar}}
=\frac{192M^2}{r^6}-\frac{384MQ^2}{r^7}+\frac{224Q^4}{r^8},
\qquad K_{\mathrm{planar}}\big|_{r=r_*}=\frac{20480M^8}{Q^{12}}.
$$

These benchmark expressions are not inserted into the curvature engine.

## Verification and scope

The supplied evaluated notebook records passing checks for connection compatibility, Born-involution compatibility, torsion, curvature symmetries, full curvature transport, divergence transport and horizon regularity in both examples.

The regular-frame calculation recomputes curvature from the transformed connection and anchor using the same curvature routine. It tests frame consistency; it is not a second implementation of the curvature definition. An additional ordinary Riemann contraction provides a separate check of the zero-flux factor-four identity.

The symbolic Buscher calculation assumes $f(r)\ne0$. Horizon regularity is then examined in the raw transported frame, with both one-sided limits, frame and metric non-degeneracy, and denominator checks for the rational component expressions. The angular chart is restricted to $0<\mathrm{th}<\pi$ in the spherical case.

The conclusions concern **local extension of the transported generalised geometry and its Born-selected connection**. They do not assert regularity of the ordinary dual metric or dilaton, global causal regularity, or geodesic completeness. At the horizon the extended dual generalised metric need not admit a semi-Riemannian metric parametrisation because its restricted anchors lose invertibility.

The driver is tailored to these invariant, locally untwisted examples. It does not implement arbitrary raw-flux T-duality, a general prescribed-divergence input, or a comparison with an independently selected dual metric-divergence connection. The examples have vanishing transported Born torsion and do not test a non-trivial torsion correction for arbitrary Born data.

## Citation

If these computations contribute to your work, please cite this repository and the associated thesis when its bibliographic details are available. Specify the release tag or commit hash used, together with your Mathematica version, so the calculation can be reproduced. Citation metadata can be supplied through `CITATION.cff` once the repository and thesis details are finalised.

## Licence

MIT License. See `LICENSE` for the full terms. The licence applies to this repository's code and accompanying documentation; Mathematica is a separate dependency governed by its own licence.
