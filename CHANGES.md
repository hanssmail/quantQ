# Library fixes found by the unit tests

Every change below fixes a defect found by `tests/` (854 tests). The patched library has **not been run
yet**. Run `q tests/runAll.q` to confirm. The full diff is in `quantQ-lib-fixes.diff`.

## Crashes: undefined names, wrong namespaces, closures

| Module | Fix |
| --- | --- |
| `complex` | Defines `.quantQ.complex.pi`, which the quantum gates use (`.quantQ.pi` is kept). |
| `covarianceShocks` | `covChange` keys use `` `eigenvalues `` (the backtick was missing). The in-sample PCA vectors are now the eigenvectors of the in-sample covariance; previously `pca` was applied to the covariance matrix as if it were data. The debug print `![0N;iP]` in `covTabChange` is removed. |
| `hopf` | `getHopfieldNet2S` uses `state1`, not the undefined `state`. `getGroundState` passes `minE` into its inner lambda, uses `.quantQ.hopf.getMinEnergyHopfield` instead of `.quantQ.ephys...`, starts from random states, and stops after 10000 tries. |
| `knn` | `errorLocal` uses `weights`, not `weight`. |
| `mcts` | `evaluateMCTS` no longer reads a global `bucket`. It accepts the tree (evaluated for player 1) or the MCTS bucket (the player to move comes from `movesState`). |
| `meb` | `blockShuffle` uses `data`, not `xT`. |
| `nn_deep` | `predictNNDeep` initialises `nodesT1`. |
| `ols` | `olsTab` uses the fitted `nObservations`, not `count[y]`. `byR2` computes the overall R² from the input table. |
| `svm` | `oneStep` uses `bucket[`data]`, not a global `data`. |
| `tda` | `getAllLoops` and `getAllPureLoops` pass `bucket[`arrayEdges]` into the nested lambda; they used to rely on a global `arrayEdges`. `bornAndDeathLoops` uses its parameter `tE`, not the global `tabEdges`. |
| `math` | `erlangCDF` uses its scale parameter (`K` was undefined) and returns its value (the `:` was missing). |
| `ts` | `estimate` calls `.quantQ.ts.innovations` (the leading dot was missing). New `residualsAR`, which `HannanRissanen` needs. `confInterval` uses `y`/`p`, not the undefined `d`. |
| `stats` | `.quant.stats.runs` is now `.quantQ.stats.runs`, with the old name kept as an alias. |

## `ts`: VAR and VARX name clash

The VARX section redefined `zt`, `Z`, `Y` and `covarianceResidual` with different arities, which broke
the VAR functions defined above it. The VARX versions are renamed `ztX`, `ZXlag` and
`covarianceResidualX`; the VARX `Y` duplicate was buggy and is dropped in favour of the VAR one.
`varXest` moves into `.quantQ.ts`, keeping `.quantQ.varXest` as an alias, and now computes
`Y Z' (Z Z')^-1` (the brackets were wrong). `estimateX` passes its locals into its lambdas.
`estimateNX` forecasts the exogenous series with `estimate1` and recurses; previously it returned a
lambda without calling it. One open question: `estimateX` still uses `x[;0]` for the contemporaneous
exogenous term. I left that as it was.

## Right-to-left evaluation

| Module | Fix |
| --- | --- |
| `quat` | `mult`: the Hamilton product is written out with every term in brackets. |
| `ols` | `rSquaredAdj`, `AICC` and `fit` (adjusted R², AICC) used `n-1-p`, which q reads as `n-(1-p)`. They now use `n-(p+1)`. |

## Wrong results

| Module | Fix |
| --- | --- |
| `ada` | `adaPredictionStats`: false positives and false negatives were swapped. |
| `svm` | `evalStats`: `falseDiscoveryRate` is now FP/(FP+TP) and `falseOmissionRate` is FN/(FN+TN). |
| `hmm` | `viterbi` backtracks through stored back-pointers. Previously it took the argmax of delta at each step. |
| `meb` | `BSorig`/`BS` place the sorted draws by rank (`asc[xTS]@iasc iasc xT`). Block indices are now `s..s+block-1`; they used to run one past the end and produce nulls. |
| `math` | `gammaApprox` gets its `1/x` factor. New `.quantQ.math.gamma` (Lanczos), used by `gammaPDF` and `gammaCDF`. `genLaguerrePol` has the correct recurrence and is now iterative rather than exponentially recursive. `gammaIncompleteLowerApprox` uses the convergent power series, and `gammaIncompleteUpperApprox` = Γ − lower; the Laguerre series did not converge for z>1/2. The Chebyshev polynomials return their value for n>1. `getChebyshevSecondRoots` uses k=1..n. `normalCDF` integrates the correct density with the midpoint rule, vectorised. `tetrachoricExpansionBivar` sums k=1..n with He_(k-1). |
| `ols` | `RMSE` is `sqrt avg e*e`. |
| `pois` | The L2 penalty is added to the minimised objective. `intensity` is exp(θ'x), so `yVarPRED` is the mode of Poisson(λ) and the log-likelihood is Σ y log λ − λ. |
| `quantum` | `U1` = diag(1, e^{iλ}). `U2` includes the 1/√2 factor. `Rz` uses `theta`. |
| `stats` | `pValueTwoSided` uses `abs x`. The t-tests use `sdev`. Wilcoxon ranks start at 1. `pHolm` uses α/(m+1−k) and is step-down. New `pairStats` enumerates each pair once (repeated observations were double-counted); it is used by Kendall, Somers and Kruskal. |
| `stats` (bioInf) | `NeedlemanWunsch` starts the traceback at `(Nseq2;Nseq1)`, so sequences of different lengths work. `LevenshteinDistance` returns `[Nseq2;Nseq1]`, which also fixes empty sequences. |
| `ta` | `msd` is `mdev`; the extra `sqrt` was wrong. This also fixes the Bollinger bands. |
| `infth` | `proteinGetSubstring` and `proteinInnerProduct` use `N` instead of a hard-coded 3. |
| `ts` | `simMA` is correct for any q. `forecastAR` multiplies phi1 by the latest value. `forecastMA[data;period;q]` now matches how `forecast` calls it: it gives h-step forecasts from the innovations estimates, and the old one-step predictor is kept as `forecastMAOneStep`. `vec` is `raze flip a`. |
| `upd` | `-3 sublist` instead of `-3#`, which wrapped around when fewer than 3 ticks were cached. |
| `mat` | `det` of a 1x1 matrix returns a number. |
| `hopf` | Hebbian descent builds W = target ⊗ input, so that `W mmu input` gives the target. The old orientation could not learn the mapping. |
