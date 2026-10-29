# quantQ unit tests

One test file per library module (`test_quantQ_<module>.q` ↔ `lib/quantQ_<module>.q`), a small
dependency-free harness (`quantQ_test.q`) and a runner (`runAll.q`). 42 modules, 843 assertions.

## Running

```
q tests/runAll.q                 # whole suite, from the repository root
q tests/test_quantQ_mat.q        # one module
cd tests && q runAll.q           # also works from inside tests/
```

Each run prints `FAIL` lines as they happen and ends with a per-module table and a list of failures. The
process exits with status 0 when everything passes and 1 otherwise, so it can be used in CI.

`test_quantQ_io.q` writes a scratch database to `tests/qtest_tmp_io` and removes it with `rm -rf`
(Linux/macOS).

## Writing tests

Every assertion takes a name and a lambda with no arguments, so an error inside the code under test is
caught and reported as a failure. The rest of the file still runs.

| Assertion | Passes when |
| --- | --- |
| `.qt.eq[name; {expr}; expected]` | `expected ~ expr` (exact, including type) |
| `.qt.approx[name; {expr}; expected; tol]` | every number is within `tol`; works on atoms, vectors, matrices, dictionaries and tables, and nulls must line up |
| `.qt.true[name; {expr}]` | `expr` returns `1b` |
| `.qt.throws[name; {expr}]` | `expr` signals an error |
| `x:.qt.fixture[name; {expr}]` | builds test data; an error is recorded as a failure instead of aborting the file |

Two q gotchas to watch for:

* A test lambda must not mention `y` or `z`, not even as a column name in a table literal, because that
  makes it a 2- or 3-argument function. Put such tables in a global first, for example `adaTab:([] y:...)`.
* `.qt.setModule` resets the random seed (`\S 42`), so the stochastic tests are repeatable.

## Where the expected values come from

Expected values were worked out independently of the library: numpy/scipy (linear algebra, distributions,
FFT, OLS, MLE fits), direct Python implementations of the textbook algorithms (DTW, Viterbi by brute
force, Needleman–Wunsch/Smith–Waterman, Levenshtein, Durbin–Levinson, the innovations algorithm, LOF),
closed-form identities (`i*j=k`, `H H = I`, `exp(i pi) = -1`), or invariants such as "probabilities sum
to 1" and "rank order is preserved". For stochastic routines (bootstraps, SGD, MCTS, random search) the
tests check properties rather than exact draws.

## Status

First run (kdb+, 29 Oct 2026): 841 tests, 115 failures. Six of those were mistakes in the tests and have
been fixed: the `dftTs` index column name, two `infth` tests that used `"a"` where q needs `enlist "a"`,
the Bollinger band column names, the starting counter of the `mcts` bandit, and an unequal-length
Needleman–Wunsch case. Every other failure is a confirmed library defect, listed below.
`covTabChange` also prints a stray number (`![0N;iP]` is a leftover debug print).

## Known library defects (each has a failing test)

**Errors: undefined names, wrong namespaces, or locals that q lambdas cannot see**
- `hopf.getHopfieldNet2S` uses `state` (should be `state1`), which breaks both `runHebbDescent*`. A
  Python simulation suggests that even after the rename the descent cannot learn the mapping, because
  the outer product is the wrong way round (`input *\: target`).
- `hopf.getGroundState`: the inner lambda uses `dimH`, a local of the outer function, and calls
  `.quantQ.ephys.getMinEnergyHopfield`.
- `tda.getAllLoops`, `getAllPureLoops*` and `nLoopsThr` read a global `arrayEdges` from inside a nested
  lambda. The TDA notebook only works because it defines that global. `tda.bornAndDeathLoops` reads a
  global `tabEdges` (its parameter is `tE`).
- `knn.errorLocal` uses `weight` (should be `weights`).
- `mcts.evaluateMCTS` reads a global `bucket`.
- `meb.blockShuffle` uses `xT` (should be `data`).
- `nn.predictNNDeep` never initialises `nodesT1`.
- `ols.olsTab` uses `count[y]`. `ols.byR2` refers to `y`/`yPred` after aggregating them away.
- `svm.oneStep` reads a global `data` (should be `bucket[`data]`), which breaks `fitSoftMargin` and
  `cvFit`.
- `rmt.covChange` builds its result with `eigenvalues` instead of `` `eigenvalues ``.
- `math.erlangCDF` uses an undefined `K`.
- `quantum`: every parametrised gate defaults to `.quantQ.complex.pi`, which does not exist
  (`quantQ_complex.q` defines `.quantQ.pi`).
- `ts.HannanRissanen` calls `.quantQ.ts.residualsAR`, which is not defined. `ts.estimate` calls
  `quantQ.ts.innovations` (missing the leading dot). `ts.forecast` calls `forecastMA` with 3 arguments
  (it takes 2). The VARX section redefines `zt`, `Z`, `Y` and `covarianceResidual` with different
  arities, which breaks `varpest`, `estimate1` and the other VAR functions. `ts.varXest` is in the
  wrong namespace, and `ts.confInterval` uses an undefined `d` (neither is tested).

**Wrong results from right-to-left evaluation**
- `quat.mult`: `r-i-j-k` is `r-(i-(j-k))`, so the real, j and k parts are wrong (e.g. `j*j` returns +1).
- `ols.fit`, `ols.rSquaredAdj`, `ols.AICC`: `n-1-p` is `n-(1-p)`, so the adjusted R² and AICC use the
  wrong degrees of freedom.

**Other wrong results**
- `ada.adaPredictionStats`: false positives and false negatives are swapped, which makes recall,
  precision and specificity wrong.
- `svm.evalStats`: `falseDiscoveryRate` is TP/(TP+FP) (that is precision) and `falseOmissionRate` is
  TN/(TN+FN) (that is NPV).
- `hmm.viterbi`: the state path is the per-step argmax of delta, with no backtracking.
- `meb.BSorig`/`meb.BS` reorder with `asc[xTS]@iasc xT` instead of by rank. The block bootstraps also
  produce nulls: block starts run one past the end of the series.
- `math`: `gammaApprox` returns Γ(x+1) (missing `1/x`), which also breaks `gammaPDF`/`gammaCDF`.
  `genLaguerrePol` is off by one in the recurrence. `first/secondChebyshevPol` return nothing for n>1.
  `getChebyshevSecondRoots` uses k=0..n-1 instead of 1..n. `normalCDF` passes sigma where the variance
  is expected. `tetrachoricExpansionBivar` includes the k=0 term and uses He_k instead of He_(k-1).
  The Laguerre series in `gammaIncompleteUpperApprox` does not converge for z>1/2.
- `ols.RMSE` uses `t wavg t`, which gives inf or null instead of the root mean square.
- `pois.logLPoissonNegReg` subtracts the L2 penalty. `intensity` stores θ'x rather than exp(θ'x), so
  `yVarPRED` is wrong.
- `stats`: `pValueTwoSided` is wrong for negative statistics. The t-tests use the population sd (`dev`).
  `wilcoxonTest` ranks from 0. `pHolm` is off by one and not step-down. `somersDBinaryXIntegerY`
  double-counts repeated (x,y) pairs.
- `stats.NeedlemanWunsch` only works for equal lengths, because the traceback starts at
  `(Nseq1;Nseq2)` on a grid indexed `[seq2;seq1]`. `stats.LevenshteinDistance` returns `()` when one
  sequence is empty.
- `ta.msd` takes `sqrt` of `mdev`, so the Bollinger bands are wrong too.
- `infth.proteinGetSubstring` ignores `N` and always takes 3 characters.
- `ts.simMA` is wrong for q≥2. `ts.forecastAR` multiplies phi1 by the oldest lag.
- `upd.trades`: `-3#` wraps around when a symbol has fewer than 3 ticks, so the first VWAPs are wrong
  (use `-3 sublist`).
- `mat.det` returns a 1x1 matrix, not a number, for 1x1 input.

**Not tested:** `.quant.stats.runs` in `quantQ_stats.q` is in the namespace `.quant`, not `.quantQ`.
