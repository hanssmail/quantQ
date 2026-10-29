# quantQ unit tests

One test file per library module (`test_quantQ_<module>.q` ↔ `lib/quantQ_<module>.q`), a small
dependency-free harness (`quantQ_test.q`) and a runner (`runAll.q`). 42 modules, 854 assertions.

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

- **First run:** 841 tests, 115 failures. Six were mistakes in the tests and have been fixed.
- **Second run:** 843 tests, 109 failures. All 109 were genuine library defects.
- **Now:** every one of those defects is fixed in `lib/`, and 11 tests cover the new helpers. That
  makes 854 tests. The patched library has not been run yet, so expect a few follow-ups.

`../CHANGES.md` lists every fix by module, with the tests each one should turn green.
