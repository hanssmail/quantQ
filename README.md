# quantQ

The repository for the Machine Learning and Big Data with kdb+/q book by Novotny et al.

Order the book at https://www.wiley.com/en-us/Machine+Learning+and+Big+Data+with+KDB%2B+Q-p-9781119404750.

## Getting Started

```
$ cd quantQ/lib
$ q quantQ.q -p 5000
```

The naming convention for each .q file reflects the corresponding book chapter and the functions it defines reside under a homonymous namespace, as outlined in Section 3.1.0.1.

## Errata

### Chapter 7: Joins

| Section       | Note                
| ------------- |--------------------
| 7.1           | The example of comma join corresponds to ```t1,t6``` and not ```t5,t6```.


| Section       | Note                
| ------------- |--------------------
| 7.2.10        | In the text we state that *we aim to aggregate the data from the table dataSet2 over a window starting 1 minute prior to the trade and ending at the time of the trade.*; however, the ```window``` is defined as starting at the time of the trade and ending 1 minute after the trade. The example should read ```window: (-00:01:00;0) +\: exec time from dataSet1;```.

### Chapter 14: Time Series Econometrics

| Section       | Note                
| ------------- |--------------------
| 14.1.6.1      | Ordering of ```phi``` vector inside implementation of  ```.quantQ.ts.simAR``` should be reversed to be in line with definition 14.1 and example on page 276. Addressed in the repo. Variation of the function using adverbs also provided, under ```.quantQ.tse.simAR```


### Chapter 15: Fourier Transform

| Section       | Note                
| ------------- |--------------------
| 15.3          | Example following implementation of the Hamilton product (15.29-15.32) should read as ```.quantQ.quat.mult[quat1;quat3]```. Definition of ```.quantQ.quat.mult``` has a typo which is fixed in the repo.    


### Chapter 22: Neural Networks

| Section       | Note                
| ------------- |--------------------
| 22.2          | Missing functions ```.quantQ.nn.funcNN``` and ```.quantQ.nn.funcErrNN``` have been added to the repo.

## Extension beyond the book

### Mathematical functions

We have added ```.quantQ.math``` namespace with various mathematical functions, constants and identities. Currently, there are constants, hyperbolic functions, number of special functions and polynomials (defined in the real domain), and the most frequently used PDF and CDF. In order to obtain the multivariate normal distribution, we have included tetrachoric expansion.

### Biostatistics

We have added into the ```.quantQ.stats``` namespace functions to work with contingency tables, namely Exact Fisher test and Barnard test.

### Optimization

We have added a native implementation of Nelder-Mead (amoeba) optimisation method. Functions can be found in the ```.quantQ.amoeba``` namespace. The illustration shows a solution to Rosenbrock's function.  

### Bioinformatics

We have implemented the Needleman-Wunsch algorithm developed in bioinformatics which can be used to align two sequences (nucleotide sequences or general finite-length sequences) using principles of dynamic programming. Functions are within the ```.quantQ.stats``` namespace. The local matching using the Smith-Waterman algorithm is available as well. Levenshtein distance calculated with Wagner-Fischer algorithm has been added.

### Dynamic Time Warp

We have added the algorithms to perform the Dynamic Time Warp calculations in the ```.quantQ.dtw``` namespace. Details of the algorithm can be found [here](data/dtw_intro.pdf).

### Deep Neural Networks

We have added deep neural networks. The features implemented include the dropout, batches, and different functional forms of annealing in learning and regularisation. The architecture is specified in the input dictionary; the rest is fully automated. More examples can be found [here](quantQ_nn_deep.md). The implementation is part of the ```.quantQ.nn``` namespace.

### Stochastic Optimisation

We have implemented the stochastic optimisation, which can be used to minimise the provided n-dimensional function using random search. The provided method is an iterative procedure which explores at every step a neighbourhood comprising of random points lying on the n-sphere, where the radius of the sphere is shrinking with occasional attempts to investigate further points. The implementation is part of the ```.quantQ.so```  namespace.

### Support Vector Machine

We have added the Support Vector Machine for binary classification using the Soft Margin and Stochastic Gradient Descent (using one observation per step). Functions are specified with a set of the default setup, which is customised based on the provided dataset. The Soft Margin specification is provided with regularisation parameter, which can be optimised using the built-in n-fold cross-validation method.

The implementation includes a function which calculates more than 20 statistics used to evaluate the binary classification, including specificity, accuracy, precision, or F1. More examples can be found [here](quantQ_svm.md). The implementation is part of the ```.quantQ.svm``` namespace.

### Poisson Regression

We have added the library with the Poisson distribution and the Poisson regression to estimate the integer-valued Poisson process. The functions are based on the maximum likelihood optimised using routines from ```.quantQ.opt``` library. The library also contains the ```L2```-regularised version with ```n```-fold cross-validation.

Examples and basic usage of the library can be found in [here](quantQ_pois.md). The implementation is part of the ```.quantQ.pois``` namespace.


### Quantum Computing

We have added name space ```.quantQ.quantum``` which contains a set of routines to set up and perform quantum computing using quibits. The library is not connected to any actual quantum computer and is for demonstration purposes only.

Examples and basic library usage can be found in [here](quantQ_quantum.md).

### Topological Data Analysis

We have added the library ```.quantQ.tda``` which allows us to perform the topological data analysis for a cloud on ```n```-dimensional points. The routines provided include calculation of the distance between points, identification of the all Vietoris-Rips complexes given the provided threshold, the routine to calculate all unique loops, and several analytical functions to get insight into data using the TDA.

Examples and basic library usage can be found in [here](quantQ_tda.md).

### Utility Library

We have added the library ```.quantQ.util``` with various utility functions.

### Hopfield Neural Networks

We have added the library ```.quantQ.hopf``` to define and work with Hopfield networks. In addition, we have added Hebbian descent to calculate train the network for a provided input/output.

### Technical Analysis

We have added the library ```.quantQ.ta``` which comprises a growing list of technical indicators from what is known as Technical Analysis. The details of the library can be found in [here](quantQ_ta.md).

### RSA Cryptoanalysis

We have added the library ```.quantQ.rsa```. The library contains the implementation of the RSA algorithm to encrypt message using public and private key. The purpose of the library is to illustrate the RSA concept. The library contains several additional functions to work with prime numbers and perform factorisation.  

### Random Matrix Theory

We have added the library ```.quantQ.rmt```. The library contains implementation of Excess Out-of-Sample Risk and Fleeting Modes paper by Bouchaud et al. [here](https://arxiv.org/abs/2205.01012). The algorithm performs the comparison of two datasets in terms of their empirical covariance matrices using the Random Matrix Theory. In order to implement the algorithm, we have added the Denman-Bevaers and power series algorithm to calculate the square-root of the matrix along with several utility functions. Further, ```.quantQ.stats``` has been extended by generalised binomial number calculator.

### Maximum Entropy Block Bootstrap

We have added he library ```.quantQ.meb```. The library implements the block bootstrap using the Maximum Entropy as proposed in Bergamelli, Novotny and Urga, 2015 [here](https://econpapers.repec.org/article/risactuec/0115.htm). The method is intended to be used to create a bootstrapped samples for non-stationary time series, where the temporal component is preserved by the blocks with the overall time series being replicated by using the Maximum Entropy bootstrap. For completeness, we provide an old implementation by Vinod and Lopez for reference. 

## Changes in v2.0

In v2.0 more 854 unit tests were added. Every change below fixes a defect found by `tests/` (854 tests). Run `q tests/runAll.q` to confirm.

### Crashes: undefined names, wrong namespaces, closures

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

### `ts`: VAR and VARX name clash

The VARX section redefined `zt`, `Z`, `Y` and `covarianceResidual` with different arities, which broke
the VAR functions defined above it. The VARX versions are renamed `ztX`, `ZXlag` and
`covarianceResidualX`; the VARX `Y` duplicate was buggy and is dropped in favour of the VAR one.
`varXest` moves into `.quantQ.ts`, keeping `.quantQ.varXest` as an alias, and now computes
`Y Z' (Z Z')^-1` (the brackets were wrong). `estimateX` passes its locals into its lambdas.
`estimateNX` forecasts the exogenous series with `estimate1` and recurses; previously it returned a
lambda without calling it. One open question: `estimateX` still uses `x[;0]` for the contemporaneous
exogenous term. I left that as it was.

### Right-to-left evaluation

| Module | Fix |
| --- | --- |
| `quat` | `mult`: the Hamilton product is written out with every term in brackets. |
| `ols` | `rSquaredAdj`, `AICC` and `fit` (adjusted R², AICC) used `n-1-p`, which q reads as `n-(1-p)`. They now use `n-(p+1)`. |

### Wrong results

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
