/ Unit tests for lib/quantQ_stats.q -- descriptive statistics, histograms, tests, rank correlation (.quantQ.stats)
/ Reference values from scipy.stats.
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`stats;

/ exponential moving averages: s0 = x0, s_t = (1-l) s_(t-1) + l x_t
.qt.approx["expma";{.quantQ.stats.expma[0.5;1 2 3 4f]};1 1.5 2.25 3.125;1e-12];
.qt.approx["expma1 equals expma";{.quantQ.stats.expma1[0.3;3 1 4 1 5 9 2 6f]};.quantQ.stats.expma[0.3;3 1 4 1 5 9 2 6f];1e-12];
.qt.approx["expma of a constant is constant";{.quantQ.stats.expma[0.2;5#7f]};5#7f;1e-12];

/ histogram bin rules
.qt.eq["histGrid returns the bins unchanged";{.quantQ.stats.histGrid[0 1 2f;10?1f]};0 1 2f];
.qt.approx["histWidth";{.quantQ.stats.histWidth[0.5;0 1 2f]};0 0.5 1 1.5 2;1e-12];
.qt.approx["histBin";{.quantQ.stats.histBin[4;0 1 2f]};0 0.5 1 1.5 2;1e-12];
.qt.approx["histSquareRoot: ceil(sqrt n) bins";{.quantQ.stats.histSquareRoot "f"$til 9};0 8 16 24%3;1e-12];
.qt.approx["histSturge: 1+ceil(log2 n) bins";{.quantQ.stats.histSturge "f"$til 10};0 1.8 3.6 5.4 7.2 9;1e-12];
.qt.approx["histRice: ceil(2 n^(1/3)) bins";{.quantQ.stats.histRice "f"$til 20};(19%6)*til 7;1e-12];
.qt.approx["histScott: width ceil(3.5 sd / n^(1/3))";{.quantQ.stats.histScott "f"$til 8};0 5 10f;1e-12];
.qt.eq["getHistogram: counts add up to the sample size";{sum exec histogram from .quantQ.stats.getHistogram[0.1*til 50;.quantQ.stats.histBin[5]]};50];
.qt.eq["getHistogram: one row per bin edge";{count .quantQ.stats.getHistogram[0.1*til 50;.quantQ.stats.histBin[5]]};6];

/ overview statistics (nulls are excluded but counted)
statsOv:.qt.fixture["statsOv";{.quantQ.stats.overviewStats 1 2 3 4 10 0n}];
.qt.approx["overviewStats: mean";{statsOv[`values]0};4f;1e-12];
.qt.approx["overviewStats: standard deviation (population, as dev)";{statsOv[`values]1};3.1622776601683795;1e-12];
.qt.approx["overviewStats: adjusted skewness (scipy bias=False)";{statsOv[`values]2};1.6970562748477143;1e-12];
.qt.approx["overviewStats: excess kurtosis (scipy bias=True)";{statsOv[`values]3};-0.21199999999999974;1e-12];
.qt.approx["overviewStats: min, max, count, missing";{statsOv[`values]4 5 6 7};1 10 6 1;0f];

/ p-values from the normal table (coarse lookup)
.qt.approx["pValueTwoSided beyond the table";{.quantQ.stats.pValueTwoSided 3.5};0.001;1e-12];
.qt.approx["pValueTwoSided at 3.0";{.quantQ.stats.pValueTwoSided 3.0};0.01;1e-12];
.qt.approx["pValueTwoSided is symmetric";{.quantQ.stats.pValueTwoSided -3.5};0.001;1e-12];
.qt.approx["pValueRightSided";{.quantQ.stats.pValueRightSided 3.0};0.005;1e-12];
.qt.approx["pValueRightSided for a negative statistic";{.quantQ.stats.pValueRightSided -3.0};0.995;1e-12];
.qt.approx["pValueLeftSided";{.quantQ.stats.pValueLeftSided -3.0};0.005;1e-12];
.qt.approx["pValueLeftSided for a positive statistic";{.quantQ.stats.pValueLeftSided 3.0};0.995;1e-12];

/ t statistics
.qt.approx["tStatsCorr";{.quantQ.stats.tStatsCorr[0.5;27]};2.886751345948129;1e-12];
.qt.approx["tTestOneSample (scipy ttest_1samp)";{.quantQ.stats.tTestOneSample[5.1 4.9 5.6 5.8 6.0 5.5 5.3;5f]};3.1278490197651534;1e-10];
.qt.approx["tTestTwoSample, Welch (scipy ttest_ind equal_var=False)";{.quantQ.stats.tTestTwoSample[5.1 4.9 5.6 5.8 6.0 5.5 5.3;4.8 5.0 5.2 4.7 5.1 4.9]};3.0753513809518913;1e-10];

/ binomial
.qt.approx["factorial";{.quantQ.stats.factorial each 0 1 5 10};1 1 120 3628800f;0f];
.qt.approx["coeffBinomial";{.quantQ.stats.coeffBinomial[10;3]};120f;1e-9];
.qt.approx["pdfBinomial (scipy binom.pmf)";{.quantQ.stats.pdfBinomial[3;10;0.3]};0.2668279319999998;1e-12];
.qt.approx["pdfBinomial sums to 1";{sum .quantQ.stats.pdfBinomial[;12;0.37] each til 13};1f;1e-12];
statsCB:.qt.fixture["statsCB";{.quantQ.stats.characteristicsBinomial[0;10;0.3]}];
.qt.approx["characteristicsBinomial: mean, variance";{statsCB`mean`variance};3 2.1;1e-12];
.qt.approx["characteristicsBinomial: skewness (scipy)";{statsCB`skewness};0.27602622373694163;1e-12];
.qt.approx["characteristicsBinomial: excess kurtosis (scipy)";{statsCB`kurtosis};-0.12380952380952381;1e-12];
.qt.approx["binomialGen with an integer argument";{.quantQ.stats.binomialGen[5f;2]};10f;1e-12];
.qt.approx["binomialGen(1/2, 2) = -1/8";{.quantQ.stats.binomialGen[0.5;2]};-0.125;1e-12];
.qt.approx["binomialGen(x, 0) = 1";{.quantQ.stats.binomialGen[0.5;0]};1f;1e-12];

/ sign test and Wilcoxon signed-rank
.qt.approx["pValSignTest (scipy binom)";{.quantQ.stats.pValSignTest[2;10]};0.0546875 0.9892578125 0.109375;1e-12];
statsW:.qt.fixture["statsW";{.quantQ.stats.wilcoxonTest[11 22 33 40f;10 20 30 44f]}];
.qt.approx["wilcoxonTest: n excludes ties";{.quantQ.stats.wilcoxonTest[1 2 3f;1 0 0f]`n};2;0f];
.qt.approx["wilcoxonTest: W = sum of signed ranks (ranks start at 1)";{statsW`w};2f;1e-12];
.qt.approx["wilcoxonTest: sigma = sqrt(n(n+1)(2n+1)/6)";{statsW`sigma};sqrt 30f;1e-12];

/ rank correlation
.qt.approx["kendallTauRank (scipy kendalltau, no ties)";{.quantQ.stats.kendallTauRank[1 2 3 4 5;3 1 2 5 4]};0.4;1e-12];
.qt.approx["kendallTauRank of identical rankings is 1";{.quantQ.stats.kendallTauRank[1 2 3 4;10 20 30 40]};1f;1e-12];
.qt.approx["kendallTauRank of reversed rankings is -1";{.quantQ.stats.kendallTauRank[1 2 3 4;4 3 2 1]};-1f;1e-12];
.qt.eq["concordanceRoutine: concordant pair";{.quantQ.stats.concordanceRoutine[1 1;2 2]};100b];
.qt.eq["concordanceRoutine: discordant pair";{.quantQ.stats.concordanceRoutine[1 2;2 1]};010b];
.qt.eq["concordanceRoutine: tie";{.quantQ.stats.concordanceRoutine[1 1;1 2]};001b];
.qt.approx["somersD without ties in x equals Kendall's tau";{.quantQ.stats.somersD[3 1 2 5 4;1 2 3 4 5]};0.4;1e-12];
.qt.approx["somersDBinaryX = (1 + D)/2";{.quantQ.stats.somersDBinaryX[3 1 2 5 4;1 2 3 4 5]};0.7;1e-12];
.qt.approx["somersDBinaryXIntegerY";{.quantQ.stats.somersDBinaryXIntegerY[1 2 2 3 3 1;0 0 1 1 1 0]};0.8888888888888888;1e-12];
.qt.approx["kruskallGamma ignores ties";{.quantQ.stats.kruskallGamma[1 2 2 3 4;1 3 2 3 5]};1f;1e-12];
.qt.approx["kruskallGamma without ties equals Kendall's tau";{.quantQ.stats.kruskallGamma[1 2 3 4 5;3 1 2 5 4]};0.4;1e-12];
.qt.true["tKruskall has the sign of gamma";{(signum .quantQ.stats.tKruskall[1 2 3 4 5;3 1 2 5 4])=1}];

/ multiple testing
.qt.approx["pBonferroni";{.quantQ.stats.pBonferroni[0.05;5]};0.01;1e-12];
.qt.approx["pSidak";{.quantQ.stats.pSidak[0.05;5]};0.010206218313011495;1e-12];
.qt.eq["pHolm: simple case";{.quantQ.stats.pHolm[0.001 0.9;0.05]};enlist 1];
.qt.eq["pHolm: smallest p is compared with alpha/m";{.quantQ.stats.pHolm[0.0124 0.5 0.6 0.7;0.05]};enlist 1];
.qt.eq["pHolm is step-down: stop at the first non-rejection";{.quantQ.stats.pHolm[0.001 0.02 0.022 0.024;0.05]};enlist 1];

.qt.finish[];
