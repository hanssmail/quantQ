/ Unit tests for lib/quantQ_covarianceShocks.q -- random-matrix comparison of covariance matrices (.quantQ.rmt)
/ Reference: Bouchaud et al., "Excess out-of-sample risk and fleeting modes" (arXiv:2205.01012).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`covarianceShocks;

/ rangeDTh -- edges of the eigenvalue spectrum of D (formula from the paper, evaluated in Python)
.qt.approx["rangeDTh: lamMax, q_in=0.2, q_out=0.25";{.quantQ.rmt.rangeDTh[10;50;10;40]`lamMax};4.1639235376052355;1e-10];
.qt.approx["rangeDTh: lamMin, q_in=0.2, q_out=0.25";{.quantQ.rmt.rangeDTh[10;50;10;40]`lamMin};0.2110764623947627;1e-10];
.qt.approx["rangeDTh: blur term (2.7e-3*N)^(-2/3)";{.quantQ.rmt.rangeDTh[10;50;10;40]`blur};11.111111111111109;1e-9];
.qt.approx["rangeDTh: no zero eigenvalues when T_out > N_out";{.quantQ.rmt.rangeDTh[10;50;10;40]`fracZero};0f;1e-12];
.qt.approx["rangeDTh: fraction of zero eigenvalues 1-T_out/N_out";{.quantQ.rmt.rangeDTh[20;100;20;10]`fracZero};0.5;1e-12];
.qt.approx["rangeDTh: lamMax, q_out=2";{.quantQ.rmt.rangeDTh[20;100;20;10]`lamMax};8.567627457812103;1e-9];
.qt.true["rangeDTh: lamMin < 1 < lamMax";{r:.quantQ.rmt.rangeDTh[10;200;10;200]; (r[`lamMin]<1) and 1<r`lamMax}];

/ covChange -- deterministic 3-feature data, in- and out-of-sample
rmtGen:{[n;s] flip `a`b`c!(sin s+0.7*til n; cos s+1.3*til n; (sin s+2.9*til n)+0.3*cos 0.4*til n)};
rmtIn:rmtGen[60;0.1];
rmtOut:rmtGen[60;5.3];
.qt.eq["covChange returns the documented keys";{key .quantQ.rmt.covChange[rmtIn;rmtOut]};`eigenvalues`lambdaMax`lambdaMaxNCor`overlap`inSampleVecs`fleetingVecs];
.qt.eq["covChange: one eigenvalue per feature";{count .quantQ.rmt.covChange[rmtIn;rmtOut]`eigenvalues};3];
.qt.approx["covChange: identical samples give D = I (all eigenvalues 1)";{.quantQ.rmt.covChange[rmtIn;rmtIn]`eigenvalues};1 1 1f;1e-6];
.qt.approx["covChange: the overlap array ends at 1 (fleeting mode is a unit vector)";{last .quantQ.rmt.covChange[rmtIn;rmtOut]`overlap};1f;1e-6];

/ covTabChange -- rolling version; one result row per input row
rmtTab:update t:i from rmtGen[40;0.5];
rmtTab:`t xcols rmtTab;
.qt.eq["covTabChange keeps one row per input row";{count .quantQ.rmt.covTabChange[`Tin`Tout!(15;15);rmtTab]};40];

.qt.finish[];
