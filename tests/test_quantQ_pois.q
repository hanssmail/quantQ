/ Unit tests for lib/quantQ_pois.q -- Poisson distribution and Poisson regression (.quantQ.pois)
/ Maximum-likelihood reference values from scipy.optimize (BFGS) on the same data.
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`pois;

/ distribution
.qt.approx["poissonPDF(3 | 0.3)";{.quantQ.pois.poissonPDF[0.3;3]};0.00333368199306773;1e-14];
.qt.approx["poissonPDF(0 | lambda) = exp(-lambda)";{.quantQ.pois.poissonPDF[1.7;0]};exp -1.7;1e-14];
.qt.eq["poissonPDFTab has maxN rows starting at 0";{exec n from .quantQ.pois.poissonPDFTab[2f;enlist[`maxN]!enlist 5]};til 5];
.qt.approx["poissonPDFTab probabilities sum to ~1";{exec sum prob from .quantQ.pois.poissonPDFTab[2f;()!()]};1f;1e-10];
.qt.approx["poissonCDFTab ends at ~1";{last exec prob from .quantQ.pois.poissonCDFTab[2f;()!()]};1f;1e-10];
.qt.true["poissonCDFTab is non-decreasing";{all 0<=1_deltas exec prob from .quantQ.pois.poissonCDFTab[2f;()!()]}];

/ log-likelihood (terms not depending on theta are dropped)
poisSmall:([] yVar:1 2; xVar:(1 0f;1 1f));
.qt.approx["logLPoissonNeg";{.quantQ.pois.logLPoissonNeg[poisSmall;0.1 0.2]};1.7550297256516507;1e-12];
.qt.approx["logLPoissonNegReg with lambda=0 equals logLPoissonNeg";{.quantQ.pois.logLPoissonNegReg[poisSmall;0f;enlist[`addConstant]!enlist 1b;0.1 0.2]};1.7550297256516507;1e-12];
/ the penalised objective that is minimised must grow with the penalty: negLL + n*lambda*|theta_1|^2
.qt.approx["logLPoissonNegReg adds the L2 penalty (constant not penalised)";{.quantQ.pois.logLPoissonNegReg[poisSmall;0.5;enlist[`addConstant]!enlist 1b;0.1 0.2]};1.7550297256516507+2*0.5*0.04;1e-12];

/ regression: y counts, xVar includes the constant
poisX:0 0.1 0.25 0.4 0.5 0.6 0.75 0.9 1 1.2 1.4 1.5 1.7 1.9 2f;
poisY:1 2 1 2 3 2 4 3 5 4 6 7 8 10 9;
poisTab:([] yVar:poisY; xVar:1f,'poisX);
poisInit:enlist[`initPoint]!enlist (0 0f;0.5 0f;0 0.5);
poisModel:{.quantQ.pois.regPoisson[poisTab;poisInit,enlist[`precision]!enlist 1e-10]};
.qt.approx["regPoisson: maximum-likelihood coefficients";{poisModel[]`theta};0.36722149336504645 0.9896223958667972;1e-3];
.qt.approx["regPoisson: log-likelihood at the optimum";{poisModel[]`logLikelihood};45.63075252571209;1e-5];
.qt.approx["regPoisson: intensity is exp(theta'x)";{exec intensity from poisModel[]`output};exp 0.36722149336504645+0.9896223958667972*poisX;1e-2];
.qt.eq["predPoisson: output columns";{asc cols .quantQ.pois.predPoisson[poisTab;poisModel[]]};asc `yVar`xVar`intensity`yVarLIST`logL`yVarPRED`corr];
.qt.approx["predPoisson: predicted count is the mode of Poisson(exp(theta'x))";{exec yVarPRED from .quantQ.pois.predPoisson[poisTab;poisModel[]]};floor exp 0.36722149336504645+0.9896223958667972*poisX;0f];

/ regularised regression (constant added by the function)
poisTabNC:([] yVar:poisY; xVar:enlist each poisX);
poisModelReg:{[lam] .quantQ.pois.regPoissonReg[poisTabNC;lam;poisInit,enlist[`precision]!enlist 1e-10]};
.qt.approx["regPoissonReg with lambda=0 equals the MLE";{poisModelReg[0f]`theta};0.36722149336504645 0.9896223958667972;1e-3];
.qt.approx["regPoissonReg with lambda=0.05 (scipy, penalising the slope only)";{poisModelReg[0.05]`theta};0.44556447851770353 0.9301470977744946;1e-3];
.qt.true["regPoissonReg: a penalty shrinks the slope";{(abs last poisModelReg[0.05]`theta)<abs last poisModelReg[0f]`theta}];
.qt.eq["regPoissonReg records lambda and addConstant";{poisModelReg[0.05]`lambda`addConstant};(0.05;1b)];
.qt.eq["predPoissonReg returns one row per observation";{count .quantQ.pois.predPoissonReg[poisTabNC;poisModelReg[0f]]};15];
.qt.eq["regPoissonCV returns one value per fold";{count .quantQ.pois.regPoissonCV[poisTabNC;poisInit;3;0f]};3];

.qt.finish[];
