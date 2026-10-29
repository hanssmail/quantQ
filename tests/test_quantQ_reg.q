/ Unit tests for lib/quantQ_reg.q -- ridge regression and least-angle regression (.quantQ.reg)
/ Ridge reference values: numpy.linalg.solve(X'X + lambda I, X'y) (the constant is penalised too).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`reg;

regXv:1 2 3 4 5 6 7 8 9 10f;
regY:3.3 4.8 7.1 8.6 11.25 12.95 15.15 16.7 19.2 21.0;
regX:(10#1f;regXv);
regYexact:1+2*regXv;

/ ridge
.qt.approx["ridge with reg=0 is OLS";{.quantQ.reg.ridge[regY;regX;0f]`beta};1.0233333333333308 1.9966666666666677;1e-10];
.qt.approx["ridge with reg=2.5";{.quantQ.reg.ridge[regY;regX;2.5]`beta};0.629209621993123 2.0397250859106535;1e-10];
.qt.true["ridge shrinks the coefficient vector";{(sum t*t:.quantQ.reg.ridge[regY;regX;10f]`beta)<sum u*u:.quantQ.reg.ridge[regY;regX;0f]`beta}];
.qt.eq["ridge reports the sample and the parameter";{.quantQ.reg.ridge[regY;regX;2.5][`nObservations`dfn`regParameter]};(10;2;2.5)];
.qt.eq["ridge drops missing observations";{.quantQ.reg.ridge[@[regY;0;:;0n];regX;1f]`nObservations};9];
.qt.approx["ridgeInOut: exact data, reg=0, zero out-of-sample error";{.quantQ.reg.ridgeInOut[regYexact;regX;0f;0.7]`RMSE};0f;1e-9];
.qt.eq["ridgeInOut: out-of-sample size";{.quantQ.reg.ridgeInOut[regYexact;regX;0f;0.7]`nObservationsOUT};3];
.qt.approx["ridgeCV: one RMSE per fold";{.quantQ.reg.ridgeCV[regYexact;regX;0f;5]};5#0f;1e-9];
.qt.true["ridgeCV: a large penalty increases the error";{(avg .quantQ.reg.ridgeCV[regYexact;regX;100f;5])>avg .quantQ.reg.ridgeCV[regYexact;regX;0f;5]}];

/ least-angle regression on y = 3 x1 - 2 x2 (+ small deterministic noise)
regN:40;
regX1:sin 0.37*til regN;
regX2:cos 1.11*til regN;
regTab:([] y:(3*regX1)-(2*regX2)-0.05*sin 2.9*til regN; x1:regX1; x2:regX2);
regLar:{.quantQ.reg.LAR[regTab;0.01]};
.qt.eq["LAR returns the path with arc length, intercept and one column per regressor";{cols regLar[]};`arcLength`beta0`beta1`beta2];
.qt.approx["LAR starts from beta = 0";{first each value flip select beta1, beta2 from regLar[]};0 0f;0f];
.qt.approx["LAR intercept is the mean of y";{first exec beta0 from regLar[]};avg regTab`y;1e-12];
.qt.true["LAR arc length is non-decreasing";{all 0<=1_deltas exec arcLength from regLar[]}];
.qt.true["LAR enters the most correlated regressor first";{p:regLar[]; (p[`beta1]1)<>0}];
/ LAR stops once no direction has correlation above 1/sqrt(n), so the end point is not the OLS fit
.qt.true["LAR ends with the signs of the true coefficients";{p:last regLar[]; (p[`beta1]>0) and p[`beta2]<0}];
.qt.true["LAR: both regressors are active at the end";{p:last regLar[]; (0<>p`beta1) and 0<>p`beta2}];
.qt.true["LARInOut adds an out-of-sample RMSE column";{`RMSE in cols .quantQ.reg.LARInOut[regTab;0.05;0.75]}];
.qt.true["LARInOut: the error falls along the path";{r:exec RMSE from .quantQ.reg.LARInOut[regTab;0.05;0.75]; (last r)<first r}];

.qt.finish[];
