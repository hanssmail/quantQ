/ Unit tests for lib/quantQ_ols.q -- ordinary least squares and diagnostics (.quantQ.ols)
/ Reference values from numpy (lstsq, classical OLS standard errors).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`ols;

/ y = 1 + 2x + e with a fixed perturbation e
olsXv:1 2 3 4 5 6 7 8 9 10f;
olsY:3.3 4.8 7.1 8.6 11.25 12.95 15.15 16.7 19.2 21.0;
olsX:(10#1f;olsXv);
olsFit:{.quantQ.ols.fit[olsY;olsX]};

/ significance labels from |t|
.qt.eq["signLevels";{.quantQ.ols.signLevels 0.5 1.5 1.7 2.5 3.5 -2.0};("-";"* (90% c.l.)";"** (95% c.l.)";"*** (99% c.l.)";"**** (99.9% c.l.)";"** (95% c.l.)")];

/ goodness of fit
olsYhat:1.0233333333333308+1.9966666666666677*olsXv;
.qt.approx["rSquared";{.quantQ.ols.rSquared[olsY;olsYhat]};0.9984325856240578;1e-12];
.qt.approx["rSquared of a perfect fit is 1";{.quantQ.ols.rSquared[olsY;olsY]};1f;1e-12];
.qt.approx["rSquaredAdj";{.quantQ.ols.rSquaredAdj[olsY;olsYhat;1]};0.998236658827065;1e-12];
.qt.approx["logL (Gaussian, MLE variance)";{.quantQ.ols.logL[10;0.5163333333333354]};0.6285537690903737;1e-10];
.qt.approx["AIC = 2p - 2logL";{.quantQ.ols.AIC[2;0.6285537690903737]};2.7428924618192525;1e-12];
.qt.approx["AICC adds the small-sample correction";{.quantQ.ols.AICC[2;10;0.6285537690903737]};4.457178176104967;1e-12];

/ fit
.qt.approx["fit: coefficients";{olsFit[]`beta};1.0233333333333308 1.9966666666666677;1e-10];
.qt.approx["fit: standard errors";{olsFit[]`errors};0.17354954463911615 0.027970041693269523;1e-10];
.qt.approx["fit: t statistics";{olsFit[]`tStats};5.896491030623441 71.38590240811578;1e-8];
.qt.approx["fit: covariance matrix";{olsFit[]`varcov};(0.030119444444444567 -0.004302777777777795;-0.004302777777777795 0.0007823232323232354);1e-12];
.qt.approx["fit: TSS";{olsFit[]`TSS};329.41724999999997;1e-9];
.qt.approx["fit: RSS";{olsFit[]`RSS};0.5163333333333354;1e-10];
.qt.approx["fit: R squared";{olsFit[]`Rsquared};0.9984325856240578;1e-12];
.qt.approx["fit: adjusted R squared";{olsFit[]`adjRsquared};0.998236658827065;1e-12];
.qt.approx["fit: log-likelihood";{olsFit[]`logL};0.6285537690903737;1e-10];
.qt.approx["fit: AIC";{olsFit[]`aic};2.7428924618192525;1e-10];
.qt.approx["fit: AICC";{olsFit[]`aicc};4.457178176104967;1e-10];
.qt.eq["fit: observations and parameters";{olsFit[][`nObservations`dfn]};10 2];
.qt.eq["fit: significance labels";{olsFit[]`signLevels};("**** (99.9% c.l.)";"**** (99.9% c.l.)")];
.qt.eq["fit drops observations with a missing y";{.quantQ.ols.fit[@[olsY;3;:;0n];olsX]`nObservations};9];
.qt.eq["fit drops observations with a missing regressor";{.quantQ.ols.fit[olsY;(10#1f;@[olsXv;0 5;:;0n])]`nObservations};8];

/ tabular output
olsTabIn:([] y:olsY; const:10#1f; xv:olsXv);
.qt.eq["olsTab: one row per regressor plus summary rows";{count .quantQ.ols.olsTab olsTabIn};12];
.qt.approx["olsTab: coefficients";{2#exec coeff from .quantQ.ols.olsTab olsTabIn};1.0233333333333308 1.9966666666666677;1e-10];

/ nested-model F test (constant-only vs full)
.qt.approx["fStatistics";{.quantQ.ols.fStatistics[.quantQ.ols.fit[olsY;enlist 10#1f];olsFit[]]};5095.947062621025;1e-6];

/ prediction and errors
olsRm1:([] y:1 2 3 4f; yhat:0 3 2 5f);
olsRm2:([] y:1 2 3f; yhat:1 2 3f);
.qt.approx["predict: fitted values";{exec yPredicted from .quantQ.ols.predict[olsFit[];([] const:10#1f; xv:olsXv)]};olsYhat;1e-10];
.qt.approx["RMSE of residuals +-1 is 1";{.quantQ.ols.RMSE olsRm1};1f;1e-12];
.qt.approx["RMSE of a perfect prediction is 0";{.quantQ.ols.RMSE olsRm2};0f;1e-12];

/ in/out-of-sample and cross-validation on an exact line y = 1 + 2x
olsYexact:1+2*olsXv;
.qt.approx["inOut: exact data are predicted out-of-sample with zero error";{.quantQ.ols.inOut[olsYexact;olsX;0.6]`RMSE};0f;1e-9];
.qt.eq["inOut: out-of-sample size";{.quantQ.ols.inOut[olsYexact;olsX;0.6]`nObservationsOUT};4];
.qt.approx["inOut with a boolean out-of-sample flag";{.quantQ.ols.inOut[olsYexact;olsX;0101010101b]`RMSE};0f;1e-9];
.qt.approx["cv: one RMSE per fold, all zero on exact data";{.quantQ.ols.cv[olsYexact;olsX;5]};5#0f;1e-9];

/ R squared by group
olsGrp:([] y:olsY,olsY; yPred:olsYhat,olsY+0.5*(-1 1)10#0 1; groupBy:(10#0),10#1);
.qt.eq["byR2 has one row per group";{count .quantQ.ols.byR2 olsGrp};2];
.qt.eq["byR2 sorts groups by their R squared";{exec groupBy from .quantQ.ols.byR2 olsGrp};1 0];

.qt.finish[];
