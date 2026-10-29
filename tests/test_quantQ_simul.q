/ Unit tests for lib/quantQ_simul.q -- simulation: Cholesky, Box-Muller, multivariate normals (.quantQ.simul)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`simul;

/ intraday volatility pattern
simPar:`tIni`tFin`c1`cOpen`aOpen`cClose`aClose!(09:30:00.000;16:00:00.000;0.5;1.2;3f;0.8;2f);
.qt.approx["stochVolIntraPatt is 1 when the pattern is switched off";{.quantQ.simul.stochVolIntraPatt[12:00:00.000;0b;simPar]};1f;0f];
.qt.approx["stochVolIntraPatt at the open";{.quantQ.simul.stochVolIntraPatt[09:30:00.000;1b;simPar]};0.5+1.2+0.8*exp -2f;1e-12];

/ Cholesky decomposition
simA:(4 2f;2 3f);
simA3:(4 2 0.4;2 3 0.5;0.4 0.5 2f);
.qt.approx["choleskyChB of a 2x2 matrix";{.quantQ.simul.choleskyChB simA};(2 0f;1,sqrt 2f);1e-12];
.qt.approx["choleskyChB: L L' = A (3x3)";{l:.quantQ.simul.choleskyChB simA3; l mmu flip l};simA3;1e-12];
.qt.approx["choleskyChB: L is lower triangular";{l:.quantQ.simul.choleskyChB simA3; (l[0;1 2];l[1;2])};0 0 0f;0f];

/ Box-Muller
.qt.approx["radius is the squared norm";{.quantQ.simul.radius 3 4f};25f;1e-12];
.qt.approx["BoxMuller transform of a point in the unit disc";{.quantQ.simul.BoxMuller 0.3 0.4};0.3 0.4*sqrt -2*(log 0.25)%0.25;1e-12];
.qt.eq["genBoxMuller returns a pair of floats";{(type r;count r:.quantQ.simul.genBoxMuller[])};(9h;2)];
.qt.eq["genBoxMuller1 returns a pair of floats";{(type r;count r:.quantQ.simul.genBoxMuller1[])};(9h;2)];
.qt.eq["getNormalVariate: odd lengths are supported";{count .quantQ.simul.getNormalVariate 7};7];
simZ:.qt.fixture["simZ";{.quantQ.simul.getNormalVariate 20000}];
.qt.true["getNormalVariate: sample mean ~ 0";{0.03>abs avg simZ}];
.qt.true["getNormalVariate: sample sd ~ 1";{0.03>abs 1-dev simZ}];
.qt.true["getNormalVariate: ~95% within 1.96";{0.01>abs 0.95-avg 1.96>abs simZ}];

/ correlated normals
simMN:.qt.fixture["simMN";{.quantQ.simul.multiNormal[2;20000;1 -2f;(1 0.5;0.5 2f)]}];
.qt.eq["multiNormal: columns var1..varN";{cols simMN};`var1`var2];
.qt.true["multiNormal: means";{all 0.05>abs (avg each value flip simMN)-1 -2f}];
.qt.true["multiNormal: variances";{all 0.08>abs (var each value flip simMN)-1 2f}];
.qt.true["multiNormal: covariance";{0.05>abs 0.5-cov[simMN`var1;simMN`var2]}];

/ exponential variates
simE:.qt.fixture["simE";{.quantQ.simul.expRandomVariate each 20000#2f}];
.qt.true["expRandomVariate is positive";{all simE>0}];
.qt.true["expRandomVariate has mean 1/lambda";{0.02>abs 0.5-avg simE}];

.qt.finish[];
