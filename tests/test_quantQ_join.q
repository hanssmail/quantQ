/ Unit tests for lib/quantQ_join.q -- helper statistics used with joins (.quantQ.join)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`join;

.qt.approx["asymmetry = fraction of ones";{.quantQ.join.asymmetry 1 0 1 1};0.75;1e-12];
.qt.approx["asymmetry of no ones is 0";{.quantQ.join.asymmetry 0 0 -1};0f;1e-12];
.qt.approx["tStatsCorr: rho*sqrt((N-2)/(1-rho^2))";{.quantQ.join.tStatsCorr[0.5;27]};2.886751345948129;1e-12];
.qt.approx["tStatsCorr of zero correlation is 0";{.quantQ.join.tStatsCorr[0f;100]};0f;1e-12];
.qt.approx["tStatsCorr is odd in rho";{.quantQ.join.tStatsCorr[-0.5;27]};-2.886751345948129;1e-12];

.qt.finish[];
