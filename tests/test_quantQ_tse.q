/ Unit tests for lib/quantQ_tse.q -- AR simulation with scan (.quantQ.tse)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`tse;

tseEps:1 0.5 -0.3 0.2 0 0.4;

/ x_t = phi1 x_(t-1) + ... + phip x_(t-p) + eps_t, starting from zeros
.qt.approx["simAR, AR(2)";{.quantQ.tse.simAR[0.5 0.2;tseEps]};1 1 0.4 0.6 0.38 0.71;1e-12];
.qt.approx["simAR, AR(1)";{.quantQ.tse.simAR[enlist 0.5;1 0 0 0f]};1 0.5 0.25 0.125;1e-12];
.qt.approx["simAR with phi=0 returns the shocks";{.quantQ.tse.simAR[0 0f;tseEps]};tseEps;1e-12];
.qt.eq["simAR keeps the length of the shocks";{count .quantQ.tse.simAR[0.5 0.2 0.1;100?1f]};100];
/ the README says this is equivalent to .quantQ.ts.simAR; .quantQ.ts.simAR starts p periods later
.qt.approx["simAR matches .quantQ.ts.simAR up to the initial p zeros";{.quantQ.tse.simAR[0.5 0.2;tseEps]};2_.quantQ.ts.simAR[0.5 0.2;tseEps,0 0f];1e-12];

.qt.finish[];
