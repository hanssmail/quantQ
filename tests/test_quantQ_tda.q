/ Unit tests for lib/quantQ_tda.q -- topological data analysis: Vietoris-Rips complexes and loops (.quantQ.tda)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`tda;

/ point cloud: unit square (0..3) plus a far-away point (4)
tdaPts:update index:i from ([] point:(0 0f;1 0f;1 1f;0 1f;5 5f));
tdaMetric:.quantQ.tda.distance[.quantQ.tda.Lp[2]];
tdaEdges:.qt.fixture["tdaEdges";{.quantQ.tda.getAllEdges[tdaPts;tdaMetric]}];
tdaSet:{asc asc each x};
tdaL1:{sum abs x-y};

/ norms
.qt.approx["Lp, p=1";{.quantQ.tda.Lp[1;1 2 3f;3 -2 -10f]};19f;1e-12];
.qt.approx["Lp, p=2";{.quantQ.tda.Lp[2;1 2 3f;3 -2 -10f]};sqrt 189f;1e-12];
.qt.approx["Lp, p=inf";{.quantQ.tda.Lp[0w;1 2 3f;3 -2 -10f]};13f;1e-12];
.qt.approx["Lp, p=1.5";{.quantQ.tda.Lp[1.5;0 0f;1 1f]};xexp[2f;1%1.5];1e-12];
.qt.approx["distance wraps a metric";{.quantQ.tda.distance[tdaL1;1 2f;4 6f]};7f;1e-12];

/ combinatorics
.qt.eq["upperTriangTri lists all ordered triplets";{.quantQ.tda.upperTriangTri 10 20 30 40};(10 20 30;10 20 40;10 30 40;20 30 40)];

/ edges
.qt.eq["getAllEdges: one edge per pair of points";{count tdaEdges};10];
.qt.approx["getAllEdges: side and diagonal lengths";{exec distance from tdaEdges where edge in (0 1;0 2)};1,sqrt 2f;1e-12];

/ triangles and complexes
.qt.eq["getTriangles: no triangles among the square's sides";{count .quantQ.tda.getTriangles[tdaEdges;1.0]};0];
.qt.eq["getTriangles: with the diagonals the square has 4 triangles";{tdaSet .quantQ.tda.getTriangles[tdaEdges;1.5]};tdaSet (0 1 2;0 1 3;0 2 3;1 2 3)];
.qt.eq["getAllVR: the square becomes one 4-point complex";{asc each .quantQ.tda.getAllVR .quantQ.tda.getTriangles[tdaEdges;1.5]};enlist 0 1 2 3];
.qt.eq["getAllVRUnique: one 4-point complex";{asc each .quantQ.tda.getAllVRUnique .quantQ.tda.getTriangles[tdaEdges;1.5]};enlist 0 1 2 3];

/ loops
.qt.eq["getAllLoops: the square's sides form one loop";{count .quantQ.tda.getAllLoops exec edge from tdaEdges where distance<=1.0};1];
.qt.true["getAllLoops: the loop is closed and visits all 4 corners";{l:first .quantQ.tda.getAllLoops exec edge from tdaEdges where distance<=1.0; ((first l)=last l) and (asc distinct l)~0 1 2 3}];
.qt.eq["getAllPureLoops: the empty square is a pure loop";{count .quantQ.tda.getAllPureLoops[tdaEdges;1.0]};1];
.qt.eq["getAllPureLoops: filling the square with triangles kills the loop";{count .quantQ.tda.getAllPureLoops[tdaEdges;1.5]};0];
.qt.eq["getAllPureLoops_NaiveFilter agrees";{count each .quantQ.tda.getAllPureLoops_NaiveFilter[tdaEdges;] each 1.0 1.5};1 0];

/ summaries over thresholds
.qt.eq["persistentGraphPoint: nodes attached to an edge";{exec deathPoint from .quantQ.tda.persistentGraphPoint[tdaEdges;1.0 1.5 10]};4 4 5];
.qt.eq["nComplexesHistThr: one 4-point complex at 1.5";{select orderVR, numberVR from .quantQ.tda.nComplexesHistThr[tdaEdges;enlist 1.5]};([] orderVR:enlist 4; numberVR:enlist 1)];
.qt.eq["nComplexesHistThr: placeholder row when there are no triangles";{select orderVR, numberVR from .quantQ.tda.nComplexesHistThr[tdaEdges;enlist 1.0]};([] orderVR:enlist 0; numberVR:enlist 0)];
.qt.eq["nLoopsThr";{exec loopCount from .quantQ.tda.nLoopsThr[tdaEdges;enlist 1.0]};enlist 1];
.qt.true["bornAndDeathLoops: the square's loop is born at 1.0 and dies before 1.5";{t:0!.quantQ.tda.bornAndDeathLoops[tdaEdges;1.0 1.2 1.5]; (1=count t) and (1f=first t`born) and 1.2=first t`death}];

.qt.finish[];
