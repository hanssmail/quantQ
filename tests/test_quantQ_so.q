/ Unit tests for lib/quantQ_so.q -- stochastic optimisation by random search on n-spheres (.quantQ.so)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`so;

soPi:acos -1;
soNorm:{[v] sqrt v wsum v};

/ spherical -> Euclidean coordinates
.qt.approx["sph2Eucl in 2-d, angle 0";{.quantQ.so.sph2Eucl (2f;enlist 0f)};2 0f;1e-12];
.qt.approx["sph2Eucl in 2-d, angle pi/2";{.quantQ.so.sph2Eucl (1f;enlist soPi%2)};0 1f;1e-12];
.qt.approx["sph2Eucl in 3-d (polar angle pi/2, azimuth pi/2)";{.quantQ.so.sph2Eucl (1f;(soPi%2;soPi%2))};0 0 1f;1e-12];
.qt.approx["sph2Eucl preserves the radius in 5-d";{soNorm .quantQ.so.sph2Eucl (1.5;0 0.2 1.3 2f)};1.5;1e-12];

/ random points on spheres
.qt.approx["randNSphere: unit 3-sphere";{soNorm .quantQ.so.randNSphere[1f;3]};1f;1e-12];
.qt.approx["randNSphere: radius 2 in 2-d";{soNorm .quantQ.so.randNSphere[2f;2]};2f;1e-12];
.qt.eq["randNSphere: dimension";{count .quantQ.so.randNSphere[1f;6]};6];
.qt.eq["nSphereNeighbours: number of points";{count .quantQ.so.nSphereNeighbours[2 2 2f;1f;5]};5];
.qt.approx["nSphereNeighbours: all at distance r from the pivot";{soNorm each .quantQ.so.nSphereNeighbours[2 2 2f;0.5;5]-\:2 2 2f};5#0.5;1e-12];

/ minimisation
soBucket:`func`x0!({x wsum x};1 2 3f);
.qt.true["minimise reduces a quadratic in 3-d";{r:last .quantQ.so.minimise soBucket; 1>r[`func] r`x0}];
.qt.true["minimise in 1-d";{r:last .quantQ.so.minimise `func`x0!({x wsum x};enlist 45f); 5>abs first r`x0}];
.qt.eq["minimise without monitoring returns the final state";{99h=type .quantQ.so.minimise soBucket,enlist[`monitorConvergence]!enlist 0};1b];
.qt.true["minimise respects maxEpochs";{r:.quantQ.so.minimise soBucket,`maxEpochs`nL2!(5;1000); 7>=count r}];
.qt.true["oneStep never makes the objective worse";{b:.quantQ.so.minimise soBucket,enlist[`maxEpochs]!enlist 20; f:{x wsum x} each exec x0 from b; all 0>=1_deltas f}];

.qt.finish[];
