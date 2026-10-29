/ Unit tests for lib/quantQ_lof.q -- local outlier factor (.quantQ.lof)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`lof;

/ 1-d points 0 1 2 3 and an outlier at 10; observations are given row-wise in column x
lofData:([] x:enlist each 0 1 2 3 10f);

.qt.approx["L2distance: observations are given row-wise";{.quantQ.lof.L2distance[0 0f;(1 0f;3 4f)]};1 5f;1e-12];
.qt.approx["getKNNProperties: distance to the 2nd neighbour (the point itself counts)";{.quantQ.lof.getKNNProperties[2;lofData;enlist 0f]`distanceToK};1f;1e-12];
.qt.approx["getKNNProperties: number of points within that distance";{.quantQ.lof.getKNNProperties[2;lofData;enlist 0f]`numberToK};2;0f];
.qt.approx["getKNNProperties: ties are all counted";{.quantQ.lof.getKNNProperties[2;lofData;enlist 2f]`numberToK};3;0f];
.qt.approx["reachabilityDistance = max(d(X,Y), k-distance(Y))";{.quantQ.lof.reachabilityDistance[enlist 0f;4;lofData;2]};10f;1e-12];
.qt.approx["reachabilityDistance is at least the k-distance of Y";{.quantQ.lof.reachabilityDistance[enlist 1.1;1;lofData;2]};1f;1e-12];

/ LOF values evaluated in Python for this implementation's neighbourhood conventions
.qt.approx["LOFdata: LOF of every point";{exec LOF from .quantQ.lof.LOFdata[lofData;2]};1.5 1.125 1.125 1.5 7.5;1e-12];
.qt.eq["LOFdata: the outlier has the largest LOF";{exec first i from .quantQ.lof.LOFdata[lofData;2] where LOF=max LOF};4];
.qt.true["LOF of a point inside a uniform cluster is about 1";{l:.quantQ.lof.LOF[enlist 1.5;([] x:enlist each 0 1 2 3f);2]; (0.5<l) and l<1.5}];
.qt.true["localReachabilityDensity is lower far from the data";{.quantQ.lof.localReachabilityDensity[enlist 20f;lofData;2]<.quantQ.lof.localReachabilityDensity[enlist 1.5;lofData;2]}];

.qt.finish[];
