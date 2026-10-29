/ Unit tests for lib/quantQ_dtw.q -- dynamic time warping (.quantQ.dtw)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`dtw;

/ reference values from a direct Python implementation of the DTW recursion with |a-b| cost
.qt.approx["optimalWarpingPath: identical series have distance 0";{.quantQ.dtw.optimalWarpingPath[1 2 3f;1 2 3f]};0f;1e-12];
.qt.approx["optimalWarpingPath: a repeated value costs nothing";{.quantQ.dtw.optimalWarpingPath[1 2 3f;1 2 2 3f]};0f;1e-12];
.qt.approx["optimalWarpingPath: 0 1 2 vs 0 2";{.quantQ.dtw.optimalWarpingPath[0 1 2f;0 2f]};1f;1e-12];
.qt.approx["optimalWarpingPath: 1 2 3 4 2 vs 1 3 4";{.quantQ.dtw.optimalWarpingPath[1 2 3 4 2f;1 3 4f]};3f;1e-12];
.qt.approx["optimalWarpingPath is symmetric";{.quantQ.dtw.optimalWarpingPath[1 3 4f;1 2 3 4 2f]};3f;1e-12];

.qt.approx["cost: full accumulated-cost matrix";{.quantQ.dtw.cost[1 2 3f;1 2 2 3f]};(0 0w 0w 0w 0w;0w 0 1 2 4;0w 1 0 0 1;0w 3 1 1 0f);1e-12];
.qt.approx["cost: last cell equals optimalWarpingPath";{last last .quantQ.dtw.cost[1 2 3 4 2f;1 3 4f]};3f;1e-12];

/ a window at least as wide as the series is the unconstrained problem
.qt.approx["optimalWarpingPathWindow with a wide window equals DTW";{.quantQ.dtw.optimalWarpingPathWindow[1 2 3 4 2f;1 3 4f;10]};3f;1e-12];
.qt.approx["optimalWarpingPathWindow with a wide window, equal lengths";{.quantQ.dtw.optimalWarpingPathWindow[0 1 2 1f;0 2 2 1f;10]};1f;1e-12];

/ warpingPath returns the 1-based indices of the optimal path
.qt.eq["warpingPath: identical series follow the diagonal";{.quantQ.dtw.warpingPath[1 2 3f;1 2 3f]};(1 2 3;1 2 3)];
.qt.eq["warpingPath: a repeated value in b is matched twice";{.quantQ.dtw.warpingPath[1 2 3f;1 2 2 3f]};(1 2 2 3;1 2 3 4)];

.qt.finish[];
