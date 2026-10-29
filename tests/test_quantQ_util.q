/ Unit tests for lib/quantQ_util.q -- utility functions (.quantQ.util)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`util;

utilTab:([] a1:til 5; a2:10+til 5; a3:20+til 5);

.qt.eq["selectCol with one column";{.quantQ.util.selectCol[utilTab;`a1]};([] a1:til 5)];
.qt.eq["selectCol with several columns";{.quantQ.util.selectCol[utilTab;`a3`a1]};([] a3:20+til 5; a1:til 5)];
.qt.eq["selectCol by reference";{.quantQ.util.selectCol[`utilTab;`a2]};([] a2:10+til 5)];
.qt.eq["deleteCol with one column";{.quantQ.util.deleteCol[utilTab;`a1]};([] a2:10+til 5; a3:20+til 5)];
.qt.eq["deleteCol with several columns";{.quantQ.util.deleteCol[utilTab;`a1`a3]};([] a2:10+til 5)];
.qt.eq["prev[n] lags by n with nulls";{.quantQ.util.prev[2] til 5};0N 0N 0 1 2];
.qt.eq["next[n] leads by n with nulls";{.quantQ.util.next[2] til 5};2 3 4 0N 0N];
.qt.eq["prev[0] is the identity";{.quantQ.util.prev[0] til 5};til 5];

.qt.finish[];
