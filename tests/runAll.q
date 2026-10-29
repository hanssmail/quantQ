/ runAll.q -- run every test_quantQ_*.q file in this directory and print a summary.
/   q tests/runAll.q          (from the repository root)
/   cd tests && q runAll.q    (from this directory)
/ Exits with status 0 if all tests pass, 1 otherwise.

system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q";

.qt.testFiles:asc {x where x like "test_quantQ_*.q"} string key hsym `$.qt.dir;

/ a test file that fails to load (e.g. a top-level error) is recorded as a failure and the run continues
{[f]
    -1 "";
    @[system;"l ",.qt.dir,"/",f;{[f;e] .qt.module:`$f; .qt.record["file loads";0b;"error while loading: ",e]}[f]]
 } each .qt.testFiles;

.qt.nFailed:.qt.summary[];
exit $[.qt.nFailed>0;1;0];
