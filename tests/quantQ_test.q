/ quantQ_test.q -- a small, dependency-free unit-test harness for the quantQ library.
/ ---
/ Every test_quantQ_<module>.q file loads this harness (once) and then the whole
/ library from ../lib. Run a single module:   q tests/test_quantQ_mat.q
/ Run the whole suite:                      q tests/runAll.q
/ Both exit with status 0 when every test passes and 1 otherwise.
/ ---
/ Assertions (each takes a test name and a niladic lambda, so errors are trapped):
/   .qt.eq[name; {expr}; expected]          -- expr ~ expected
/   .qt.approx[name; {expr}; expected; tol] -- numeric match, element-wise |a-e|<=tol
/                                              (works on atoms, vectors, matrices,
/                                              dictionaries and tables; nulls must
/                                              line up)
/   .qt.true[name; {expr}]                  -- expr returns 1b
/   .qt.throws[name; {expr}]                -- expr signals an error
/ Fixtures that call the library: x:.qt.fixture["x";{expr}] (records an error instead of aborting).

/ ---------------------------------------------------------------------------
/ locations
.qt.dirOf:{[s] $[count w:where s="/";(last w)#s;"."]};
.qt.dir:.qt.dirOf string .z.f;
.qt.exists:{[p] 0<count key hsym `$p};
.qt.libCandidates:{x where .qt.exists each x,\:"/quantQ_mat.q"} (.qt.dir,"/../lib";"lib";"../lib";"quantQ/lib");
if[0=count .qt.libCandidates; -2 "quantQ tests: cannot find the lib directory (looked next to ",.qt.dir,")"; exit 2];
.qt.libDir:first .qt.libCandidates;
.qt.isRunAll:.z.f like "*runAll.q";

/ ---------------------------------------------------------------------------
/ results
.qt.results:([] module:`symbol$(); name:(); ok:`boolean$(); msg:());
.qt.module:`;

.qt.setModule:{[m]
    .qt.module:m;
    / every module starts from the same random seed so stochastic tests are repeatable
    system "S 42";
    -1 "== ",string m;
 };

.qt.str:{[x] s:-3!x; $[200<count s;(200#s),"...";s]};
.qt.errStr:{[e] $[10h=abs type e;e;.qt.str e]};

.qt.record:{[name;ok;msg]
    `.qt.results upsert ([] module:enlist .qt.module; name:enlist name; ok:enlist ok; msg:enlist msg);
    if[not ok; -1 "  FAIL ",name," :: ",msg];
    ok
 };

/ evaluate a niladic lambda; returns (1b;result) or (0b;error message)
.qt.run:{[f] @[{(1b;x[])};f;{(0b;x)}]};

/ ---------------------------------------------------------------------------
/ numeric helpers

/ flatten atoms, vectors, nested lists, dictionaries and tables into one list
.qt.flat:{[x]
    t:type x;
    $[t<0h; enlist x;
      t within 1 19h; x;
      t=98h; .z.s value flip x;
      t=99h; .z.s value x;
      0=count x; x;
      raze .z.s each x]
 };

.qt.close:{[a;e;tol]
    a:"f"$.qt.flat a;
    e:"f"$.qt.flat e;
    if[not count[a]=count e; :0b];
    if[not (null a)~null e; :0b];
    w:where not null e;
    all tol>=abs a[w]-e[w]
 };

/ ---------------------------------------------------------------------------
/ assertions
.qt.eq:{[name;f;expected]
    r:.qt.run f;
    if[not first r; :.qt.record[name;0b;"error: ",.qt.errStr last r]];
    ok:expected~last r;
    .qt.record[name;ok;$[ok;"";"expected ",.qt.str[expected],"  got ",.qt.str last r]]
 };

.qt.approx:{[name;f;expected;tol]
    r:.qt.run f;
    if[not first r; :.qt.record[name;0b;"error: ",.qt.errStr last r]];
    ok:@[{.qt.close . x};(last r;expected;tol);0b];
    .qt.record[name;ok;$[ok;"";"expected ",.qt.str[expected]," (tol ",string[tol],")  got ",.qt.str last r]]
 };

.qt.true:{[name;f]
    r:.qt.run f;
    if[not first r; :.qt.record[name;0b;"error: ",.qt.errStr last r]];
    ok:1b~last r;
    .qt.record[name;ok;$[ok;"";"expected 1b  got ",.qt.str last r]]
 };

.qt.throws:{[name;f]
    r:.qt.run f;
    .qt.record[name;not first r;$[first r;"expected an error, got ",.qt.str last r;""]]
 };

/ build a fixture without aborting the test file: on error the failure is recorded and the
/ error string is returned, so the tests that use the fixture fail instead of the whole file
.qt.fixture:{[name;f]
    r:.qt.run f;
    if[not first r; .qt.record["fixture ",name;0b;"error: ",.qt.errStr last r]];
    last r
 };

/ ---------------------------------------------------------------------------
/ reporting
.qt.summary:{[]
    -1 "";
    show select tests:count i, passed:sum ok, failed:sum not ok by module from .qt.results;
    nf:sum not .qt.results`ok;
    -1 "";
    -1 string[count .qt.results]," tests, ",string[nf]," failed";
    if[nf>0;
        -1 "";
        -1 "Failures:";
        {-1 "  ",string[x`module],": ",x[`name]," :: ",x`msg} each select from .qt.results where not ok
    ];
    nf
 };

/ called at the end of each test file: prints a summary and exits when the file is run on its own
.qt.finish:{[] if[not .qt.isRunAll; nf:.qt.summary[]; exit $[nf>0;1;0]]; };

/ ---------------------------------------------------------------------------
/ load the library (every quantQ_*.q file, in the same order as lib/quantQ.q)
.qt.loadLib:{[]
    files:asc {x where x like "quantQ_*.q"} string key hsym `$.qt.libDir;
    {system "l ",.qt.libDir,"/",x} each files;
 };
.qt.loadLib[];
