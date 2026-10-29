/ Unit tests for lib/quantQ_io.q -- splayed / partitioned table I/O (.quantQ.io)
/ Writes to a scratch directory next to the tests (removed at the end; Linux/macOS shell).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`io;

ioDir:.qt.dir,"/qtest_tmp_io";
ioPath:hsym `$ioDir;
@[system;"rm -rf ",ioDir;::];

ioTab:([] sym:`c`a`b`a`c; px:5 1 2 3 4f; qty:50 10 20 30 40);
ioMore:([] sym:`b`a; px:6 7f; qty:60 70);

.qt.true["saveSplayedTab writes the table";{.quantQ.io.saveSplayedTab[ioPath;`sym;`ioTab]; 0<count key ` sv ioPath,`ioTab}];
.qt.true["saveSplayedTab: reloaded rows are sorted by the parted column";{t:get ` sv ioPath,`ioTab,`; (`$string t`sym)~`a`a`b`c`c}];
.qt.true["saveSplayedTab: values survive the round trip";{t:get ` sv ioPath,`ioTab,`; e:`sym xasc ioTab; (t[`px]~e`px) and t[`qty]~e`qty}];
.qt.eq["saveSplayedTab applies the parted attribute";{attr (get ` sv ioPath,`ioTab,`)`sym};`p];
.qt.true["appendSplayedTab appends and keeps the table sorted";{.quantQ.io.appendSplayedTab[ioPath;`sym;`ioTab;ioMore]; t:get ` sv ioPath,`ioTab,`; (7=count t) and (`$string t`sym)~`a`a`a`b`b`c`c}];

/ dpftsAppend: like .Q.dpft, but amends the partition with a given operator (here: append)
ioTrade:([] sym:`x`y`x; px:1 2 3f);
.qt.eq["dpftsAppend writes a partition";{.quantQ.io.dpftsAppend[ioPath;2020.01.01;`sym;`ioTrade;`sym;,]; count get ` sv ioPath,`2020.01.01,`ioTrade,`};3];
.qt.eq["dpftsAppend with ',' appends to an existing partition";{.quantQ.io.dpftsAppend[ioPath;2020.01.01;`sym;`ioTrade;`sym;,]; count get ` sv ioPath,`2020.01.01,`ioTrade,`};6];

@[system;"rm -rf ",ioDir;::];
.qt.finish[];
