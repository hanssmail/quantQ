/ Unit tests for lib/quantQ_meb.q -- maximum-entropy (block) bootstrap (.quantQ.meb)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`meb;

mebX:0.23 0.24 0.12 0.09 0.02 0.45;
mebLong:{(sin 0.2*x)+0.01*x} "f"$til 60;
mebRank:{iasc iasc x};

/ original Vinod-Lopez maximum entropy bootstrap
.qt.eq["BSorig keeps the length";{count .quantQ.meb.BSorig mebX};6];
.qt.eq["BSorig preserves the rank order of the original series";{mebRank .quantQ.meb.BSorig mebX};mebRank mebX];
.qt.true["BSorig draws within the trimmed range";{b:.quantQ.meb.BSorig mebX; m:avg abs 1_deltas asc mebX; all (b>=min[mebX]-m) and b<=max[mebX]+m}];

/ exponential-tail version
.qt.eq["BS keeps the length";{count .quantQ.meb.BS[()!();mebX]};6];
.qt.eq["BS preserves the rank order";{mebRank .quantQ.meb.BS[()!();mebX]};mebRank mebX];
.qt.eq["BS preserves the rank order of a longer series";{mebRank .quantQ.meb.BS[()!();mebLong]};mebRank mebLong];
.qt.true["BS replicates differ from each other";{not .quantQ.meb.BS[()!();mebLong]~.quantQ.meb.BS[()!();mebLong]}];

/ block shuffles
.qt.eq["blockShuffle keeps the length";{count .quantQ.meb.blockShuffle[()!();mebLong]};60];
.qt.true["blockShuffle only reuses observed values";{all .quantQ.meb.blockShuffle[()!();mebLong] in mebLong}];
.qt.eq["blockShuffleUni keeps the length";{count .quantQ.meb.blockShuffleUni[()!();mebLong]};60];
.qt.true["blockShuffleUni without sewing and with tiny noise stays on the data";{b:.quantQ.meb.blockShuffleUni[`block`uniSpread`sew!(5;1e-12;0b);mebLong]; all {[d;v] 1e-9>min abs d-v}[mebLong] each b}];
.qt.true["blockShuffleUni: blocks are runs of consecutive observations";{b:.quantQ.meb.blockShuffleUni[`block`uniSpread`sew!(5;1e-12;0b);mebLong]; i:{[d;v] first where 1e-9>abs d-v}[mebLong] each b; all raze 1=1_'deltas each 5 cut i}];
.qt.eq["blockMeb keeps the length";{count .quantQ.meb.blockMeb[()!();mebLong]};60];
.qt.true["blockMeb returns finite floats";{b:.quantQ.meb.blockMeb[()!();mebLong]; (9h=type b) and (not any null b) and all (abs b)<0w}];

.qt.finish[];
