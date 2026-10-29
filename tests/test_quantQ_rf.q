/ Unit tests for lib/quantQ_rf.q -- random forests built on .quantQ.trees (.quantQ.rf)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`rf;

/ two features, the class is decided by the first one
rfX1:0.1*til 20;
rfX2:(20#0 1f);
rfY:?[rfX1>0.95;1;0];
rfS:`x`y!((rfX1;rfX2);rfY);

/ bootstrap sample
rfSmp:.qt.fixture["rfSmp";{.quantQ.rf.sampleTree[rfS;20]}];
.qt.eq["sampleTree returns x, y and the in/out-of-bag indices";{key rfSmp};`x`y`oobi`ibi];
.qt.eq["sampleTree: bootstrap sample has n draws";{count rfSmp`ibi};20];
.qt.true["sampleTree: out-of-bag indices are exactly those not drawn";{(asc rfSmp`oobi)~asc (til 20) except rfSmp`ibi}];
.qt.true["sampleTree: y is resampled with x";{rfSmp[`y]~rfY rfSmp`ibi}];

/ forest
rfPar:`x`y`classes`rule`m`n`B`maxdepth!((rfX1;rfX2);rfY;0 1;>;1;20;5;5);
.qt.eq["randomForest: one ibi entry per tree";{count (.quantQ.rf.randomForest rfPar)`ibi};5];
.qt.true["randomForest: trees are tagged with their bootstrap index";{asc[distinct exec B from (.quantQ.rf.randomForest rfPar)`tree]~til 5}];
.qt.eq["predictOnRF: a clear class-0 point";{(.quantQ.rf.predictOnRF[rfY;.quantQ.rf.randomForest rfPar;0.2 0f])`prediction};0];
.qt.eq["predictOnRF: a clear class-1 point";{(.quantQ.rf.predictOnRF[rfY;.quantQ.rf.randomForest rfPar;1.8 1f])`prediction};1];
.qt.eq["predictOnRF reports the out-of-bag error";{key .quantQ.rf.predictOnRF[rfY;.quantQ.rf.randomForest rfPar;1.8 1f]};`prediction`mean_error`dev_error];

.qt.finish[];
