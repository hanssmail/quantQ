/ Unit tests for lib/quantQ_apriori.q -- Apriori frequent itemsets (.quantQ.apriori)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`apriori;

/ makeDictionary
.qt.eq["makeDictionary maps each value to its one-hot vector";{.quantQ.apriori.makeDictionary `a`b`c};`a`b`c!(100b;010b;001b)];

/ createDummies
aprTab:([] id:1 2 3; col:`x`y`x);
.qt.eq["createDummies: column names are <col>_<k>";{cols .quantQ.apriori.createDummies[aprTab;0]};`col_0`col_1];
.qt.eq["createDummies: one-hot encoding of the sorted distinct values";{.quantQ.apriori.createDummies[aprTab;0]};([] col_0:101b; col_1:010b)];

/ isIn
.qt.true["isIn: every item of the first list is in the second";{.quantQ.apriori.isIn[`a`b;`a`b`c]}];
.qt.true["isIn: false when an item is missing";{not .quantQ.apriori.isIn[`a`d;`a`b`c]}];

/ aprioriOneRun -- baskets: a in 3/4, b in 2/4, c in 1/4, {a,b} in 2/4; threshold 40%
aprDummies:([] id:1 2 3 4; a:1110b; b:1100b; c:0001b);
aprInit:`step`largeTabDummies`thresh`colDummiesStepPrev`colDummiesStep1`colDummiesStepALL`colDummiesStepALLX`colDummiesStepPrevX`canImprove!(0;aprDummies;0.4;();();();();();1b);
.qt.eq["aprioriOneRun step 1: frequent single items";{(.quantQ.apriori.aprioriOneRun aprInit)`colDummiesStep1};(enlist `a;enlist `b)];
.qt.eq["aprioriOneRun step 1: infrequent single items";{(.quantQ.apriori.aprioriOneRun aprInit)`colDummiesStepALLX};enlist enlist `c];
.qt.true["aprioriOneRun step 1: can improve";{(.quantQ.apriori.aprioriOneRun aprInit)`canImprove}];
.qt.eq["aprioriOneRun step 2: frequent pairs";{(.quantQ.apriori.aprioriOneRun .quantQ.apriori.aprioriOneRun aprInit)`colDummiesStepPrev};enlist `a`b];
.qt.true["aprioriOneRun step 2: {a,b} added to the frequent itemsets";{any (`a`b)~/:(.quantQ.apriori.aprioriOneRun .quantQ.apriori.aprioriOneRun aprInit)`colDummiesStepALL}];

.qt.finish[];
