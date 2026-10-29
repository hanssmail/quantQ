/ Unit tests for lib/quantQ_trees.q -- decision trees (ID3-style information gain) (.quantQ.trees)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`trees;

/ helpers
.qt.eq["classify bins values against sorted breakpoints";{.quantQ.trees.classify[3 1 2f;1.5 2.5 0.5 3.5]};1 2 0 3];
.qt.approx["entropy of a 50/50 split is 1 bit";{.quantQ.trees.entropy[0 0 1 1;0 1]};1f;1e-12];
.qt.approx["entropy of a pure node is 0";{.quantQ.trees.entropy[0 0 0 0;0 1]};0f;1e-12];
.qt.approx["entropy of three equiprobable classes";{.quantQ.trees.entropy[0 1 2;0 1 2]};2 xlog 3;1e-12];
.qt.approx["infogain of a perfect split is the parent entropy";{.quantQ.trees.infogain[0 0 1 1;(0 0;1 1);0 1]};1f;1e-12];
.qt.approx["infogain of an uninformative split is 0";{.quantQ.trees.infogain[0 0 1 1;(0 1;0 1);0 1]};0f;1e-12];
.qt.eq["runRule applies rule[x i; j]";{.quantQ.trees.runRule[>;1;2.5;1 3 2f]};1b];
.qt.eq["applyRule returns the two child bitmaps";{(.quantQ.trees.applyRule[4;0 1 2 3;1 2 3 4f;>;2.5])`bitmap};(1100b;0011b)];

/ learning a tree on a one-feature, separable problem
treeP:`x`y`classes`rule!(enlist 1 2 3 4 5 6f;0 0 0 1 1 1;0 1;>);
treeT:{.quantQ.trees.learnTree treeP};
.qt.true["learnTree returns a table with the path columns";{all `i`p`path in cols treeT[]}];
.qt.true["learnTree: the root sees every observation";{all first[treeT[]]`bitmap}];
.qt.true["learnTree: every leaf is pure";{t:treeT[]; l:.quantQ.trees.leaves t; all {1=count distinct treeP[`y] where x} each l`bitmap}];
.qt.approx["learnTree: the first split has information gain 1 bit";{(treeT[]`infogain)1};1f;1e-12];
.qt.eq["predictOnTree: a point on the left";{l:.quantQ.trees.predictOnTree[treeT[];enlist 1.5]; distinct treeP[`y] where first l`bitmap};enlist 0];
.qt.eq["predictOnTree: a point on the right";{l:.quantQ.trees.predictOnTree[treeT[];enlist 5.5]; distinct treeP[`y] where first l`bitmap};enlist 1];
.qt.true["learnTree respects maxdepth";{t:.quantQ.trees.learnTree treeP,enlist[`maxdepth]!enlist 1; all 1>=count each t`rulepath}];

.qt.finish[];
