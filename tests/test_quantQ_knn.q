/ Unit tests for lib/quantQ_knn.q -- k-nearest neighbours, prototypes, kernel regression (.quantQ.knn)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`knn;

/ distanceLX: the sample is given column-wise (one vector per feature)
.qt.approx["distanceLX L2";{.quantQ.knn.distanceLX[0 0f;(1 3f;0 4f);2]};1 5f;1e-12];
.qt.approx["distanceLX L1";{.quantQ.knn.distanceLX[0 0f;(1 3f;0 4f);1]};1 7f;1e-12];

knnTab:([] x1:1 3 0.5 -2f; x2:0 4 0.5 0f; f:`a`b`a`b);
.qt.eq["kNNneighbours returns the k closest classes";{(.quantQ.knn.kNNneighbours[0 0f;knnTab;2;2])`f};`a`a];
.qt.approx["kNNneighbours returns the sorted distances";{(.quantQ.knn.kNNneighbours[0 0f;knnTab;2;2])`dist};(sqrt 0.5;1f);1e-12];
.qt.eq["kNNclassifyMajorityVote";{.quantQ.knn.kNNclassifyMajorityVote `f`dist!(`a`b`b;0.1 0.5 0.6)};`b];
.qt.eq["kNNclassifyMajorityVote on the neighbours";{.quantQ.knn.kNNclassifyMajorityVote .quantQ.knn.kNNneighbours[0 0f;knnTab;2;3]};`a];
.qt.eq["kNNclassifyDistanceWeighted";{.quantQ.knn.kNNclassifyDistanceWeighted `f`dist!(`a`b`b;0.1 0.5 0.6)};`b];

/ prototype kNN: two well separated clusters
knnClusters:([] x:(0 0f;0.1 0f;0 0.1;0.1 0.1;5 5f;5.1 5f;5 5.1;5.1 5.1); f:`a`a`a`a`b`b`b`b);
.qt.eq["prototypeKNN: R centres per class";{count .quantQ.knn.prototypeKNN[knnClusters;2;1;3]};2];
.qt.eq["prototypeKNN: one centre for each class";{asc exec f from .quantQ.knn.prototypeKNN[knnClusters;2;1;3]};`a`b];
.qt.true["prototypeKNN: centres stay inside their cluster";{c:.quantQ.knn.prototypeKNN[knnClusters;2;1;3]; all 0.2>abs raze (exec centre from c)-(0.05 0.05;5.05 5.05)(`a`b?exec f from c)}];
.qt.eq["getPrototypeKNN classifies a point near cluster a";{.quantQ.knn.getPrototypeKNN[.quantQ.knn.prototypeKNN[knnClusters;2;1;3];0.2 0.3]};`a];
.qt.eq["getPrototypeKNN classifies a point near cluster b";{.quantQ.knn.getPrototypeKNN[.quantQ.knn.prototypeKNN[knnClusters;2;1;3];4.5 5.5]};`b];

/ weighted kernel regression
.qt.approx["weightedL2norm with unit weights";{.quantQ.knn.weightedL2norm[1 2f;0 0f;1 1f]};5f;1e-12];
.qt.approx["weightedL2norm: weights scale the coordinates";{.quantQ.knn.weightedL2norm[1 2f;0 0f;2 0f]};4f;1e-12];
.qt.approx["kNNMeanFunction: exp(-d/beta)-weighted mean";{.quantQ.knn.kNNMeanFunction[0 0f;(1 0f;0 2f);1 3f;1 1f;1f]};1.0948517463551337;1e-12];
.qt.approx["kNNMeanFunction of a constant is that constant";{.quantQ.knn.kNNMeanFunction[0.3 0.1;(1 0f;0 2f;3 3f);3#7f;1 1f;0.5]};7f;1e-12];
.qt.approx["errorLocal is 0 when the weights are all 1";{.quantQ.knn.errorLocal[0 0f;(1 0f;0 2f);1 3f;1 1f;1f]};0f;1e-12];
.qt.approx["aFunc";{.quantQ.knn.aFunc[0 0f;1 0f;0 2f;2f;1 1f]};0.0820849986238988;1e-12];
.qt.approx["uFunc";{.quantQ.knn.uFunc[0 0f;1 0f;0 2f;1 2f]};1 8f;1e-12];
.qt.eq["gradFunc returns one component per feature";{count .quantQ.knn.gradFunc[0 0f;(1 0f;0 2f;1 1f);1 3 2f;1 1f;1f]};2];

/ feature selection: runs nIterations steps and keeps one weight per feature
knnFS:`coordFull`valueFull`learningRatio`decay`kNN`beta`fractionS`nIterations!(
    (0 0f;1 0f;0 1f;1 1f;0.5 0.2;0.2 0.7;0.9 0.4;0.3 0.3);
    0 1 0 1 0.5 0.2 0.9 0.3;
    0.01; 0.9; 3; 1f; 0.5; 2);
.qt.eq["kNNFeatureSelection: one state per iteration plus the initial one";{count .quantQ.knn.kNNFeatureSelection knnFS};3];
.qt.eq["kNNFeatureSelection: weight vector has one entry per feature";{count last exec weight from .quantQ.knn.kNNFeatureSelection knnFS};2];

.qt.finish[];
