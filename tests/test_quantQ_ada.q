/ Unit tests for lib/quantQ_ada.q -- AdaBoost with decision stumps (.quantQ.ada)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`ada;

/ fixtures: one feature, y in {-1,1}; the point at f1=4 is "mislabelled" so no stump is perfect
adaTab:([] y:-1 -1 1 -1 1; weights:5#0.2; f1:1 2 3 4 5f);
adaEmpty:([] i:`long$(); alpha:`float$(); stump:());

/ stump
.qt.eq["stump returns +1 above the threshold, -1 otherwise";{.quantQ.ada.stump[1 5 3 7f;4f]};-1 1 -1 1];
.qt.eq["stump: equality counts as below";{.quantQ.ada.stump[1 2 3f;2f]};-1 -1 1];

/ addNegFeatureTable
adaNegTab:([] y:1 -1; weights:0.5 0.5; f1:1 2f; f2:3 4f);
.qt.eq["addNegFeatureTable appends a negated copy of each feature";{cols .quantQ.ada.addNegFeatureTable adaNegTab};`y`weights`f1`f2`f1n`f2n];
.qt.eq["addNegFeatureTable: negated values";{exec f2n from .quantQ.ada.addNegFeatureTable adaNegTab};-3 -4f];

/ weightedError
adaWE1:([] y:1 -1 1 -1; weights:0.1 0.2 0.3 0.4);
adaWE2:([] y:1 -1; weights:0.5 0.5);
.qt.approx["weightedError sums the weights of misclassified points";{.quantQ.ada.weightedError[adaWE1;1 1 -1 -1]};0.5;1e-12];
.qt.approx["weightedError is 0 for a perfect prediction";{.quantQ.ada.weightedError[adaWE2;1 -1]};0f;1e-12];

/ runOneAdaBoost: best stump has weighted error 0.2, so alpha = 0.5*log(0.8/0.2) = log 2
.qt.approx["runOneAdaBoost: alpha of the best stump";{first exec alpha from (.quantQ.ada.runOneAdaBoost `tab`ada!(adaTab;adaEmpty))`ada};log 2f;1e-12];
.qt.approx["runOneAdaBoost: error of the best stump";{last first exec stump from (.quantQ.ada.runOneAdaBoost `tab`ada!(adaTab;adaEmpty))`ada};0.2;1e-12];
.qt.approx["runOneAdaBoost: updated weights still sum to 1";{sum exec weights from (.quantQ.ada.runOneAdaBoost `tab`ada!(adaTab;adaEmpty))`tab};1f;1e-12];
.qt.approx["runOneAdaBoost: misclassified point gets weight 0.5, others 0.125";{asc exec weights from (.quantQ.ada.runOneAdaBoost `tab`ada!(adaTab;adaEmpty))`tab};0.125 0.125 0.125 0.125 0.5;1e-12];
.qt.eq["runOneAdaBoost appends one row to the ada table";{count (.quantQ.ada.runOneAdaBoost `tab`ada!(adaTab;adaEmpty))`ada};1];

/ adaBoostPrediction with a single known stump (feature f1, split 2)
adaOne:([] i:enlist 1; alpha:enlist log 2f; stump:enlist (`f1;2f;0.2));
.qt.approx["adaBoostPrediction: predictions from a single stump";{exec prediction from .quantQ.ada.adaBoostPrediction `tab`ada!(adaTab;adaOne)};-1 -1 1 1 1f;0f];
.qt.approx["adaBoostPrediction: isPredicted = prediction*y";{exec isPredicted from .quantQ.ada.adaBoostPrediction `tab`ada!(adaTab;adaOne)};1 1 1 -1 1f;0f];

/ adaPredictionStats -- standard confusion-matrix definitions:
/   TP: y=1,pred=1   TN: y=-1,pred=-1   FP: y=-1,pred=1   FN: y=1,pred=-1
adaPred:([] y:1 1 1 -1 -1 -1; prediction:1 -1 -1 1 -1 -1);
.qt.approx["adaPredictionStats: true positives";{.quantQ.ada.adaPredictionStats[adaPred]`truePositive};1;0f];
.qt.approx["adaPredictionStats: true negatives";{.quantQ.ada.adaPredictionStats[adaPred]`trueNegative};2;0f];
.qt.approx["adaPredictionStats: false positives (y=-1, predicted 1)";{.quantQ.ada.adaPredictionStats[adaPred]`falsePositive};1;0f];
.qt.approx["adaPredictionStats: false negatives (y=1, predicted -1)";{.quantQ.ada.adaPredictionStats[adaPred]`falseNegative};2;0f];
.qt.approx["adaPredictionStats: accuracy";{.quantQ.ada.adaPredictionStats[adaPred]`accuracy};0.5;1e-12];
.qt.approx["adaPredictionStats: recall = TP/(TP+FN)";{.quantQ.ada.adaPredictionStats[adaPred]`recall};1%3;1e-12];
.qt.approx["adaPredictionStats: precision = TP/(TP+FP)";{.quantQ.ada.adaPredictionStats[adaPred]`precision};0.5;1e-12];
.qt.approx["adaPredictionStats: specificity = TN/(TN+FP)";{.quantQ.ada.adaPredictionStats[adaPred]`specificity};2%3;1e-12];
.qt.approx["adaPredictionStats: F1";{.quantQ.ada.adaPredictionStats[adaPred]`F1score};0.4;1e-12];

.qt.finish[];
