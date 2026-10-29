/ Unit tests for lib/quantQ_svm.q -- soft-margin SVM with SGD and classification statistics (.quantQ.svm)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`svm;

/ soft margin of a single observation: lambda*|w|^2 + max(0, 1 - y(w.x - b))
.qt.approx["softMargin: correctly classified outside the margin";{.quantQ.svm.softMargin[1f;1 1f;0f;1 1f;0f]};0f;1e-12];
.qt.approx["softMargin: hinge loss plus penalty";{.quantQ.svm.softMargin[1f;0 0f;0f;1 1f;0.5]};2f;1e-12];
.qt.approx["softMargin: misclassified point";{.quantQ.svm.softMargin[-1f;1 1f;0f;1 1f;0f]};3f;1e-12];
.qt.approx["softMarginGradient inside the margin";{.quantQ.svm.softMarginGradient[1f;0 0f;0f;1 1f;0.5]};1 1 1f;1e-12];
.qt.approx["softMarginGradient outside the margin is the penalty gradient only";{.quantQ.svm.softMarginGradient[1f;2 2f;0f;1 1f;0.5]};0 1 1f;1e-12];
.qt.approx["softMarginGradient matches a finite difference in w";{h:1e-6; f:{.quantQ.svm.softMargin[1f;0.3 -0.2;0.1;x;0.2]}; ((f (0.5+h),0f)-f (0.5-h),0f)%2*h};.quantQ.svm.softMarginGradient[1f;0.3 -0.2;0.1;0.5 0f;0.2]1;1e-6];

/ synthetic, linearly separable data: y = 1 if x1 + x2 > 1 (points near the boundary removed)
svmGrid:(0.05*til 21) cross 0.05*til 21;
svmGrid:svmGrid where 0.15<abs -1+sum each svmGrid;
svmTab:([] yVar:?[1<sum each svmGrid;1f;-1f]; xVar:svmGrid);

/ evaluation with known parameters (the true separating line w=(1,1), b=1)
.qt.approx["evalModel: raw score w.x - b";{exec yVarPredRaw from .quantQ.svm.evalModel[`w`b!(1 1f;1f);([] yVar:1 -1f; xVar:(0.8 0.9;0.1 0.2))]};0.7 -0.7;1e-12];
.qt.true["evalModel: the true line classifies the grid perfectly";{d:.quantQ.svm.evalModel[`w`b!(1 1f;1f);svmTab]; all d[`yVar]=d`yVarPred}];

/ classification statistics: TP=3, FN=1, FP=1, TN=2
svmDM:([] yVar:1 1 1 1 -1 -1 -1; yVarPred:1 1 1 -1 1 -1 -1);
svmSt:.qt.fixture["svmSt";{.quantQ.svm.evalStats svmDM}];
.qt.approx["evalStats: confusion counts";{svmSt`truePositive`trueNegative`falsePositive`falseNegative};3 2 1 1;0f];
.qt.approx["evalStats: accuracy";{svmSt`accuracy};5%7;1e-12];
.qt.approx["evalStats: sensitivity";{svmSt`sensitivity};0.75;1e-12];
.qt.approx["evalStats: specificity";{svmSt`specificity};2%3;1e-12];
.qt.approx["evalStats: precision";{svmSt`precision};0.75;1e-12];
.qt.approx["evalStats: negative predictive value";{svmSt`negPredictiveValue};2%3;1e-12];
.qt.approx["evalStats: miss rate";{svmSt`missRate};0.25;1e-12];
.qt.approx["evalStats: fall-out";{svmSt`fallout};1%3;1e-12];
.qt.approx["evalStats: false discovery rate = FP/(FP+TP)";{svmSt`falseDiscoveryRate};0.25;1e-12];
.qt.approx["evalStats: false omission rate = FN/(FN+TN)";{svmSt`falseOmissionRate};1%3;1e-12];
.qt.approx["evalStats: prevalence threshold";{svmSt`prevalenceThreshold};0.4;1e-12];
.qt.approx["evalStats: threat score";{svmSt`threatScore};0.6;1e-12];
.qt.approx["evalStats: balanced accuracy";{svmSt`balancedAccuracy};0.5*0.75+2%3;1e-12];
.qt.approx["evalStats: F1";{svmSt`F1};0.75;1e-12];
.qt.approx["evalStats: Matthews correlation";{svmSt`MatthewsCorrelation};0.4166666666666667;1e-12];
.qt.approx["evalStats: Fowlkes-Mallows";{svmSt`FowlkesMallowsIndex};0.75;1e-12];
.qt.approx["evalStats: informedness";{svmSt`informedness};0.75+(2%3)-1;1e-12];
.qt.approx["evalStats: markedness";{svmSt`markedness};0.75+(2%3)-1;1e-12];

/ fitting with SGD
svmFitPar:`maxIters`softMarginDecayError`learningRateParams!(20000;-0w;`method`learningRate!(`const;0.05));
.qt.true["fitSoftMargin separates the training data";{f:.quantQ.svm.fitSoftMargin[svmTab;svmFitPar]; d:.quantQ.svm.evalModel[f;svmTab]; 0.95<=avg d[`yVar]=d`yVarPred}];
.qt.eq["fitSoftMargin stops at maxIters when convergence checks are off";{(.quantQ.svm.fitSoftMargin[svmTab;@[svmFitPar;`maxIters;:;200]])`stopReason};`nIterExceeded];
.qt.eq["cvFit returns one soft margin per fold";{count .quantQ.svm.cvFit[svmTab;@[svmFitPar;`maxIters;:;500];4;0.01]};4];

.qt.finish[];
