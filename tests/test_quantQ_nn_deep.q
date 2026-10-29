/ Unit tests for lib/quantQ_nn_deep.q -- deep networks with batches, dropout and annealing (.quantQ.nn)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`nn_deep;

/ batching helpers
.qt.eq["getBatch: inclusive index range";{.quantQ.nn.getBatch[10+til 10;2 4]};12 13 14];
.qt.eq["getBatch: single index";{.quantQ.nn.getBatch[10+til 10;3 3]};enlist 13];
.qt.eq["setBatchIndices: first/last index of each batch";{.quantQ.nn.setBatchIndices[10;2]};(0 4;5 9)];
.qt.eq["setBatchIndices: batches cover the data";{b:.quantQ.nn.setBatchIndices[9;3]; raze {x[0]+til 1+x[1]-x[0]} each b};til 9];
.qt.eq["getEpochBatch[6;6] is epoch 1, batch 0";{.quantQ.nn.getEpochBatch[6;6]};`epoch`batch!1 0];
.qt.eq["getEpochBatch[4;11] is epoch 2, batch 3";{.quantQ.nn.getEpochBatch[4;11]};`epoch`batch!2 3];

/ dropout mask: one vector per layer feeding forward, bias entries always kept
.qt.eq["getDropoutVector: one mask per input/hidden layer";{count .quantQ.nn.getDropoutVector[3 5 4 1;0.8;0.5]};3];
.qt.eq["getDropoutVector: mask lengths follow the structure";{count each .quantQ.nn.getDropoutVector[3 5 4 1;0.8;0.5]};3 5 4];
.qt.true["getDropoutVector: bias entries are 1";{all 1f=first each .quantQ.nn.getDropoutVector[3 5 4 1;0.8;0.5]}];
.qt.true["getDropoutVector: kept units are rescaled by 1/p";{m:.quantQ.nn.getDropoutVector[3 50 1;0.8;0.5]; (all (1_m 0) in 0 1.25) and all (1_m 1) in 0 2f}];
.qt.eq["getDropoutVector with p=1 keeps everything";{distinct raze .quantQ.nn.getDropoutVector[4 6 1;1f;1f]};enlist 1f];

/ learning-rate schedules
.qt.approx["learningFuncs: default is a constant 0.01";{.quantQ.nn.learningFuncs[()!();7]};0.01;1e-12];
.qt.approx["learningFuncs: constant";{.quantQ.nn.learningFuncs[`method`learningRate!(`const;0.2);100]};0.2;1e-12];
.qt.approx["learningFuncs: power law lr*c^p";{.quantQ.nn.learningFuncs[`method`learningRate`power!(`powerLaw;0.1;-0.5);4]};0.05;1e-12];
.qt.approx["learningFuncs: power law treats counter 0 as 1";{.quantQ.nn.learningFuncs[`method`learningRate`power!(`powerLaw;0.1;-0.5);0]};0.1;1e-12];
.qt.approx["learningFuncs: stepwise decay";{.quantQ.nn.learningFuncs[`method`learningRate`frac`step!(`stepwise;0.1;0.5;10);25]};0.025;1e-12];
.qt.approx["learningFuncs: power law with annealing";{.quantQ.nn.learningFuncs[`method`learningRate`power`lengthAnneal`powerAnneal!(`powerLawAnneal;0.1;-0.5;10;-1f);25]};0.004;1e-12];

/ training: two hidden layers on a separable problem
/ interleaved so that each of the two batches contains both classes
dnX:enlist each 0 0.8 0.1 0.9 0.2 0.7 0.3 1f;
dnY:enlist each 0 1 0 1 0 1 0 1f;
dnBucket:`nNeurons`func`nEpochs`nBatches`monitorConvergence`pKeepInputFunc`pKeepHiddenFunc`learningRateFunc!(4 3;`classifier;1500;2;1;{1f};{1f};{0.5});
.qt.true["modelNNDeep: error decreases";{m:.quantQ.nn.modelNNDeep[dnX;dnY;dnBucket]; e:m`Error; (avg -50#e)<0.5*avg 1_51#e}];
.qt.eq["modelNNDeep: one parameter matrix per layer";{count last[.quantQ.nn.modelNNDeep[dnX;dnY;@[dnBucket;`nEpochs;:;5]]]`parsLayers};3];
.qt.eq["modelNNDeep: layer shapes";{p:last[.quantQ.nn.modelNNDeep[dnX;dnY;@[dnBucket;`nEpochs;:;5]]]`parsLayers; {(count x;count first x)} each p};(2 4;5 3;4 1)];
.qt.eq["predictNNDeep returns one prediction per input";{m:last .quantQ.nn.modelNNDeep[dnX;dnY;dnBucket]; count .quantQ.nn.predictNNDeep[dnX;m]};8];
.qt.approx["predictNNDeep: separable data are classified correctly";{m:last .quantQ.nn.modelNNDeep[dnX;dnY;dnBucket]; .quantQ.nn.probToClasses .quantQ.nn.predictNNDeep[dnX;m]};0 1 0 1 0 1 0 1;0f];

.qt.finish[];
