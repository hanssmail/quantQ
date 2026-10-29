/ Unit tests for lib/quantQ_hmm.q -- hidden Markov models (.quantQ.hmm)
/ Model from exa/quantQ_hmm_example.q; reference values from forward/backward recursions and
/ brute-force enumeration of all state paths in Python.
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`hmm;

hmmStates:`Healthy`Fever;
hmmV:`cold`dizzy`normal;
hmmPi:0.6 0.4;
hmmA:(0.7 0.3;0.4 0.6);
hmmB:(0.5 0.4 0.1;0.1 0.3 0.6);
hmmIdx:hmmV!til 3;
hmmObs3:hmmIdx `normal`cold`dizzy;
hmmObs5:hmmIdx `normal`cold`dizzy`normal`dizzy;

/ init / forward / backward
.qt.approx["init: delta[0;i] = pi[i]*b[i;o0]";{first .quantQ.hmm.init[hmmA;hmmB;hmmPi;hmmObs3]};0.06 0.24;1e-12];
.qt.approx["forward: full alpha matrix";{.quantQ.hmm.forward[hmmA;hmmB;hmmPi;hmmObs3]};(0.06 0.24;0.069 0.0162;0.021912 0.009126);1e-12];
.qt.approx["backward: full beta matrix";{.quantQ.hmm.backward[hmmA;hmmB;hmmPi;hmmObs3]};(0.1397 0.0944;0.37 0.34;1 1f);1e-12];

/ likelihood of an observation sequence
.qt.approx["probObservations, 3 observations";{.quantQ.hmm.probObservations[hmmA;hmmB;hmmPi;hmmObs3]};0.031038;1e-12];
.qt.approx["probObservations, 5 observations";{.quantQ.hmm.probObservations[hmmA;hmmB;hmmPi;hmmObs5]};0.0031606224;1e-12];
.qt.approx["forward and backward give the same likelihood";{bwd:.quantQ.hmm.backward[hmmA;hmmB;hmmPi;hmmObs5]; sum hmmPi*hmmB[;hmmObs5 0]*bwd[0]};0.0031606224;1e-12];
.qt.approx["likelihoods of all length-2 sequences sum to 1";{sum {.quantQ.hmm.probObservations[hmmA;hmmB;hmmPi;x]} each (til 3) cross til 3};1f;1e-12];

/ Viterbi: probability and the most likely state path
.qt.approx["viterbi: probability of the best path";{first .quantQ.hmm.viterbi[hmmA;hmmB;hmmPi;hmmObs3]};0.01344;1e-12];
.qt.eq["viterbi: best path (Fever, Healthy, Healthy)";{last .quantQ.hmm.viterbi[hmmA;hmmB;hmmPi;hmmObs3]};1 0 0];
.qt.eq["viterbi: best path for the 5-step example";{last .quantQ.hmm.viterbi[hmmA;hmmB;hmmPi;hmmObs5]};1 0 0 1 1];
/ here the per-step argmax of delta differs from the backtracked optimum (brute force: 0 0 1 1 0 0)
.qt.approx["viterbi: probability, cold cold normal normal dizzy cold";{first .quantQ.hmm.viterbi[hmmA;hmmB;hmmPi;hmmIdx `cold`cold`normal`normal`dizzy`cold]};0.000381024;1e-12];
.qt.eq["viterbi: path is backtracked, not a per-step argmax";{last .quantQ.hmm.viterbi[hmmA;hmmB;hmmPi;hmmIdx `cold`cold`normal`normal`dizzy`cold]};0 0 1 1 0 0];

/ sampling
.qt.eq["weightedRandomDistribution with a degenerate distribution";{.quantQ.hmm.weightedRandomDistribution[hmmStates;1 0f]};`Healthy];
.qt.true["weightedRandomDistribution always returns one of the states";{all (.quantQ.hmm.weightedRandomDistribution[hmmStates;] each 100#enlist 0.5 0.5) in hmmStates}];
.qt.eq["generateObservations returns T observations";{count .quantQ.hmm.generateObservations[hmmA;hmmB;hmmPi;hmmStates;hmmV;25]};25];
.qt.true["generateObservations uses the vocabulary";{all .quantQ.hmm.generateObservations[hmmA;hmmB;hmmPi;hmmStates;hmmV;25] in hmmV}];
.qt.eq["generateObservations with deterministic emissions";{distinct .quantQ.hmm.generateObservations[(1 0f;0 1f);(0 0 1f;1 0 0f);1 0f;hmmStates;hmmV;10]};enlist `normal];

.qt.finish[];
