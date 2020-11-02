\l quantQ_hmm.q

states:`Healthy`Fever;
V:`cold`dizzy`normal;

startProb:(0.6f;0.4f);
transitionProb:((0.7f;0.3f);(0.4f;0.6f));
emissionProb:((0.5f;0.4f;0.1f);(0.1f;0.3f;0.6f));
dict:V!til count V;

/ 1) What is the probability that the model generated the sequence of observations obs?

obs:`normal`cold`dizzy`normal`dizzy;
0N!.quantQ.hmm.probObservations[transitionProb;emissionProb;startProb;dict obs];

/ alternatively it can also be computed using the backward function

bwd:.quantQ.hmm.backward[transitionProb;emissionProb;startProb;dict obs];
0N!sum startProb*emissionProb[;dict obs[0]]*bwd[0;];

/ 2) What sequence of states best explains a sequence of observations?

obs:`normal`cold`dizzy`normal`dizzy;
0N!states[.quantQ.hmm.viterbi[transitionProb;emissionProb;startProb;dict obs][1]];
/ How likely did we get this sequence of states?
0N!.quantQ.hmm.viterbi[transitionProb;emissionProb;startProb;dict obs][0];
