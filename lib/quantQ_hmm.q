/ A Hidden Markov Model is defined with the triplet (a,b,pi) where
/ 'a' is the transition probability matrix (size NxN). a[i;j] is the transition probability from previous state i to current state j.
/ 'b' is the emission proability matrix (size NxT). b[j;obs[t]] is probability of observing symbol obs[t] given the current state j.
/ 'pi' is the initial probability. pi[i] is the initial probability of starting at i.


/ calculate the probability that a model generated a sequence of observations.

.quantQ.hmm.probObservations:{[a;b;pi;obs]
 T:count obs;
 fwd:.quantQ.hmm.forward[a;b;pi;obs];
 :sum fwd[T-1;]
 }

/ initial probability Pr(o1,q=Si|a,b,pi)
.quantQ.hmm.init:{[a;b;pi;obs]
 T:count obs;
 N:count a;
 delta:(T;N)#0f;
 i:0;
 while[i < count a;
  delta[0;i]:pi[i]*b[i;obs[0]];
  i:i+1;
  ];
 :delta}


/ compute alpha[t;j] = Pr(o1,o2,...,ot,qt=j | a,b,pi) where qt=j means the state at time t is j.

.quantQ.hmm.forward:{[a;b;pi;obs]
 alpha:.quantQ.hmm.init[a;b;pi;obs];
 t:1;
 T:count obs;
 N:count a;
 while[t<T;
  j:0;
  while[j<N;
   alpha[t;j]:sum alpha[t-1;]*a[;j]*b[j;obs[t]];
   j:j+1;
   ];
  t:t+1;
  ];
 :alpha
 }


/ returns the likelihood and the index of the most likely hidden state sequence that generates the given observations obs.

.quantQ.hmm.viterbi:{[a;b;pi;obs]
 / obs: sequence of T observations.

 T:count obs;
 N:count a;
 delta:(T;N)#0f;
 i:0;
 / 1. initialisation
 while[i<N;
  delta[0;i]:pi[i]*b[i;obs[0]];
  i:i+1;
  ];
 / 2. recursion
 t:1;
 while[t<T;
  j:0;
  while[j<N;
   delta[t;j]:(max delta[t-1;]*a[;j])*b[j;obs[t]];
   j:j+1;
   ];
  t:t+1;
  ];
 / 3. termination
 pstar:max last delta;
 qtstar:{[t;delta]
  d:max@/:delta;
  delta[t;]?d[t]}[;delta] each til T;
 (pstar;qtstar)}

/ randomly select a state according to the weights.

.quantQ.hmm.weightedRandomDistribution:{[states;weights]
 states[binr[sums weights;first 1?1f]]
 }

/ generate T observations given a HMM defined by a,b,pi,states,V.
/ a: transition matrix
/ b: emission probability matrix
/ pi: initial probability
/ states: states of the HMM
/ V: vocabulary of the observed states.

.quantQ.hmm.generateObservations:{[a;b;pi;states;V;T]
 t:0;
 qt:states?.quantQ.hmm.weightedRandomDistribution[states;pi];
 obs:();
 while[t<T;
  obs:obs,.quantQ.hmm.weightedRandomDistribution[V;b[qt;]]; 
  qt:states?.quantQ.hmm.weightedRandomDistribution[states;a[qt;]];
  t:t+1;
  ];
 obs}
