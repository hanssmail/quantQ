// returns the likelihood and the index of the most likely hidden state sequence that generates the given observations obs.

.quantQ.hmm.viterbi:{[a;b;pi;obs]
    // a -- transition probability matrix (size NxN). a[i;j] is the transition probability from previous state i to current state j.
    // b -- emission probability matrix (size NxT). b[j;obs[t]] is the state observation likelihood of the observation symbol obs[t] given the current state j.
    // pi -- initial probability. pi[i] is the initial probability of starting at i.
    // obs -- sequence of T observations.

    T:count obs;
    N:count a;
    delta:(T;N)#0f;
    i:0;
    // 1. initialisation
    while[i<N;
        delta[0;i]:pi[i]*b[i;obs[0]];
        i:i+1;
    ];
    // 2. recursion
    t:1;
    while[t<T;
        j:0;
        while[j<N;
            delta[t;j]:(max delta[t-1;]*a[;j])*b[j;obs[t]];
            j:j+1;
        ];
        t:t+1;
    ];
    // 3. termination
    pstar:max last delta;
    qtstar:{[t;delta]
        d:max@/:delta;
        delta[t;]?d[t]
    }[;delta] each til T;
    :(pstar;qtstar);
 };
 
