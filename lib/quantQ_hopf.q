/////////////////////////////////////////////////
// Hopfield networks (another Ising model)
/////////////////////////////////////////////////

// Get a Hopfield network corresponding to a single binary state
.quantQ.hopf.getHopfieldNet:{[state]
    // state -- array of +-1
    out:{x*\:x} state;
    // remove diagonal term
    :out*"f"$not v=/:v:til count state;
 };
 
// Get energy of the Hopfield network 
.quantQ.hopf.getEnergyHopfield:{[hopnet;state]
    // hopnet -- Hopfield network
    // state -- array of +-1
    // energy:
    :0.5*neg wsum[state;] hopnet mmu state;
 };
 
// Get the minimum attainable energy for a given state 
.quantQ.hopf.getMinEnergyHopfield:{[state]
    // state -- array of +-1
    // minimum energy
    :neg sum til count state;        
 }; 
 
// Get an energy change when going from a state to another stateNew
.quantQ.hopf.deltaEnergyHopfield:{[hopnet;state;stateNew]
    // hopnet -- Hopfield network
    // state -- array of +-1
    // stateNew -- array of +-1
    // energy state -> stateNew
    :neg wsum[stateNew-state;] hopnet mmu state;
 };

// Signum function
.quantQ.hopf.signumHopfield:{[state]
    // state -- array
    // signum with >=0 decision rule
    :?[state>=0;1f;-1f];
 };
 
// Perform one transformation of a state
.quantQ.hopf.oneHopfieldJump:{[hopnet;state]
    // hopnet -- Hopfield network
    // state -- array of +-1    
    // transform input based on the Hopfield net
    :.quantQ.hopf.signumHopfield hopnet mmu state;
 }; 

// Get the ground state of the Hopfield network; random search for an attractor with min Energy
.quantQ.hopf.getGroundState:{[hopnet]
    // hopnet -- Hopfield network   
    // dimensionality
    dimH:count hopnet;
    // min energy
    minE: .quantQ.hopf.getMinEnergyHopfield[dimH#1f];
    // random search for the ground state (q lambdas cannot see outer locals, so minE is passed in)
    groundState: {x[`state]}({[minE;bucket]
        // random starting point of +-1
        startingPoint:-1 1f (count bucket[`hopnet])?2;
        // transform
        stateNew: .quantQ.hopf.oneHopfieldJump[bucket[`hopnet];startingPoint];
        // compare energy; give up after maxIter random starts
        bucket[`cnt]:bucket[`cnt]+1;
        if[(minE=.quantQ.hopf.getEnergyHopfield[bucket[`hopnet];stateNew]) or bucket[`cnt]>=bucket[`maxIter];
            bucket[`continue]:0b;
        ];
        // return updated bucket
        :(`hopnet`state`continue`cnt`maxIter)!(bucket[`hopnet];stateNew;bucket[`continue];bucket[`cnt];bucket[`maxIter]);
        }[minE]/)[{x[`continue]};(`hopnet`state`continue`cnt`maxIter)!(hopnet;dimH#1f;1b;0;10000)];
    // return ground state 
    : groundState;
 };   
// Example: .quantQ.hopf.getGroundState[hopnet]
 
// Get a Hopfield network corresponding to a pair of (different) states 
.quantQ.hopf.getHopfieldNet2S:{[state1;state2]
    // state1,state2 -- array of +-1 
    out: state1 *\: state2;
    // remove diagonal term
    :out*"f"$not v=/:v:til count state1;
 };

// exa: state:(1 -1 1 1f); state3:(1 -1 1 -1f);
// .quantQ.hopf.getHopfieldNet[state]
// .quantQ.hopf.getHopfieldNet2S[state;state]
// .quantQ.hopf.getHopfieldNet2S[state;state3]

// One transformation of a state with a bias
.quantQ.hopf.oneHopfieldJumpBias:{[hopnet;state;bias]
    // hopnet -- Hopfield network
    // state -- array of +-1    
    // bias -- array of bias
    // transform input based on the Hopfield net with bias vector
    :.quantQ.hopf.signumHopfield bias+hopnet mmu state;
 }; 

// Calculating Hopf network for stateInput ->stateTarget using Hebbian descent        
.quantQ.hopf.runHebbDescent:{[stateInput;stateTarget;params]
    // stateInput -- array of +-1    
    // stateTarget -- array of +-1      
    // params -- dictionary with params, default is params:()!()
    // initiate the network 
    // W mmu input should give target, so W = target *\: input
    hopnetInit:.quantQ.hopf.getHopfieldNet2S[stateTarget;stateInput];      
    // initiate the bias
    biasInit: count[stateInput]#0.0;
    // default learning rate, stopping rule (L2 by default)
    params:(`learningRate`relError`maxIter!(0.1;0.01;100)),params;
    // Hebbian descent iteration
    bucketOut:({[bucket]
        params:bucket[`params];    
        // calculate target
        stateImplied: .quantQ.hopf.oneHopfieldJumpBias[bucket[`hopnet];bucket[`stateInput];bucket[`bias]];
        // correction for the Hopfield net
        deltaW:params[`learningRate]*
        .quantQ.hopf.getHopfieldNet2S[bucket[`stateTarget];bucket[`stateInput]]-
        .quantQ.hopf.getHopfieldNet2S[stateImplied;bucket[`stateInput]];
        // correction for the bias 
        deltaB:params[`learningRate]*
        bucket[`stateTarget]-stateImplied;
        // update the Hopfield net
        bucket[`hopnet]:bucket[`hopnet]+deltaW;
        // update the bias
        bucket[`bias]:bucket[`bias]+deltaB;
        // update counter
        bucket[`cnt]:bucket[`cnt]+1;
        // error L2
        errorL2:sqrt {wsum[x;x]} stateImplied-bucket[`stateTarget];
        // stopping rule
        if[(bucket[`cnt]>=params[`maxIter]) or (errorL2<=params[`relError]);
            bucket[`continue]:0b
        ];
        // update bucket
        :(`stateInput`stateTarget`continue`cnt`params`hopnet`bias`errorL2)!
        (bucket[`stateInput];bucket[`stateTarget];bucket[`continue];bucket[`cnt];params;bucket[`hopnet];bucket[`bias];errorL2);
    }/)[{x[`continue]};(`stateInput`stateTarget`continue`cnt`params`hopnet`bias`errorL2)!
        (stateInput;stateTarget;1b;0;params;hopnetInit;biasInit;0wf)];
    // return output
    :(`hopnet`bias)!bucketOut[`hopnet`bias];
 };
 
// Calculating Hopf network for stateInput ->stateTarget using Hebbian descent and fixed bias
.quantQ.hopf.runHebbDescentBiased:{[stateInput;stateTarget;bias;params]
    // stateInput -- array of +-1    
    // stateTarget -- array of +-1      
    // bias -- array with fixed bias
    // params -- dictionary with params, default is params:()!()
    // initiate the network 
    // W mmu input should give target, so W = target *\: input
    hopnetInit:.quantQ.hopf.getHopfieldNet2S[stateTarget;stateInput];      
    // default learning rate, stopping rule (L2 by default)
    params:(`learningRate`relError`maxIter!(0.1;0.01;100)),params;
    // Hebbian descent iteration
    bucketOut:({[bucket]
        params:bucket[`params];    
        // calculate target
        stateImplied: .quantQ.hopf.oneHopfieldJumpBias[bucket[`hopnet];bucket[`stateInput];bucket[`bias]];
        // correction for the Hopfield net
        deltaW:params[`learningRate]*
        .quantQ.hopf.getHopfieldNet2S[bucket[`stateTarget];bucket[`stateInput]]-
        .quantQ.hopf.getHopfieldNet2S[stateImplied;bucket[`stateInput]];
        // update the Hopfield net
        bucket[`hopnet]:bucket[`hopnet]+deltaW;
        // update counter
        bucket[`cnt]:bucket[`cnt]+1;
        // error L2
        errorL2:sqrt {wsum[x;x]} stateImplied-bucket[`stateTarget];
        // stopping rule
        if[(bucket[`cnt]>=params[`maxIter]) or (errorL2<=params[`relError]);
            bucket[`continue]:0b
        ];
        // update bucket
        :(`stateInput`stateTarget`continue`cnt`params`hopnet`bias`errorL2)!
        (bucket[`stateInput];bucket[`stateTarget];bucket[`continue];bucket[`cnt];params;bucket[`hopnet];bucket[`bias];errorL2);
    }/)[{x[`continue]};(`stateInput`stateTarget`continue`cnt`params`hopnet`bias`errorL2)!
        (stateInput;stateTarget;1b;0;params;hopnetInit;bias;0wf)];
    // return output
    :(`hopnet`bias)!bucketOut[`hopnet`bias];
 };
 