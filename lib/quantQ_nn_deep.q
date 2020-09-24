// get Batch
.quantQ.nn.getBatch:{[data;y] 
    // data -- array to be batched
    // y -- pair of indices, init-end (both inclusive)
    indices:first[y]+til 1+last[y]-first[y];
    $[1<count indices;
        outt:data[indices];
        outt:data[indices]
    ];
    : outt;
 };

// set Batch indices
.quantQ.nn.setBatchIndices:{[n;nBatches] 
    // n -- size of data set
    // nBatches
    // return the beggining/end of the batch
    :{(first[x];last[x])} each (nBatches;0N)#til n;
 };

// get epoch/batch from nRun
// .quantQ.nn.getEpochBatch[6;6] epoch=1, batch=0
.quantQ.nn.getEpochBatch:{[nBatches;nRun] 
    // nRun -- run number
    // nBatches -- number of batches 
    // counting starts from 0
    epochN: floor nRun%nBatches;
    batchN:mod[nRun;nBatches];
    // return dictionary
    :(`epoch`batch)!(epochN;batchN);
 };


// deep neural network with batches
.quantQ.nn.modelNNDeep:{[input;output;bucket]
    // input -- array of input values
    // output -- array of output values
    // bucket -- parameters
    // determine epochs/batches; default 100 hundred iteration with all data at once
    bucket:((`nEpochs`nBatches`typeTerminal`relError)!(100;1;`count;1)),bucket;
    nRuns: bucket[`nEpochs]*bucket[`nBatches];
    // determine batch indices
    batchIndices:.quantQ.nn.setBatchIndices[count[input];bucket[`nBatches]]; 
    bucket[`batchIndices]:batchIndices;    
    // initial data manipulation
    input:1.0,'input;
    // counting dimension of the input including intercept
    nInput: count flip input;
    // counting dimension of the output
    nOutput: count first output;
    // initialization of layer: Input -> Hidden
    layerI2H:.quantQ.nn.weightNNInit[nInput;first bucket[`nNeurons]];
    // initialization of deep network if present
    layerDeep:();
    if[1<cntN:count bucket[`nNeurons];
        layerDeep:{[nNeurons;k] .quantQ.nn.weightNNInit[1+nNeurons[k];nNeurons[k+1]]}[bucket[`nNeurons];] each til cntN-1;  
    ];
    // initialization of layer: Hidden (with bias neuron) -> Output
    layerH2O:.quantQ.nn.weightNNInit[1 + last[bucket[`nNeurons]];nOutput];
    // all layers initialised
    layers:enlist[layerI2H],layerDeep,enlist[layerH2O];
    // decide between two specified terminal conditions
    if[bucket[`typeTerminal]=`count;terminus:nRuns];
    if[bucket[`typeTerminal]=`relativeError;terminus:{[argTerminal;param] param[`Error]>argTerminal}[bucket[`relError];]];
    // monitoring convergence
    bucket:(enlist[`monitorConvergence]!(enlist[0])),bucket;
    // training the network and return output (depends on monitoring the convergence):
    $[bucket[`monitorConvergence]=1;
        :(.quantQ.nn.forwardBackwardNNDeepwErrorDrop[input;output]\)[terminus;`out`parsLayers`Error`bucket`Counter!(0,();layers;1f;bucket;0)];
        :(.quantQ.nn.forwardBackwardNNDeepwErrorDrop[input;output]/)[terminus;`out`parsLayers`Error`bucket`Counter!(0,();layers;1f;bucket;0)]
    ];
 };

// dropout function to generate the dropout mask
.quantQ.nn.getDropoutVector:{[nnStructure;pKeepInput;pKeepHidden]
    // nnStructure -- structure of the NN
    // pKeepInput -- probability to keep input
    // pKeepHidden -- probability to keep hidden layer
    vInp:enlist {[p;x] "f"$($[x=0;1.0;(1%p)*(p > first 1?1.0)])  }[pKeepInput;] each til nnStructure[0];
    vHid:{[p;x] {[p;x] "f"$($[x=0;1.0;(1%p)*(p > first 1?1.0)])  }[p;] each til x }[pKeepHidden;]  each {-1_1_x} nnStructure;
    // vDropout: raze (vInp;vHid);
    :raze (vInp;vHid);
 };

// one step
.quantQ.nn.forwardBackwardNNDeepwErrorDrop:{[inputValue;outputValue;params]
    // inputValue -- array of input values
    // outputValue -- array of output values
    // params -- bucket of parameters to be updated
    // initialisation
    func:params[`bucket][`func];
    // infer counter for epochs -- used in the 
    epochN:{x[`epoch]} nEpBa:.quantQ.nn.getEpochBatch[params[`bucket][`nBatches];params[`Counter]]; 
    // learning rate; default    
    learningRate: $[null params[`bucket][`learningRateFunc];0.015;params[`bucket][`learningRateFunc][epochN]];
    // pKeepInputFunc -- dropout annealing; default
    pKeepInput:$[null params[`bucket][`pKeepInputFunc];0.8;params[`bucket][`pKeepInputFunc][epochN]];
    // pKeepHiddenFunc -- dropout annealing; default
    pKeepHidden:$[null params[`bucket][`pKeepHiddenFunc];0.5;params[`bucket][`pKeepHiddenFunc][epochN]];
    // counting dimension of the input including intercept
    nInput: count flip inputValue;
    // counting dimension of the output
    nOutput: count first outputValue;
    // define dropout mask
    nnStructure:raze (nInput; 1+params[`bucket][`nNeurons]; nOutput);
    vDropout:.quantQ.nn.getDropoutVector[nnStructure;pKeepInput;pKeepHidden];
    // arrays of parameters
    layers:params[`parsLayers];
    cntL:count layers;
    // define the updated hidden nodes
    // layersT1:layers;
    nodesT1:(cntL)#enlist[()];
    // Batches
    inputValue:.quantQ.nn.getBatch[inputValue;params[`bucket][`batchIndices][nEpBa[`batch]]];
    outputValue:.quantQ.nn.getBatch[outputValue;params[`bucket][`batchIndices][nEpBa[`batch]]];
    // Forward stage
    // update input layer
    nodesT1[0]:1.0,/:.quantQ.nn.sigmoid[(vDropout[0]*/:inputValue) mmu params[`parsLayers][0]];
    // update layer by layer deep in the network if present
    c:0;
    while[c<cntL-2;
        nodesT1[c+1]:1.0,/:.quantQ.nn.sigmoid[(vDropout[c+1]*/:nodesT1[c]) mmu params[`parsLayers][c+1]];
        c+:1;
    ];
    // update output layer (model outcome)
    nodesT1[cntL-1]:.quantQ.nn.funcNN[func][(vDropout[cntL-1]*/:nodesT1[cntL-2]) mmu last params[`parsLayers]];
    // Backward stage
    deltaLayerDeep: (cntL)#enlist ();
    // error Hidden -> Output
    // deltaLayerH2O:(outputValue-nodesT1[cntL-1]);
    deltaLayerDeep[cntL-1]:(outputValue-nodesT1[cntL-1]);
    // iterate through hidden layers
    c:0;
    while[c<cntL-1;
        deltaLayerDeep[cntL-(2+c)]:1_/: vDropout[cntL-(1+c)]*/:$[deltaLayerDeep[cntL-(1+c)];flip params[`parsLayers][cntL-(1+c)]]*nodesT1[cntL-(2+c)]*1-nodesT1[cntL-(2+c)];
        c+:1;
    ];
    // update params -- this will be new params[`parsLayers]
    paramsNew:params[`parsLayers];
    c:0;
    while[c<cntL-1;
        paramsNew[c+1]:params[`parsLayers][c+1]+learningRate*
            flip[nodesT1[c]] mmu deltaLayerDeep[c+1];
        c+:1;
    ];
    // input layer
    paramsNew[0]:params[`parsLayers][0]+learningRate*flip[inputValue] mmu 
    deltaLayerDeep[0];
    // error Input -> Hidden
    //deltaLayerI2H:1_/:$[deltaLayerH2O;flip params[`parsHid2Out]]*hiddenLayer*1-hiddenLayer;
    // updating the model
    :`out`parsLayers`Error`bucket`Counter!(nodesT1[cntL-1];paramsNew;.quantQ.nn.funcErrNN[func][outputValue;nodesT1[cntL-1]];params[`bucket];params[`Counter]+1);
 };

// get prediction
.quantQ.nn.predictNNDeep:{[inputData;params]
    // inputData -- input data
    // params -- predicted model
    // func -- output layer function
    layers:params[`parsLayers];
    cntL:count layers;
    // update input layer
    nodesT1[0]:1.0,/:.quantQ.nn.sigmoid[(1.0,/:inputData) mmu params[`parsLayers][0]];
    // update layer by layer deep in the network if present
    c:0;
    while[c<cntL-2;
        nodesT1[c+1]:1.0,/:.quantQ.nn.sigmoid[nodesT1[c] mmu params[`parsLayers][c+1]];
        c+:1;
    ];
    // update output layer (model outcome)
    outputData:.quantQ.nn.funcNN[params[`bucket][`func]][nodesT1[cntL-2] mmu last params[`parsLayers]];
    // prediction
    :outputData;
 };   


// additional functions for annealing
.quantQ.nn.learningFuncs:{[bucket;counter]
    // bucket -- dictionary with parameter
    // counter -- number of iterations passed
    // default, when empty dictionary provided
    bucket:(`method`learningRate!(`const;0.01)),bucket;
    outt:0.0;
    // constant annealing
    if[bucket[`method]=`const;
        outt:bucket[`learningRate]
    ];
    // power-law decay
    if[bucket[`method]=`powerLaw;
        outt:bucket[`learningRate]*xexp[(1|counter);bucket[`power]]
    ];
    // stepwise decay
    if[bucket[`method]=`stepwise;
        outt:bucket[`learningRate]*xexp[bucket[`frac];floor[counter%bucket[`step]]]
    ];
    // power-law with annealing
    if[bucket[`method]=`powerLawAnneal;
        outt:bucket[`learningRate]*xexp[(1|counter);bucket[`power]]*xexp[(1|counter)-bucket[`lengthAnneal]*floor (1|counter)%bucket[`lengthAnneal];bucket[`powerAnneal]]
    ];
   // output
   :outt;
 };


