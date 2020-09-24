.quantQ.nn.encoderIndex:{[arg]
    // arg -- list of classes to be encoded
    :arg!til count arg;
 };

.quantQ.nn.encoderMulti:{[arg]
    // arg -- list of classes to be encoded
    :arg!`float$arg=/:arg;
 };

.quantQ.nn.probToClasses:{[arg]
    // arg -- the probabilities of belonging to a particular class as returned by `multiClassifier
    raze $[1=count arg[0];`int$({x>=0.5}')arg;(til count arg[0])where each(arg=')max each arg]
 };

.quantQ.nn.accuracy:{[true;pred]
    :(sum true=pred)%count true;
 };

.quantQ.nn.sigmoid:{[arg]
    // arg -- argument of the function
    :1%1+exp neg arg;
 };

.quantQ.nn.sigmoidErr:{[argX;argY]
    // argX -- array of true values
    // argY -- array of predicted values
    :neg sum sum flip (argX*log argY) + ((1.0-argX)*log[1.0 - argY]);
 };

.quantQ.nn.softmax:{[arg]
    // arg -- input argument (array with entry for every class)
    :exp[arg]%sum flip exp[arg];
 };

.quantQ.nn.softmaxErr:{[argX;argY]
    // argX -- true value
    // argY -- predicted values
    :neg sum sum flip argX*log[argY];
 };

.quantQ.nn.linearNN:{[arg]
    // arg -- input argument
    :arg;
 };

.quantQ.nn.linearNNErr:{[argX;argY]
    // argX -- true value
    // argY -- predicted values
    :sum sum z*z:(argX-argY);
 };

.quantQ.nn.funcNN: `classifier`multiClassifier`nonlinearReg!
                   (.quantQ.nn.sigmoid;.quantQ.nn.softmax;.quantQ.nn.linearNN);
                   
.quantQ.nn.funcErrNN: `classifier`multiClassifier`nonlinearReg!
                      (.quantQ.nn.sigmoidErr;.quantQ.nn.softmaxErr;.quantQ.nn.linearNNErr);

.quantQ.nn.weightNNInit:{[argIn;argOut]
    // argIn -- number of input arguments
    // argOut -- number of output arguments
    :flip flip[weights]-avg weights:{[x;y]x?1.0}[argOut;]each til argIn;
 };

.quantQ.nn.forwardBackwardNNwError:{[inputValue;outputValue;learningRate;func;params]
    // inputValue -- array of input values
    // outputValue -- array of output values
    // learningRate -- learning rate for update
    // func -- purpose of the NN: `classifier`multiClassifier`nonlinearReg
    // params -- bucket of parameters to be updated
    // Forward stage
    // mapping Input -> Hidden using current parameters
    hiddenLayer:1.0,/:.quantQ.nn.sigmoid[inputValue mmu params[`parsIn2Hid]];
    // mapping Hidden -> Output using current parameters
    outputLayer:.quantQ.nn.funcNN[func][hiddenLayer mmu params[`parsHid2Out]];
    // Backward stage
    // error Hidden -> Output
    deltaLayerH2O:(outputValue-outputLayer);
    // error Input -> Hidden
    deltaLayerI2H:1_/:$[deltaLayerH2O;flip params[`parsHid2Out]]*hiddenLayer*1-hiddenLayer;
    // updating the model
    :`out`parsHid2Out`parsIn2Hid`Error`Counter!(outputLayer;params[`parsHid2Out]+learningRate*
        flip[hiddenLayer] mmu deltaLayerH2O;params[`parsIn2Hid]+learningRate*flip[inputValue] mmu 
        deltaLayerI2H;.quantQ.nn.funcErrNN[func][outputValue;outputLayer];params[`Counter]+1);
 };

.quantQ.nn.modelNN:{[input;output;nNeurons;typeTerminal;argTerminal;learningRate;monitorConvergence;func]
    // input -- array of input values
    // output -- array of output values
    // nNeurons -- number of neurons in hidden layer
    // typeTerminal -- type of terminal condition:`count`relativeError
    // argTerminal -- parameter for terminal condition
    // learningRate -- learning rate for update
    // monitorConvergence -- binary variable to monitor convergence
    // func -- purpose of the NN: `classifier`multiClassifier`nonlinearReg
    // adding bias/intercept to the input data
    input:1.0,'input;
    // counting dimension of the input including intercept
    nInput: count flip input;
    // counting dimension of the output
    nOutput: count first output;
    // initialization of layer: Input -> Hidden
    layerI2H:.quantQ.nn.weightNNInit[nInput;nNeurons];
    // initialization of layer: Hidden (with bias neuron) -> Output
    layerH2O:.quantQ.nn.weightNNInit[1+nNeurons;nOutput];
    // decide between two specified terminal conditions
    $[typeTerminal=`count;terminus:argTerminal;];
        $[typeTerminal=`relativeError;terminus:{[argTerminal;param] param[`Error]>argTerminal}[argTerminal;];
    ];
    // training the network and return output (depends on monitoring the convergence):
    $[monitorConvergence=1;
        :(.quantQ.nn.forwardBackwardNNwError[input;output;learningRate;func]\)[terminus;`out`parsHid2Out`parsIn2Hid`Error`Counter!(0,();layerH2O;layerI2H;1f;0)];
        :(.quantQ.nn.forwardBackwardNNwError[input;output;learningRate;func]/)[terminus;`out`parsHid2Out`parsIn2Hid`Error`Counter!(0,();layerH2O;layerI2H;1f;0)]
    ];
 };

.quantQ.nn.predict:{[input;func;model]
    // input -- array of input values
    // func -- purpose of the NN: `classifier`multiClassifier`nonlinearReg
    // model -- the model, as output by .quantQ.nn.modelNN
    model:$[98h=type model;last model;model];
    input:1f,'input;
    hiddenLayer:1.0,/:.quantQ.nn.sigmoid[input mmu model[`parsIn2Hid]];
    outputLayer:.quantQ.nn.funcNN[func][hiddenLayer mmu model[`parsHid2Out]];
    outputLayer
 };
