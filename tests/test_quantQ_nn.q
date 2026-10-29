/ Unit tests for lib/quantQ_nn.q -- single-hidden-layer neural network (.quantQ.nn)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`nn;

/ encoders
.qt.eq["encoderIndex";{.quantQ.nn.encoderIndex `a`b`c};`a`b`c!0 1 2];
.qt.eq["encoderMulti one-hot encodes";{.quantQ.nn.encoderMulti `a`b};`a`b!(1 0f;0 1f)];
.qt.approx["probToClasses: argmax for multi-class output";{.quantQ.nn.probToClasses (0.2 0.7 0.1;0.6 0.3 0.1;0.1 0.1 0.8)};1 0 2;0f];
.qt.approx["probToClasses: threshold 0.5 for a single output";{.quantQ.nn.probToClasses (enlist 0.7;enlist 0.2;enlist 0.5)};1 0 1;0f];
.qt.approx["accuracy";{.quantQ.nn.accuracy[1 0 1;1 1 1]};2%3;1e-12];

/ activation and error functions
.qt.approx["sigmoid(0) = 1/2";{.quantQ.nn.sigmoid 0f};0.5;1e-12];
.qt.approx["sigmoid(x) + sigmoid(-x) = 1";{.quantQ.nn.sigmoid[1.7]+.quantQ.nn.sigmoid -1.7};1f;1e-12];
.qt.approx["sigmoidErr is the binary cross-entropy";{.quantQ.nn.sigmoidErr[(enlist 1f;enlist 0f);(enlist 0.8;enlist 0.3)]};neg (log 0.8)+log 0.7;1e-12];
.qt.approx["softmax normalises each row";{.quantQ.nn.softmax (1 2 3f;0 0 0f)};((exp 1 2 3f)%sum exp 1 2 3f;3#1%3);1e-12];
.qt.approx["softmaxErr is the categorical cross-entropy";{.quantQ.nn.softmaxErr[(0 1 0f;1 0 0f);(0.2 0.5 0.3;0.6 0.2 0.2)]};neg (log 0.5)+log 0.6;1e-12];
.qt.eq["linearNN is the identity";{.quantQ.nn.linearNN 1 2 3f};1 2 3f];
.qt.approx["linearNNErr is the sum of squared errors";{.quantQ.nn.linearNNErr[(1 2f;3 4f);(1 1f;3 2f)]};5f;1e-12];
.qt.eq["funcNN dispatch";{key .quantQ.nn.funcNN};`classifier`multiClassifier`nonlinearReg];
.qt.eq["funcErrNN dispatch";{key .quantQ.nn.funcErrNN};`classifier`multiClassifier`nonlinearReg];

/ weight initialisation: argIn rows, argOut columns, each column centred
.qt.eq["weightNNInit shape";{w:.quantQ.nn.weightNNInit[3;2]; (count w;count first w)};3 2];
.qt.approx["weightNNInit: each column has zero mean";{avg .quantQ.nn.weightNNInit[5;3]};3#0f;1e-12];

/ training on a separable 1-d classification problem
nnX:enlist each 0 0.1 0.2 0.8 0.9 1f;
nnY:enlist each 0 0 0 1 1 1f;
.qt.true["modelNN: training error decreases";{m:.quantQ.nn.modelNN[nnX;nnY;3;`count;3000;1.0;1;`classifier]; (last m`Error)<0.5*(m`Error)1}];
.qt.approx["modelNN: separable data are classified correctly";{m:.quantQ.nn.modelNN[nnX;nnY;3;`count;3000;1.0;0;`classifier]; .quantQ.nn.probToClasses .quantQ.nn.predict[nnX;`classifier;m]};0 0 0 1 1 1;0f];
.qt.eq["modelNN: counter equals the number of iterations";{(.quantQ.nn.modelNN[nnX;nnY;3;`count;50;0.5;0;`classifier])`Counter};50];
.qt.eq["predict accepts the monitored (table) output";{m:.quantQ.nn.modelNN[nnX;nnY;3;`count;10;0.5;1;`classifier]; count .quantQ.nn.predict[nnX;`classifier;m]};6];
.qt.true["modelNN regression: fit improves on y=2x";{rx:enlist each 0.1*til 10; ry:2*rx; m:.quantQ.nn.modelNN[rx;ry;4;`count;2000;0.05;1;`nonlinearReg]; (last m`Error)<0.1*(m`Error)1}];

.qt.finish[];
