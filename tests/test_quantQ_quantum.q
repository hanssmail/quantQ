/ Unit tests for lib/quantQ_quantum.q -- quantum-circuit simulator (.quantQ.quantum)
/ Gate definitions checked against their standard matrices (e.g. Nielsen & Chuang / Qiskit conventions).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`quantum;

/ helpers: a state is a table of complex amplitudes (re, im); a matrix is a list of such rows
qs:{[r;i] ([] re:"f"$r; im:"f"$i)};
qm:{[r;i] {[r;i] ([] re:"f"$r; im:"f"$i)}'[r;i]};
qKet0:qs[1 0;0 0];
qKet1:qs[0 1;0 0];
qH:.quantQ.quantum.H;
qApply:{[g;s] .quantQ.quantum.multMatVec[g;s]};
qNorm2:{[s] sum ((s`re)*s`re)+(s`im)*s`im};
qProbs:{[s] exec probs from .quantQ.quantum.measurement s};
qPi:acos -1;
qR:sqrt 0.5;

/ single-qubit gates
.qt.approx["PauliX |0> = |1>";{qApply[.quantQ.quantum.PauliX;qKet0]};qKet1;1e-12];
.qt.approx["PauliY |0> = i|1>";{qApply[.quantQ.quantum.PauliY;qKet0]};qs[0 0;0 1];1e-12];
.qt.approx["PauliZ |1> = -|1>";{qApply[.quantQ.quantum.PauliZ;qKet1]};qs[0 -1;0 0];1e-12];
.qt.approx["One is the identity";{qApply[.quantQ.quantum.one;qs[0.6 0;0 0.8]]};qs[0.6 0;0 0.8];1e-12];
.qt.approx["Null annihilates";{qApply[.quantQ.quantum.null;qKet0]};qs[0 0;0 0];1e-12];
.qt.approx["H |0> = (|0>+|1>)/sqrt 2";{qApply[qH;qKet0]};qs[qR,qR;0 0];1e-12];
.qt.approx["H H = I";{qApply[qH;qApply[qH;qs[0.6 0;0 0.8]]]};qs[0.6 0;0 0.8];1e-12];
.qt.approx["S |1> = i|1>";{qApply[.quantQ.quantum.S;qKet1]};qs[0 0;0 1];1e-12];
.qt.approx["Sconj S = I";{qApply[.quantQ.quantum.Sconj;qApply[.quantQ.quantum.S;qKet1]]};qKet1;1e-12];
.qt.approx["T |1> = exp(i pi/4)|1>";{qApply[.quantQ.quantum.T;qKet1]};qs[0,qR;0,qR];1e-12];
.qt.approx["T T = S";{qApply[.quantQ.quantum.T;qApply[.quantQ.quantum.T;qKet1]]};qApply[.quantQ.quantum.S;qKet1];1e-12];
.qt.approx["Tconj T = I";{qApply[.quantQ.quantum.Tconj;qApply[.quantQ.quantum.T;qKet1]]};qKet1;1e-12];

/ parametrised gates
.qt.true["Rx works with default parameters";{4=count raze .quantQ.quantum.Rx[()!()]}];
.qt.approx["Rx(pi) |0> = -i|1>";{qApply[.quantQ.quantum.Rx[enlist[`theta]!enlist qPi];qKet0]};qs[0 0;0 -1];1e-12];
.qt.approx["Ry(pi) |0> = |1>";{qApply[.quantQ.quantum.Ry[enlist[`theta]!enlist qPi];qKet0]};qKet1;1e-12];
.qt.approx["Rz(pi) |0> = exp(-i pi/2)|0>";{qApply[.quantQ.quantum.Rz[enlist[`theta]!enlist qPi];qKet0]};qs[0 0;-1 0];1e-12];
.qt.approx["U1(lambda) |0> = |0>";{qApply[.quantQ.quantum.U1[enlist[`lambda]!enlist qPi];qKet0]};qKet0;1e-12];
.qt.approx["U1(pi) |1> = -|1>";{qApply[.quantQ.quantum.U1[enlist[`lambda]!enlist qPi];qKet1]};qs[0 -1;0 0];1e-12];
.qt.approx["U2 is unitary: U2|0> has norm 1";{qNorm2 qApply[.quantQ.quantum.U2[`phi`lambda!0.3 1.1];qKet0]};1f;1e-12];
.qt.approx["U2(0,pi) = H";{qApply[.quantQ.quantum.U2[`phi`lambda!(0f;qPi)];qKet0]};qApply[qH;qKet0];1e-12];
.qt.eq["singleDict lists all single-qubit gates";{key .quantQ.quantum.singleDict};`PauliX`PauliY`PauliZ`One`Null`U2`U1`H`S`Sconj`T`Tconj`Rx`Ry`Rz];

/ products
.qt.approx["multConstMat: i*X";{.quantQ.quantum.multConstMat[`re`im!0 1f;.quantQ.quantum.PauliX]};qm[(0 0;0 0);(0 1;1 0)];1e-12];
.qt.approx["tensorProd of identities is the 4x4 identity";{.quantQ.quantum.tensorProd[.quantQ.quantum.one;.quantQ.quantum.one]};qm[(1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1);4#enlist 0 0 0 0];1e-12];
.qt.approx["tensorProd[X;I] flips the most significant qubit: |00> -> |10>";{qApply[.quantQ.quantum.tensorProd[.quantQ.quantum.PauliX;.quantQ.quantum.one];qs[1 0 0 0;0 0 0 0]]};qs[0 0 1 0;0 0 0 0];1e-12];

/ computational basis
.qt.eq["getCompState lists the basis states";{.quantQ.quantum.getCompState 2};(0 0;0 1;1 0;1 1)];
.qt.approx["q2number reads bits most-significant first";{.quantQ.quantum.q2number 1 0 1};5f;0f];
.qt.approx["q2nMap numbers the basis states 0..2^n-1";{value .quantQ.quantum.q2nMap 3};til 8;0f];
.qt.approx["setSystem starts in |00>";{.quantQ.quantum.setSystem 2};qs[1 0 0 0;0 0 0 0];1e-12];
.qt.approx["setSystem on one qubit";{.quantQ.quantum.setSystem 1};qKet0;1e-12];

/ the bench (circuit)
.qt.eq["setBench starts with one identity matrix";{count .quantQ.quantum.setBench[2]`bench};1];
.qt.eq["setBench view has one column per qubit";{cols .quantQ.quantum.setBench[2]`view};`index`quibit0`quibit1];
qCirc:{[n;gates] b:.quantQ.quantum.setBench n; b:{[b;g] .quantQ.quantum.appendSingleQuibitMat[b;g 0;g 1;()!()]}/[b;gates]; b};
.qt.eq["appendSingleQuibitMat adds a matrix and a view row";{b:qCirc[2;enlist (`PauliX;0)]; (count b`bench;count b`view)};2 2];
.qt.approx["X on qubit 0 (least significant): |00> -> |01>";{qProbs .quantQ.quantum.evolveState[.quantQ.quantum.setSystem 2;qCirc[2;enlist (`PauliX;0)]]};0 1 0 0f;1e-12];
.qt.approx["X on qubit 1: |00> -> |10>";{qProbs .quantQ.quantum.evolveState[.quantQ.quantum.setSystem 2;qCirc[2;enlist (`PauliX;1)]]};0 0 1 0f;1e-12];
.qt.approx["X on qubit 2 of 3: |000> -> |100>";{qProbs .quantQ.quantum.evolveState[.quantQ.quantum.setSystem 3;qCirc[3;enlist (`PauliX;2)]]};0 0 0 0 1 0 0 0f;1e-12];
.qt.approx["X on qubit 1 of 3: |000> -> |010>";{qProbs .quantQ.quantum.evolveState[.quantQ.quantum.setSystem 3;qCirc[3;enlist (`PauliX;1)]]};0 0 1 0 0 0 0 0f;1e-12];
.qt.approx["H on every qubit gives the uniform superposition";{qProbs .quantQ.quantum.evolveState[.quantQ.quantum.setSystem 3;qCirc[3;((`H;0);(`H;1);(`H;2))]]};8#0.125;1e-12];
.qt.approx["an unknown gate name defaults to the identity";{qProbs .quantQ.quantum.evolveState[.quantQ.quantum.setSystem 2;qCirc[2;enlist (`NotAGate;0)]]};1 0 0 0f;1e-12];

/ CNOT and entanglement
.qt.eq["cnotInd: control 1, target 0 swaps |10> and |11>";{.quantQ.quantum.cnotInd[2;1;0]};(2 3;3 2)];
.qt.approx["CNOT flips the target when the control is set";{b:.quantQ.quantum.appendCNOTMat[qCirc[2;enlist (`PauliX;1)];1;0]; qProbs .quantQ.quantum.evolveState[.quantQ.quantum.setSystem 2;b]};0 0 0 1f;1e-12];
.qt.approx["CNOT does nothing when the control is clear";{b:.quantQ.quantum.appendCNOTMat[.quantQ.quantum.setBench 2;1;0]; qProbs .quantQ.quantum.evolveState[.quantQ.quantum.setSystem 2;b]};1 0 0 0f;1e-12];
.qt.approx["H then CNOT makes a Bell state";{b:.quantQ.quantum.appendCNOTMat[qCirc[2;enlist (`H;1)];1;0]; qProbs .quantQ.quantum.evolveState[.quantQ.quantum.setSystem 2;b]};0.5 0 0 0.5;1e-12];
.qt.eq["appendCNOTMat labels control and target in the view";{b:.quantQ.quantum.appendCNOTMat[.quantQ.quantum.setBench 2;1;0]; value last b`view};(1;`tCNOT;`cCNOT)];
.qt.eq["doubleDict lists CNOT";{key .quantQ.quantum.doubleDict};enlist `CNOT];

/ measurement
.qt.approx["measurement probabilities sum to 1";{sum qProbs .quantQ.quantum.evolveState[.quantQ.quantum.setSystem 2;qCirc[2;((`H;0);(`T;0);(`H;1))]]};1f;1e-12];
.qt.eq["measurement lists the basis states";{exec state from .quantQ.quantum.measurement .quantQ.quantum.setSystem 2};(0 0;0 1;1 0;1 1)];

.qt.finish[];
