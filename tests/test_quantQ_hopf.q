/ Unit tests for lib/quantQ_hopf.q -- Hopfield networks and Hebbian descent (.quantQ.hopf)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`hopf;

hopS:1 -1 1 1f;
hopNoisy:1 -1 1 -1f;
hopT:-1 1 1 -1f;

/ network construction: W = s s' with zero diagonal
.qt.approx["getHopfieldNet is the outer product with zero diagonal";{.quantQ.hopf.getHopfieldNet 1 -1 1f};(0 -1 1f;-1 0 -1f;1 -1 0f);0f];
.qt.true["getHopfieldNet is symmetric";{w:.quantQ.hopf.getHopfieldNet hopS; all raze w=flip w}];
.qt.approx["getHopfieldNet2S of a state with itself equals getHopfieldNet";{.quantQ.hopf.getHopfieldNet2S[hopS;hopS]};.quantQ.hopf.getHopfieldNet hopS;0f];
.qt.approx["getHopfieldNet2S: s1 *\\: s2 with zero diagonal";{.quantQ.hopf.getHopfieldNet2S[1 -1 1f;1 1 -1f]};(0 1 -1f;-1 0 1f;1 1 0f);0f];

/ energies
.qt.approx["energy of the stored pattern is -n(n-1)/2";{.quantQ.hopf.getEnergyHopfield[.quantQ.hopf.getHopfieldNet hopS;hopS]};-6f;1e-12];
.qt.approx["getMinEnergyHopfield equals the energy of the stored pattern";{.quantQ.hopf.getMinEnergyHopfield hopS};-6;1e-12];
.qt.approx["deltaEnergyHopfield is 0 for no change";{.quantQ.hopf.deltaEnergyHopfield[.quantQ.hopf.getHopfieldNet hopS;hopS;hopS]};0f;1e-12];
.qt.approx["deltaEnergyHopfield for one flipped spin";{.quantQ.hopf.deltaEnergyHopfield[.quantQ.hopf.getHopfieldNet hopS;hopS;hopNoisy]};6f;1e-12];

/ dynamics
.qt.eq["signumHopfield uses >=0 -> +1";{.quantQ.hopf.signumHopfield -1 0 2f};-1 1 1f];
.qt.eq["the stored pattern is a fixed point";{.quantQ.hopf.oneHopfieldJump[.quantQ.hopf.getHopfieldNet hopS;hopS]};hopS];
.qt.eq["one update recovers the pattern from a one-bit corruption";{.quantQ.hopf.oneHopfieldJump[.quantQ.hopf.getHopfieldNet hopS;hopNoisy]};hopS];
.qt.eq["oneHopfieldJumpBias with zero bias equals oneHopfieldJump";{.quantQ.hopf.oneHopfieldJumpBias[.quantQ.hopf.getHopfieldNet hopS;hopNoisy;4#0f]};hopS];
.qt.eq["oneHopfieldJumpBias: a large bias dominates";{.quantQ.hopf.oneHopfieldJumpBias[.quantQ.hopf.getHopfieldNet hopS;hopS;-10 10 -10 10f]};-1 1 -1 1f];

/ ground state: an attractor with the minimum energy
.qt.approx["getGroundState reaches the minimum energy";{w:.quantQ.hopf.getHopfieldNet hopS; .quantQ.hopf.getEnergyHopfield[w;.quantQ.hopf.getGroundState w]};-6f;1e-12];

/ Hebbian descent: learn a network mapping input -> target
.qt.eq["runHebbDescent learns input -> target";{r:.quantQ.hopf.runHebbDescent[hopS;hopT;()!()]; .quantQ.hopf.oneHopfieldJumpBias[r`hopnet;hopS;r`bias]};hopT];
.qt.eq["runHebbDescent returns hopnet and bias";{key .quantQ.hopf.runHebbDescent[hopS;hopT;()!()]};`hopnet`bias];
.qt.eq["runHebbDescentBiased keeps the bias fixed";{(.quantQ.hopf.runHebbDescentBiased[hopS;hopT;0.1 0.2 0.3 0.4;()!()])`bias};0.1 0.2 0.3 0.4];
.qt.eq["runHebbDescentBiased learns input -> target";{b:0.1 0.2 0.3 0.4; r:.quantQ.hopf.runHebbDescentBiased[hopS;hopT;b;()!()]; .quantQ.hopf.oneHopfieldJumpBias[r`hopnet;hopS;b]};hopT];

.qt.finish[];
