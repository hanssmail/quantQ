/ Unit tests for lib/quantQ_art.q -- generative art on a canvas (.quantQ.art)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`art;

/ 4x4 canvas, all colour 1
artCanvas:update colour:1 from flip `x`y!flip (til 4) cross til 4;

/ a shape much larger than the canvas covers every cell wherever it is placed, so every colour flips
.qt.true["placeCircle with a huge diameter inverts every cell";{all -1=exec colour from .quantQ.art.placeCircle[100;artCanvas]}];
.qt.true["placeSquare with a huge side inverts every cell";{all -1=exec colour from .quantQ.art.placeSquare[100;artCanvas]}];
.qt.true["placeCircleGradient with a huge diameter inverts every cell";{all -1=exec colour from .quantQ.art.placeCircleGradient[100;artCanvas]}];
.qt.eq["placeCircle keeps the canvas size";{count .quantQ.art.placeCircle[1;artCanvas]};16];
.qt.true["placeCircle only produces colours +-1";{all (exec colour from .quantQ.art.placeCircle[1;artCanvas]) in -1 1}];
.qt.approx["placeCircle with diameter 0 flips exactly one (integer-positioned) cell";{sum 1=exec colour from .quantQ.art.placeCircle[0;artCanvas]};15;0f];

/ pivotColours
artPiv:([] x:0 0 1 1; y:0 1 0 1; colour:1 -1 -1 1);
.qt.eq["pivotColours splits y by colour and fills gaps with -1";{.quantQ.art.pivotColours artPiv};([] x:0 1 0 1; colour1:0 1 -1 -1; colour2:-1 -1 1 0)];

.qt.finish[];
