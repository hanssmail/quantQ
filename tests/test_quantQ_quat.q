/ Unit tests for lib/quantQ_quat.q -- quaternions (.quantQ.quat)
/ Hamilton product reference: (a1+b1i+c1j+d1k)(a2+b2i+c2j+d2k).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`quat;

qt:{[r;i;j;k] `r`i`j`k!"f"$(r;i;j;k)};
qtOne:qt[1;0;0;0]; qtI:qt[0;1;0;0]; qtJ:qt[0;0;1;0]; qtK:qt[0;0;0;1];
qtNeg:{[q] neg q};

/ Hamilton product on the basis: i^2 = j^2 = k^2 = ijk = -1
.qt.approx["i*i = -1";{.quantQ.quat.mult[qtI;qtI]};qtNeg qtOne;1e-12];
.qt.approx["j*j = -1";{.quantQ.quat.mult[qtJ;qtJ]};qtNeg qtOne;1e-12];
.qt.approx["k*k = -1";{.quantQ.quat.mult[qtK;qtK]};qtNeg qtOne;1e-12];
.qt.approx["i*j = k";{.quantQ.quat.mult[qtI;qtJ]};qtK;1e-12];
.qt.approx["j*k = i";{.quantQ.quat.mult[qtJ;qtK]};qtI;1e-12];
.qt.approx["k*i = j";{.quantQ.quat.mult[qtK;qtI]};qtJ;1e-12];
.qt.approx["j*i = -k";{.quantQ.quat.mult[qtJ;qtI]};qtNeg qtK;1e-12];
.qt.approx["j*1 = j";{.quantQ.quat.mult[qtJ;qtOne]};qtJ;1e-12];
.qt.approx["k*1 = k";{.quantQ.quat.mult[qtK;qtOne]};qtK;1e-12];
.qt.approx["general product";{.quantQ.quat.mult[qt[1;2;3;4];qt[0.5;-1;2;0.25]]};qt[-4.5;-7.25;-1;9.25];1e-12];

/ inverse, norm, conjugate
.qt.approx["q * inverse q = 1";{q:qt[1;2;3;4]; .quantQ.quat.mult[q;.quantQ.quat.inverse q]};qtOne;1e-12];
.qt.approx["inverse q * q = 1";{q:qt[0.5;-1;2;0.25]; .quantQ.quat.mult[.quantQ.quat.inverse q;q]};qtOne;1e-12];
.qt.approx["norm";{.quantQ.quat.norm qt[1;2;3;4]};sqrt 30f;1e-12];
.qt.eq["conjugate negates the vector part";{.quantQ.quat.conjugate qt[1;2;3;4]};qt[1;-2;-3;-4]];
.qt.approx["|q|^2 = q * conj(q)";{.quantQ.quat.mult[qt[1;2;3;4];.quantQ.quat.conjugate qt[1;2;3;4]]};qt[30;0;0;0];1e-12];
.qt.approx["polarDecomp: norm and unit quaternion";{p:.quantQ.quat.polarDecomp qt[0;3;0;4]; (p 0;.quantQ.quat.norm p 1)};5 1f;1e-12];
.qt.approx["crossProduct i x j = k";{.quantQ.quat.crossProduct[qtI;qtJ]};qtK;1e-12];
.qt.approx["crossProduct is antisymmetric";{.quantQ.quat.crossProduct[qt[0;1;2;3];qt[0;-1;0.5;2]]+.quantQ.quat.crossProduct[qt[0;-1;0.5;2];qt[0;1;2;3]]};qt[0;0;0;0];1e-12];

/ exponential and logarithm
.qt.approx["exp(0) = 1";{.quantQ.quat.exp qt[0;0;0;0]};qtOne;1e-12];
.qt.approx["exp(i*pi) = -1";{.quantQ.quat.exp qt[0;acos -1;0;0]};qtNeg qtOne;1e-12];
.qt.approx["exp of a general quaternion";{.quantQ.quat.exp qt[0.5;0.2;0.3;0.1]};qt[1.5346509696798105;0.3221039005928113;0.4831558508892168;0.16105195029640565];1e-12];
.qt.approx["log(exp q) = q";{.quantQ.quat.log .quantQ.quat.exp qt[0.5;0.2;0.3;0.1]};qt[0.5;0.2;0.3;0.1];1e-12];
.qt.approx["log of a positive real";{.quantQ.quat.log qt[2;0;0;0]};qt[log 2;0;0;0];1e-12];

/ pretty printing
.qt.eq["visual";{.quantQ.quat.visual qt[1.5;-2.5;3.5;-4.5]};"1.5-2.5i+3.5j-4.5k"];

.qt.finish[];
