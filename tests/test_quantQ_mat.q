/ Unit tests for lib/quantQ_mat.q -- linear algebra: QR, eigen-decomposition, PCA, determinants, matrix roots (.quantQ.mat)
/ Reference values from numpy.linalg.
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`mat;

matI3:(1 0 0f;0 1 0f;0 0 1f);
matA:(4 1 0f;1 3 1f;0 1 2f);
matB:(2 1f;1 2f);
matM4:(2 1 0 3f;1 4 2 0f;0 1 3 1f;5 0 1 2f);
matLower:{[m] raze {[m;i] i#m i}[m] each til count m};

/ identity and rotations
.qt.eq["diagMatrix";{.quantQ.mat.diagMatrix 3};matI3];
.qt.eq["diagMatrix2 equals diagMatrix";{.quantQ.mat.diagMatrix2 3};matI3];
.qt.approx["givensRotation in 2-d";{.quantQ.mat.givensRotation[2;0;1;0.3]};(cos 0.3;neg sin 0.3;sin 0.3;cos 0.3);1e-12];
.qt.approx["givensRotation is orthogonal";{g:.quantQ.mat.givensRotation[4;1;3;0.7]; g mmu flip g};.quantQ.mat.diagMatrix 4;1e-12];

/ QR decomposition via Givens rotations
.qt.approx["QRdecomp: rotated matrix is upper triangular";{matLower .quantQ.mat.QRdecomp[matM4]`matrixRotated};6#0f;1e-10];
.qt.approx["QRdecomp: rotation matrix is orthogonal";{q:.quantQ.mat.QRdecomp[matM4]`rotationMatrix; q mmu flip q};.quantQ.mat.diagMatrix 4;1e-10];
.qt.approx["QRdecomp: Q R reproduces the input";{r:.quantQ.mat.QRdecomp matM4; flip[r`rotationMatrix] mmu r`matrixRotated};matM4;1e-10];

/ eigen-decomposition of symmetric matrices
.qt.approx["eigenSystem: eigenvalues of [[2,1],[1,2]]";{.quantQ.mat.eigenSystem[matB]`eigenvalues};3 1f;1e-8];
.qt.approx["eigenSystem: eigenvalues of a 3x3 symmetric matrix (numpy)";{.quantQ.mat.eigenSystem[matA]`eigenvalues};4.7320508075688785 3 1.2679491924311228;1e-8];
.qt.approx["eigenSystem: eigenvectors, first component positive (numpy)";{.quantQ.mat.eigenSystem[matA]`eigenvectors};(0.7886751345948129 0.5773502691896258 0.2113248654051871;0.5773502691896258 -0.5773502691896257 -0.5773502691896255;0.21132486540518702 -0.5773502691896255 0.7886751345948129);1e-6];
.qt.approx["eigenSystem: A v = lambda v";{e:.quantQ.mat.eigenSystem matA; (matA mmu/: e`eigenvectors)-e[`eigenvalues]*e`eigenvectors};9#0f;1e-6];
.qt.approx["extractEigenvalQR returns the sorted diagonal";{.quantQ.mat.extractEigenvalQR enlist[`matrix]!enlist (1 5 0f;0 3 2f;0 0 2f)};3 2 1f;1e-12];
.qt.approx["inverseIteration recovers the eigenvector of a known eigenvalue";{.quantQ.mat.inverseIteration[matB;3f;1e-10;1000]};(sqrt 0.5;sqrt 0.5);1e-6];
.qt.approx["normL1 is the maximum absolute column sum";{.quantQ.mat.normL1 (1 -2f;3 4f)};6f;1e-12];

/ covariance and PCA (population covariance, as numpy.cov(bias=True))
matData:(1 2 3 4 5 6f;2 1 4 3 6 8f;0.5 0.1 0.3 0.9 0.7 0.2);
.qt.approx["covarianceMatrixFast equals the population covariance";{.quantQ.mat.covarianceMatrixFast matData};(2.9166666666666665 3.6666666666666665 0.07499999999999998;3.6666666666666665 5.666666666666666 0f;0.07499999999999998 0 0.07916666666666666);1e-10];
.qt.approx["pca: variances are the eigenvalues of the covariance";{exec variance from .quantQ.mat.pca[matData]`tab};8.20789290482536 0.38796166749027583 0.06664542768436248;1e-6];
.qt.approx["pca: relative variances sum to 1";{sum exec relVariance from .quantQ.mat.pca[matData]`tab};1f;1e-10];
.qt.eq["pca: transformed data has one series per component";{count .quantQ.mat.pca[matData]`dataTransformed};3];
.qt.approx["covEmpTab = X'X / n (no demeaning)";{.quantQ.mat.covEmpTab ([] a:1 3 5f; b:2 4 7f)};(11.666666666666666 16.333333333333332;16.333333333333332 23f);1e-10];

/ minors and determinants
.qt.eq["minorMatrix removes a row and a column";{.quantQ.mat.minorMatrix[matA;0;1]};(1 1f;0 2f)];
.qt.approx["det of a 1x1 matrix is a number";{.quantQ.mat.det enlist enlist 5f};5f;1e-12];
.qt.true["det of a 1x1 matrix is an atom";{-9h=type .quantQ.mat.det enlist enlist 5f}];
.qt.approx["det 2x2";{.quantQ.mat.det matB};3f;1e-12];
.qt.approx["det 3x3 (Sarrus)";{.quantQ.mat.det matA};18f;1e-10];
.qt.approx["det 4x4 (Laplace)";{.quantQ.mat.det matM4};-136f;1e-9];
.qt.approx["det of the identity";{.quantQ.mat.det .quantQ.mat.diagMatrix 5};1f;1e-12];

/ powers and square roots
.qt.approx["pow[A;3] = A A A";{.quantQ.mat.pow[matB;3]};matB mmu matB mmu matB;1e-12];
.qt.approx["pow[A;0] is the identity";{.quantQ.mat.pow[matB;0]};.quantQ.mat.diagMatrix 2;1e-12];
.qt.approx["DenmanBeavers: S S = A";{s:.quantQ.mat.DenmanBeavers[()!();matA]; s mmu s};matA;1e-8];
.qt.approx["DenmanBeavers: sqrt of a diagonal matrix";{.quantQ.mat.DenmanBeavers[()!();(4 0f;0 9f)]};(2 0f;0 3f);1e-10];
.qt.approx["powerSeriesSq: S S = A for A close to I";{s:.quantQ.mat.powerSeriesSq[()!();(1.2 0.1;0.1 0.9)]; s mmu s};(1.2 0.1;0.1 0.9);1e-8];

.qt.finish[];
