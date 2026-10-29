/ Unit tests for lib/quantQ_complex.q -- complex numbers, DFT and the Mandelbrot set (.quantQ.complex)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`complex;

cx:{[r;i] `re`im!"f"$(r;i)};
cxPi:acos -1;

/ arithmetic
.qt.approx["norm |3+4i| = 5";{.quantQ.complex.norm cx[3;4]};5f;1e-12];
.qt.approx["mult (1+2i)(3+4i) = -5+10i";{.quantQ.complex.mult[cx[1;2];cx[3;4]]};cx[-5;10];1e-12];
.qt.approx["mult by i rotates by 90 degrees";{.quantQ.complex.mult[cx[0;1];cx[2;0]]};cx[0;2];1e-12];
.qt.approx["div (-5+10i)/(3+4i) = 1+2i";{.quantQ.complex.div[cx[-5;10];cx[3;4]]};cx[1;2];1e-12];
.qt.approx["div then mult round-trips";{.quantQ.complex.mult[.quantQ.complex.div[cx[1.5;-0.7];cx[0.2;2.5]];cx[0.2;2.5]]};cx[1.5;-0.7];1e-12];
.qt.eq["conjugate flips the imaginary part";{.quantQ.complex.conjugate cx[1;2]};cx[1;-2]];

/ polar form, all quadrants and the axes
.qt.approx["polar of 1+i";{.quantQ.complex.polar cx[1;1]};`radius`angle!(sqrt 2f;cxPi%4);1e-12];
.qt.approx["polar of -1+i (2nd quadrant)";{.quantQ.complex.polar cx[-1;1]};`radius`angle!(sqrt 2f;3*cxPi%4);1e-12];
.qt.approx["polar of -1-i (3rd quadrant)";{.quantQ.complex.polar cx[-1;-1]};`radius`angle!(sqrt 2f;-3*cxPi%4);1e-12];
.qt.approx["polar of 2i";{.quantQ.complex.polar cx[0;2]};`radius`angle!(2f;cxPi%2);1e-12];
.qt.approx["polar of -2i";{.quantQ.complex.polar cx[0;-2]};`radius`angle!(2f;neg cxPi%2);1e-12];
.qt.true["polar of 0 has a null angle";{null (.quantQ.complex.polar cx[0;0])`angle}];
.qt.approx["polar2Canonical inverts polar";{.quantQ.complex.polar2Canonical .quantQ.complex.polar cx[-0.3;0.8]};cx[-0.3;0.8];1e-12];

/ exponential and logarithm
.qt.approx["realExp[pi] = -1";{.quantQ.complex.realExp cxPi};cx[-1;0];1e-12];
.qt.approx["exp(i*pi) = -1";{.quantQ.complex.exp cx[0;cxPi]};cx[-1;0];1e-12];
.qt.approx["exp(1) = e";{.quantQ.complex.exp cx[1;0]};cx[exp 1;0];1e-12];
.qt.approx["log(-1), principal branch = i*pi";{.quantQ.complex.log[cx[-1;0];0]};cx[0;cxPi];1e-12];
.qt.approx["log(-1), branch k=1 = 3i*pi";{.quantQ.complex.log[cx[-1;0];1]};cx[0;3*cxPi];1e-12];
.qt.approx["exp(log z) = z";{.quantQ.complex.exp .quantQ.complex.log[cx[0.4;-1.3];0]};cx[0.4;-1.3];1e-12];

/ real2Cx
.qt.eq["real2Cx on an atom";{.quantQ.complex.real2Cx 3};cx[3;0]];
.qt.eq["real2Cx on a vector";{.quantQ.complex.real2Cx 1 2};`re`im!(1 2f;0 0f)];

/ discrete Fourier transform, compared with numpy.fft.fft
cxSeries:{[r] ([] re:"f"$r; im:count[r]#0f)};
.qt.approx["dft of 1 2 3 4 matches numpy.fft";{.quantQ.complex.dft cxSeries 1 2 3 4};([] re:10 -2 -2 -2f; im:0 2 0 -2f);1e-9];
.qt.approx["dft of an odd-length series matches numpy.fft";{.quantQ.complex.dft cxSeries 1 0 2 -1 0.5};([] re:2.5 0.3454915028125263 0.9045084971874737 0.9045084971874737 0.3454915028125263; im:0 -1.2878274987298426 3.147062175031697 -3.147062175031697 1.2878274987298426);1e-9];
.qt.approx["dft of a unit impulse is flat";{.quantQ.complex.dft cxSeries 1 0 0 0};([] re:4#1f; im:4#0f);1e-12];
/ the index column is created by `update i from ...`, which q names x
.qt.eq["dftTs adds polar coordinates";{4#cols .quantQ.complex.dftTs .quantQ.complex.dft cxSeries 1 2 3 4};`re`im`radius`angle];
.qt.eq["dftTs adds a row-index column";{last value flip .quantQ.complex.dftTs .quantQ.complex.dft cxSeries 1 2 3 4};til 4];
.qt.approx["dftTab frequency grid";{exec f from .quantQ.complex.dftTab ([] t:0 1 2 3f; xt:1 2 3 4f)};0 1 2 3%3;1e-12];
.qt.approx["dftTab real parts";{exec re from .quantQ.complex.dftTab ([] t:0 1 2 3f; xt:1 2 3 4f)};10 -2 -2 -2f;1e-9];

/ Mandelbrot iteration
.qt.eq["fractal: 0 never escapes, iteration stops at the cap";{last .quantQ.complex.fractal cx[0;0]};51];
.qt.eq["fractal: 2 escapes after one step";{last .quantQ.complex.fractal cx[2;0]};1];
.qt.true["mandelbrot_: 0 is in the set";{(.quantQ.complex.mandelbrot_[0 0;0f;0f])`res}];
.qt.true["mandelbrot_: 2 is not in the set";{not (.quantQ.complex.mandelbrot_[0 0;2f;0f])`res}];
.qt.eq["mandelbrot returns one cell per grid point";{count .quantQ.complex.mandelbrot[4;4]};16];

.qt.finish[];
