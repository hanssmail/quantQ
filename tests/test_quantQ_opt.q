/ Unit tests for lib/quantQ_opt.q -- test functions and the Nelder-Mead (amoeba) optimiser (.quantQ.opt)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`opt;

/ test functions (values quoted in the library comments)
.qt.approx["rosenbrock at (-1.5,-0.5)";{.quantQ.opt.func.rosenbrock[1f;100f;-1.5 -0.5]};762.5;1e-10];
.qt.approx["rosenbrock at (1.5,-0.5)";{.quantQ.opt.func.rosenbrock[1f;100f;1.5 -0.5]};756.5;1e-10];
.qt.approx["rosenbrock global minimum at (a,a^2)";{.quantQ.opt.func.rosenbrock[2.5;100f;2.5 6.25]};0f;1e-12];
.qt.approx["himmelblau minimum at (3,2)";{.quantQ.opt.func.himmelblau 3 2f};0f;1e-12];
.qt.approx["himmelblau minimum at (-2.805118,3.131312)";{.quantQ.opt.func.himmelblau -2.805118 3.131312};0f;1e-9];

/ Nelder-Mead
optRb:{.quantQ.opt.amoeba[1e-10;(-1 1f;-1 0f;0.5 0f);.quantQ.opt.func.rosenbrock[1f;100f]]};
.qt.eq["amoeba returns pmin, fmin and nfunc";{key optRb[]};`pmin`fmin`nfunc];
.qt.approx["amoeba finds the Rosenbrock minimum (1,1)";{optRb[]`pmin};1 1f;1e-3];
.qt.true["amoeba: fmin is close to 0 for Rosenbrock";{1e-6>optRb[]`fmin}];
.qt.approx["amoeba: fmin equals the function at pmin";{r:optRb[]; .quantQ.opt.func.rosenbrock[1f;100f;r`pmin]};optRb[]`fmin;1e-12];
.qt.approx["amoeba on a shifted quadratic";{.quantQ.opt.amoeba[1e-12;(0 0f;1 0f;0 1f);{(t*t:x[0]-0.3)+2*u*u:x[1]+1.7}]`pmin};0.3 -1.7;1e-4];
.qt.approx["amoeba finds the Himmelblau minimum near (3,2)";{.quantQ.opt.amoeba[1e-12;(2.5 1.5;3.5 1.5;3 2.5);.quantQ.opt.func.himmelblau]`pmin};3 2f;1e-4];
.qt.approx["amoeba in 3 dimensions";{.quantQ.opt.amoeba[1e-12;(0 0 0f;1 0 0f;0 1 0f;0 0 1f);{sum t*t:x-1 2 3f}]`pmin};1 2 3f;1e-4];
.qt.true["amoeba counts function evaluations";{0<optRb[]`nfunc}];

.qt.finish[];
