/ Unit tests for lib/quantQ_ts.q -- ARMA simulation, estimation, forecasting and VAR models (.quantQ.ts)
/ Reference values from numpy, using population autocovariances gamma(h) = cov(x[0..n-h-1], x[h..n-1])
/ (q's cov) and the textbook Durbin-Levinson / innovations recursions (Brockwell & Davis).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`ts;

tsD:{(sin 0.3*x)+(0.5*cos 1.7*x)+0.1*sin 5.1*x} "f"$til 200;
tsEps:1 0.5 -0.3 0.2 0 0.4;

/ simulation
.qt.approx["simAR: x_t = phi1 x_(t-1) + phi2 x_(t-2) + eps_(t-p), zero start";{.quantQ.ts.simAR[0.5 0.2;tsEps]};0 0 1 1 0.4 0.6;1e-12];
.qt.approx["simAR with phi=0 copies the shocks";{.quantQ.ts.simAR[enlist 0f;1 2 3 4f]};0 1 2 3f;1e-12];
/ MA(q) impulse responses: a unit shock at position q produces 1, alpha1, ..., alphaq
.qt.approx["simMA(1) impulse response";{.quantQ.ts.simMA[1;enlist 0.5;0 1 0 0 0f]};1 0.5 0 0 0f;1e-12];
.qt.approx["simMA(2) impulse response";{.quantQ.ts.simMA[2;0.5 0.25;0 0 1 0 0 0f]};1 0.5 0.25 0 0 0f;1e-12];
.qt.approx["simARMA(1,1)";{.quantQ.ts.simARMA[enlist 0.5;enlist 0.3;tsEps]};0 0.3 0.3 0.06 0.09 0.045;1e-12];
.qt.approx["simulate dispatches to simAR";{.quantQ.ts.simulate[`AR;0.5 0.2;tsEps]};.quantQ.ts.simAR[0.5 0.2;tsEps];1e-12];
.qt.approx["simulate dispatches to simARMA";{.quantQ.ts.simulate[`ARMA;(enlist 0.5;enlist 0.3);tsEps]};0 0.3 0.3 0.06 0.09 0.045;1e-12];
.qt.throws["simulate rejects integer parameters";{.quantQ.ts.simulate[`AR;1 2;tsEps]}];
.qt.throws["simulate rejects an unknown process";{.quantQ.ts.simulate[`GARCH;0.5 0.2;tsEps]}];

/ autocovariances and Yule-Walker
.qt.approx["gamma(0..3)";{.quantQ.ts.gamma[;tsD] each til 4};0.6262645597240284 0.46219954609568636 0.2891974451444185 0.3567328184506532;1e-12];
.qt.approx["r(1) = gamma(1)/gamma(0)";{.quantQ.ts.r[tsD;1]};0.7380260289666728;1e-12];
.qt.approx["toeplitz";{.quantQ.ts.toeplitz[tsD;3]};(1 0.7380260289666728 0.46178159158783805;0.7380260289666728 1 0.7380260289666728;0.46178159158783805 0.7380260289666728 1);1e-12];
.qt.approx["yuleWalker, p=2";{.quantQ.ts.yuleWalker[tsD;2]};0.8724003017453077 -0.18207253877857843;1e-10];
.qt.approx["DurbinLevinson, p=1";{.quantQ.ts.DurbinLevinson[tsD;1]};enlist 0.7380260289666728;1e-12];
.qt.approx["DurbinLevinson, p=2";{.quantQ.ts.DurbinLevinson[tsD;2]};0.8724003017453077 -0.18207253877857843;1e-10];
.qt.approx["DurbinLevinson, p=3";{.quantQ.ts.DurbinLevinson[tsD;3]};0.9969473974625791 -0.7788397647103619 0.6840520627261386;1e-10];
.qt.approx["DurbinLevinson agrees with Yule-Walker, p=4";{.quantQ.ts.DurbinLevinson[tsD;4]};.quantQ.ts.yuleWalker[tsD;4];1e-9];
.qt.approx["estARP: OLS on lags (intercept first)";{.quantQ.ts.estARP[tsD;2]};0.010471422337115729 0.8678799878027227 -0.17870392132582935;1e-10];
.qt.approx["innovations, q=1";{.quantQ.ts.innovations[tsD;1]};enlist 0.6781576819514514;1e-10];
.qt.approx["innovations, q=2";{.quantQ.ts.innovations[tsD;2]};0.7948329986353043 0.4228498075243186;1e-10];
.qt.approx["estimate dispatches AR to Durbin-Levinson";{.quantQ.ts.estimate[tsD;`AR;2]};0.8724003017453077 -0.18207253877857843;1e-10];
.qt.approx["estimate accepts the AR order as a one-element list";{.quantQ.ts.estimate[tsD;`AR;enlist 2]};0.8724003017453077 -0.18207253877857843;1e-10];
.qt.approx["estimate dispatches MA to the innovations algorithm";{.quantQ.ts.estimate[tsD;`MA;1]};enlist 0.6781576819514514;1e-10];
.qt.approx["estimate accepts the MA order as a one-element list";{.quantQ.ts.estimate[tsD;`MA;enlist 1]};enlist 0.6781576819514514;1e-10];
.qt.eq["HannanRissanen returns p+q coefficients";{count .quantQ.ts.HannanRissanen[tsD;1;1]};2];
.qt.throws["estimate rejects an unknown process";{.quantQ.ts.estimate[tsD;`GARCH;1]}];

/ forecasting
.qt.approx["forecastAR, p=1, three steps (re-estimated after each step)";{.quantQ.ts.forecastAR[tsD;3;1]};0.1825023456311589 0.13471177985059515 0.0994418428898235;1e-10];
.qt.approx["forecastAR, p=2: phi1 multiplies the latest value";{.quantQ.ts.forecastAR[tsD;1;2]};0.2633715363333801;1e-10];
.qt.approx["forecast dispatches AR";{.quantQ.ts.forecast[tsD;3;`AR;1]};0.1825023456311589 0.13471177985059515 0.0994418428898235;1e-10];
.qt.throws["forecast rejects a non-positive horizon";{.quantQ.ts.forecast[tsD;0;`AR;1]}];
.qt.true["forecast dispatches MA";{9h=type .quantQ.ts.forecast[tsD;1;`MA;1]}];

/ matrix helpers
.qt.eq["vec stacks the columns";{.quantQ.ts.vec (1 2;3 4)};1 3 2 4];
.qt.approx["kron (numpy.kron)";{.quantQ.ts.kron[(1 2f;3 4f);(0 1f;1 0f)]};(0 1 0 2;1 0 2 0;0 3 0 4;3 0 4 0f);0f];
.qt.eq["diag";{.quantQ.ts.diag (1 2 3;4 5 6;7 8 9)};1 5 9];
.qt.eq["eye";{.quantQ.ts.eye 2};(1 0f;0 1f)];
.qt.approx["pValueChi2 approximates the chi-square tail (scipy chi2.sf), chi=10, nu=5";{.quantQ.ts.pValueChi2[10f;5]};0.07523524614651217;1e-3];
.qt.approx["pValueChi2, chi=20, nu=9";{.quantQ.ts.pValueChi2[20f;9]};0.017912404529843298;1e-3];

/ VAR(1) on two deterministic series: OLS B = Y Z' (Z Z')^-1 with Z = (1; y_(t-1))
tsY:(60#tsD;{(cos 0.4*x)+0.2*sin 2.3*x} "f"$til 60);
.qt.approx["varpest: VAR(1) coefficients";{.quantQ.ts.varpest[tsY;1]};(-0.007267738967158783 0.7422922884266169 -0.028821176945655506;-0.023964569453215366 0.021211755660376225 0.8446693893509232);1e-8];
.qt.eq["estimate1: one-step forecast for each series";{count .quantQ.ts.estimate1[tsY;1]};2];

.qt.finish[];
