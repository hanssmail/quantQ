/ Unit tests for lib/quantQ_math.q -- constants, hyperbolic and special functions, distributions (.quantQ.math)
/ Reference values from Python's math module and scipy.special / scipy.stats.
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`math;

/ constants
.qt.approx["pi";{.quantQ.math.pi};3.141592653589793;1e-15];
.qt.approx["e";{.quantQ.math.e};2.718281828459045;1e-15];
.qt.approx["golden ratio";{.quantQ.math.goldenRato};1.618033988749895;1e-15];

/ hyperbolic functions and their inverses
.qt.approx["cosh";{.quantQ.math.cosh 0.7};1.255169005630943;1e-12];
.qt.approx["sinh";{.quantQ.math.sinh 0.7};0.7585837018395334;1e-12];
.qt.approx["tanh";{.quantQ.math.tanh 0.7};0.6043677771171636;1e-12];
.qt.approx["coth";{.quantQ.math.coth 0.7};1.654621635802629;1e-12];
.qt.approx["sech";{.quantQ.math.sech 0.7};0.796705459992875;1e-12];
.qt.approx["csch";{.quantQ.math.csch 0.7};1.3182460914662975;1e-12];
.qt.approx["asinh";{.quantQ.math.asinh 0.7};0.6526665660823557;1e-12];
.qt.approx["acosh";{.quantQ.math.acosh 1.7};1.123230982587296;1e-12];
.qt.approx["atanh";{.quantQ.math.atanh 0.5};0.5493061443340548;1e-12];
.qt.approx["acoth";{.quantQ.math.acoth 2f};0.5493061443340548;1e-12];
.qt.approx["asech";{.quantQ.math.asech 0.5};1.3169578969248166;1e-12];
.qt.approx["acsch";{.quantQ.math.acsch 2f};0.48121182505960347;1e-12];
.qt.approx["asinh(sinh x) = x";{.quantQ.math.asinh .quantQ.math.sinh -1.3};-1.3;1e-12];

/ gamma function (Euler product; converges like 1/n)
.qt.approx["gammaApprox(4) = 3! = 6";{.quantQ.math.gammaApprox[4f;100000]};6f;1e-2];
.qt.approx["gammaApprox(0.5) = sqrt(pi)";{.quantQ.math.gammaApprox[0.5;100000]};1.7724538509055159;1e-3];

/ Lanczos gamma function (math.gamma)
.qt.approx["gamma(4) = 6";{.quantQ.math.gamma 4f};6f;1e-12];
.qt.approx["gamma(0.5) = sqrt(pi)";{.quantQ.math.gamma 0.5};1.7724538509055159;1e-12];
.qt.approx["gamma(0.3) uses the reflection formula";{.quantQ.math.gamma 0.3};2.991568987687591;1e-12];
.qt.approx["gamma(7.2)";{.quantQ.math.gamma 7.2};1050.317816662683;1e-9];
.qt.approx["lower + upper incomplete gamma = gamma";{.quantQ.math.gammaIncompleteLowerApprox[2.3;1f;30]+.quantQ.math.gammaIncompleteUpperApprox[2.3;1f;30]};.quantQ.math.gamma 2.3;1e-12];

/ orthogonal polynomials
.qt.approx["genLaguerrePol n=0";{.quantQ.math.genLaguerrePol[1.2;1.1;0]};1f;1e-12];
.qt.approx["genLaguerrePol n=1";{.quantQ.math.genLaguerrePol[1.2;1.1;1]};0.9;1e-12];
.qt.approx["genLaguerrePol n=2 (scipy eval_genlaguerre)";{.quantQ.math.genLaguerrePol[1.2;1.1;2]};0.25500000000000006;1e-12];
.qt.approx["genLaguerrePol n=3 (scipy eval_genlaguerre)";{.quantQ.math.genLaguerrePol[2f;0.5;3]};-0.8958333333333335;1e-12];
.qt.approx["physHermitePol H3(0.5) = 8x^3-12x";{.quantQ.math.physHermitePol[0.5;3]};-5f;1e-12];
.qt.approx["probHermitePol He3(0.5) = x^3-3x";{.quantQ.math.probHermitePol[0.5;3]};-1.375;1e-12];
.qt.approx["probHermitePol He4(1.5)";{.quantQ.math.probHermitePol[1.5;4]};-5.4375;1e-12];
.qt.approx["firstChebyshevPol T1(x) = x";{.quantQ.math.firstChebyshevPol[0.3;1]};0.3;1e-12];
.qt.approx["firstChebyshevPol T2(0.5) = -0.5";{.quantQ.math.firstChebyshevPol[0.5;2]};-0.5;1e-12];
.qt.approx["firstChebyshevPol T3(0.3)";{.quantQ.math.firstChebyshevPol[0.3;3]};-0.792;1e-12];
.qt.approx["secondChebyshevPol U2(0.5) = 0";{.quantQ.math.secondChebyshevPol[0.5;2]};0f;1e-12];
.qt.approx["secondChebyshevPol U3(0.3)";{.quantQ.math.secondChebyshevPol[0.3;3]};-0.984;1e-12];
.qt.approx["getChebyshevFirstRoots: roots of T3";{.quantQ.math.getChebyshevFirstRoots 3};0.8660254037844387 0 -0.8660254037844387;1e-12];
.qt.approx["getChebyshevSecondRoots: roots of U3 are cos(k pi/4), k=1..3";{.quantQ.math.getChebyshevSecondRoots 3};0.7071067811865476 0 -0.7071067811865475;1e-12];

/ incomplete gamma (Laguerre expansion)
.qt.approx["gammaIncompleteUpperApprox(2.3, 1) (scipy)";{.quantQ.math.gammaIncompleteUpperApprox[2.3;1f;12]};0.9445066805236956;1e-2];

/ distributions
.qt.approx["expPDF";{.quantQ.math.expPDF[1;2]};0.2706705664732254;1e-12];
.qt.approx["expPDF is 0 for x<0";{.quantQ.math.expPDF[-1 0 1;2]};0 2 0.2706705664732254;1e-12];
.qt.approx["expCDF";{.quantQ.math.expCDF[1;2]};0.8646647167633873;1e-12];
.qt.approx["erlangPDF (scipy gamma, a=3, scale=2)";{.quantQ.math.erlangPDF[1.5;3;2f]};0.0664265464792052;1e-12];
.qt.approx["erlangCDF (scipy gamma, a=2, scale=1)";{.quantQ.math.erlangCDF[1f;2;1f;20]};0.2642411176571153;1e-3];
.qt.approx["gammaPDF (scipy gamma, a=2.5, scale=1.5)";{.quantQ.math.gammaPDF[2f;2.5;1.5;100000]};0.20352667466866572;1e-3];
.qt.approx["gammaCDF (scipy gamma, a=2.5, scale=1.5)";{.quantQ.math.gammaCDF[2f;2.5;1.5;20]};0.2487882896338786;1e-2];
.qt.approx["normalPDF standard at 0";{.quantQ.math.normalPDF[0;0;1]};0.3989422804014327;1e-12];
.qt.approx["normalPDF with variance 4";{.quantQ.math.normalPDF[1;0;4]};0.17603266338214973;1e-12];
.qt.approx["normalCDF standard at 0";{.quantQ.math.normalCDF[0;0;1]};0.5;1e-3];
.qt.approx["normalCDF standard at 1.96";{.quantQ.math.normalCDF[1.96;0;1]};0.9750021048517795;1e-3];
.qt.approx["normalCDF with variance 4 at 2 = Phi(1)";{.quantQ.math.normalCDF[2;0;4]};0.8413447460685429;1e-3];
.qt.approx["normalMultiPDF (scipy multivariate_normal)";{.quantQ.math.normalMultiPDF[0.2 0.2;0 0f;(1 0.2;0.2 1)]};0.15711152075990206;1e-10];
.qt.approx["tetrachoric expansion: bivariate normal CDF at (0,0), rho=0.3";{.quantQ.math.tetrachoricExpansionBivar[0 0f;0.3;6]};0.29849334201033917;2e-3];
.qt.approx["tetrachoric expansion: bivariate normal CDF at (0.5,-0.2), rho=0.3";{.quantQ.math.tetrachoricExpansionBivar[0.5 -0.2;0.3;8]};0.3320262544201823;2e-3];

.qt.finish[];
