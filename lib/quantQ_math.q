// collection of mathematical functions


// constants
.quantQ.math.pi:acos[-1];
.quantQ.math.e:exp[1];
.quantQ.math.eulerMascheroni:0.5772156649;
.quantQ.math.kninchin:2.6854520010;
.quantQ.math.glaisherKinkelin:1.2824271291;
.quantQ.math.goldenRato:(1+sqrt[5])%2.0;

// Hyperbolic functions
.quantQ.math.cosh:{[x] 
    :0.5*exp[x]+exp[neg x];
 };
 
.quantQ.math.sinh:{[x] 
    :0.5*exp[x]-exp[neg x];
 };
 
.quantQ.math.tanh:{[x]
    :.quantQ.math.sinh[x]%.quantQ.math.cosh[x];
 };
 
.quantQ.math.coth:{[x]
     :.quantQ.math.cosh[x]%.quantQ.math.sinh[x];
 };
 
.quantQ.math.sech:{[x]
    :1.0%.quantQ.math.cosh[x]; 
 };
 
.quantQ.math.csch:{[x]
    :1.0%.quantQ.math.sinh[x];
 };
 
.quantQ.math.asinh:{[x]
    :log[x+sqrt 1+x*x];
 };

.quantQ.math.acosh:{[x]
    :log[x+sqrt neg[1]+x*x];
 };
 
.quantQ.math.atanh:{[x]
    :0.5*log[(1.0+x)%(1.0-x)];
 }; 

.quantQ.math.acoth:{[x]
    :0.5*log[(x+1.0)%(x-1.0)];
 };
 
.quantQ.math.asech:{[x]
    :log[(1.0%x)+sqrt neg[1.0]+1.0%x*x];
 };
 
.quantQ.math.acsch:{[x]
    :log[(1.0%x)+sqrt 1.0+1.0%x*x];
 };
 
///////////////////////////////////////////// 
// Special functions -- real domain
/////////////////////////////////////////////
.quantQ.math.gammaApprox:{[x;n]
    // Euler expansion: Gamma(x) = (1/x) prod_k (1+1/k)^x / (1+x/k); converges like 1/n
    :(1.0%x)*prd {[x;n] xexp[(1.0+1.0%n);x]%(1.0+x%n) }[x;] 1+til n;               
 }; 

// Lanczos approximation of the gamma function (g=7, 9 coefficients), accurate to ~1e-15
.quantQ.math.lanczosCoef:0.99999999999980993 676.5203681218851 -1259.1392167224028 771.32342877765313 -176.61502916214059 12.507343278686905 -0.13857109526572012 9.9843695780195716e-6 1.5056327351493116e-7;

.quantQ.math.gamma:{[x]
    // x -- real argument (not a non-positive integer)
    // reflection formula for x < 1/2
    if[x<0.5; :.quantQ.math.pi%(sin .quantQ.math.pi*x)*.z.s 1-x];
    x:x-1;
    a:first[.quantQ.math.lanczosCoef]+sum (1_.quantQ.math.lanczosCoef)%x+1+til 8;
    t:x+7.5;
    :(sqrt 2*.quantQ.math.pi)*(xexp[t;x+0.5])*(exp neg t)*a;
 };
 
.quantQ.math.genLaguerrePol:{[x;alpha;n] 
    // generalised Laguerre polynomials 
    // example: .quantQ.math.genLaguerrePol[1.2;1.1;2]
    // recurrence: k L_k = (2k-1+alpha-x) L_(k-1) - (k-1+alpha) L_(k-2)
    if[n=0;:1.0];
    if[n=1;:1.0+alpha-x];
    l:(1.0;1.0+alpha-x);
    k:2;
    while[k<=n;
        l:(l 1;((((((2*k)-1)+alpha)-x)*l 1)-((k-1)+alpha)*l 0)%k);
        k+:1
    ];
    :l 1;
 };
 
.quantQ.math.gammaIncompleteLowerApprox:{[z;x;n]
     // lower incomplete gamma by its power series (converges for every x>=0), n terms
     // example: .quantQ.math.gammaIncompleteLowerApprox[2.3;1.0;12]
     :xexp[x;z]*exp[neg x]*(1.0%z)*sum prds 1.0,x%z+1+til n-1;
 };

.quantQ.math.gammaIncompleteUpperApprox:{[z;x;n]
     // upper incomplete gamma = Gamma(z) - lower incomplete gamma
     // example: .quantQ.math.gammaIncompleteUpperApprox[2.3;1.0;12]
     :.quantQ.math.gamma[z]-.quantQ.math.gammaIncompleteLowerApprox[z;x;n];
 };
 
.quantQ.math.physHermitePol:{[x;n]
    // x -- variable
    // n -- order of polynomial
    // using recurrence relation;Physics
    if[n=0;:1.0];
    if[n=1;:2*x];
    if[n>1;:(2*x*.z.s[x;n-1])-2*(n-1)*.z.s[x;n-2]];
 };
 
.quantQ.math.probHermitePol:{[x;n]
    // x -- variable
    // n -- order of polynomial
    // using recurrence relation;Probability
    if[n=0;:1.0];
    if[n=1;:x];
    if[n>1;:(x*.z.s[x;n-1])-(n-1)*.z.s[x;n-2]];
 };
 
.quantQ.math.firstChebyshevPol:{[x;n]
    // x -- variable
    // n -- order of polynomial
    if[n=0;:1.0];
    if[n=1;:x];
    if[n>1;:(2.0*x*.z.s[x;n-1])-.z.s[x;n-2]];
 };

.quantQ.math.secondChebyshevPol:{[x;n]
    // x -- variable
    // n -- order of polynomial
    if[n=0;:1.0];
    if[n=1;:2.0*x];
    if[n>1;:(2.0*x*.z.s[x;n-1])-.z.s[x;n-2]];
 }; 

.quantQ.math.getChebyshevFirstRoots:{[n]
    // n -- order of the polynomial to find the roots
    :{[n;k] cos[(.quantQ.math.pi*(k+0.5))%n]}[n;] each til n;
 }; 
 
.quantQ.math.getChebyshevSecondRoots:{[n]
    // n -- order of the polynomial to find the roots: cos(k pi/(n+1)), k=1..n
    :{[n;k] cos[.quantQ.math.pi*k%n+1]}[n;] each 1+til n;
 }; 

///////////////////////////////////////////// 
// PDF and CDF
/////////////////////////////////////////////
.quantQ.math.expPDF:{[x;lambda]
    x:"f"$x;
    :?[x>=0;lambda*exp[neg lambda*x];neg[x]*0];
 };

.quantQ.math.expCDF:{[x;lambda]
    x:"f"$x;
    :?[x>=0;1-exp[neg lambda*x];neg[x]*0];
 };
 
.quantQ.math.erlangPDF:{[x;k;K]
    // k -- shape, positive integer
    // K -- scale
    // x -- value
    :xexp[(1.0%K);k]*xexp[x;k-1]*exp[neg x%K]%.quantQ.stats.factorial[k-1];
 };
 
.quantQ.math.erlangCDF:{[l;k;x;n]
    // l -- scale
    // k -- shape, positive integer
    // x -- value
    // n -- number of series terms used
    :.quantQ.math.gammaIncompleteLowerApprox[k;x%l;n]%.quantQ.stats.factorial[k-1];
 };

.quantQ.math.gammaPDF:{[x;k;K;n]
    // k -- shape
    // K -- scale
    // x -- value
    // n -- unused (kept for compatibility); the gamma function is evaluated with .quantQ.math.gamma
    :reciprocal[.quantQ.math.gamma[k]*xexp[K;k]]*xexp[x;k-1]*exp[neg x%K];
 };

.quantQ.math.gammaCDF:{[x;k;K;n]
    // k -- shape
    // K -- scale
    // x -- value
    // n -- order of approximation used
    // n -- number of series terms used for the incomplete gamma function
    :reciprocal[.quantQ.math.gamma[k]]
        *.quantQ.math.gammaIncompleteLowerApprox[k;x%K;n]; 
 };
 
 .quantQ.math.normalPDF:{[x;mu;sigma2]
    // x -- value
    // mu -- mean
    // sigma2 -- variance
    x:"f"$x;
    : (1.0%sqrt 4*asin[1]*sigma2)*exp[neg (1.0%2*sigma2)*t*t:(x-mu)]
 };
 
 .quantQ.math.normalCDF:{[x;mu;sigma2]
    // x -- value
    // mu -- mean
    // sigma2 -- variance
    sigma:sqrt[sigma2];
    x:"f"$x;
    // hardcoded step size
    step:sigma%10000;
    // starting point, 8 standard deviations below the mean
    sPoint:mu-8*sigma;
    if[x<=sPoint;:0f];
    // the CDF (midpoint-rule integration of the density)
    :step*sum .quantQ.math.normalPDF[sPoint+step*0.5+til floor (x-sPoint)%step;mu;sigma2];
 };
       
.quantQ.math.normalMultiPDF:{[x;mu;sigma2]
    // x -- value
    // mu -- array of means
    // sigma2 -- covariance matrix
    x:"f"$x;
    // dimension of the x
    dimX: count x;
    // factor
    fact:xexp[(4*asin[1]);neg dimX%2.0]%sqrt .quantQ.mat.det[sigma2];
    // output
    : fact*exp[neg[0.5]*sum (x-mu)*inv[sigma2] mmu (x-mu)];
 };    
//  x:(0.2;0.2); mu:(0.0;0.0); sigma2:((1.0;0.2);(0.2;1.0)); .quantQ.math.normalMultiPDF[x;mu;sigma2]      

    
// Approximating Multivartiate CDF using Tetrachorix expansion
// Bivariate Tetrachoric expansion
.quantQ.math.tetrachoricExpansionBivar:{[x;rho;n]
    // x -- 2-dim array at which we evaluate the expansion
    // rho -- correlation between two N(0,1) processes
    // n -- order of expansion

    // CDF part
    out:.quantQ.math.normalCDF[x[0];0.0;1.0]*.quantQ.math.normalCDF[x[1];0.0;1.0];
    // add expansion in rho: sum_(k=1..n) rho^k/k! He_(k-1)(x) He_(k-1)(y)
    out+: .quantQ.math.normalPDF[x[0];0.0;1.0]*.quantQ.math.normalPDF[x[1];0.0;1.0]*
        sum {[x;rho;k] :(xexp[rho;k]%.quantQ.stats.factorial[k]) *
               .quantQ.math.probHermitePol[x[0];k-1] *
               .quantQ.math.probHermitePol[x[1];k-1]   
        }[x;rho;] each 1+til n;
    // output
    :out;
 };
// Example x:(0.0;0.0); rho:0.3;n:4; .quantQ.math.tetrachoricExpansion[x;rho;n]

// TODO:
// Trivariate Tetrachoric expansion

// General Tetrachoric expansion
