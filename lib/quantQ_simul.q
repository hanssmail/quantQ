.quantQ.simul.stochVolIntraPatt:{[time;isFlag;params]
    // time -- time tick
    // isFlag -- boolean flag
    // params -- dictionary with parameters
    // define intraday time
    tIntra:("f"$"j"$time-params[`tIni])%(1000*"f"$"j"$params[`tFin]-params[`tIni]);
    // return intraday volatility pattern (point at given time)
    :$[isFlag=1b;params[`c1]+(params[`cOpen]*exp[neg params[`aOpen]*tIntra])
        +(params[`cClose]*exp[neg params[`aClose]*(1.0 - tIntra)]);1];
 };

.quantQ.simul.choleskyChB:{[A]
    // A -- matrix to be decomposed
    L:0f*A;
    i:j:0;
    do[floor (N*1+N:count A)%2;
        L[i;j]: $[b:i=j; sqrt A[i;j] - {x mmu x} j#L[j] ;
            reciprocal[L[j;j]] * A[i;j] - (i#L i)mmu i#L j];
        i:i+b;
        j:(j+nb)*nb:not b;
    ];
    :L;
 };

.quantQ.simul.multiNormal:{[N;nRows;mu;Sigma]
    // N -- number of variables
    // nRows -- number of observations
    // mu -- array of means
    // Sigma -- variance-covariance matrix
    // generate N times nRows array
    multinormalVariate:{[nRows;x].quantQ.simul.getNormalVariate[nRows]}[nRows;] each til N;
    // multiply the multinormal variate by the Cholesky lower diagonal matrix
    multinormalVariateCorr:(.quantQ.simul.choleskyChB[Sigma]) mmu multinormalVariate;
    // add a mean of each process
    multinormalVariateCorr:multinormalVariateCorr+mu;
    // create a table, variables are named var1, var2,...
    :flip ({`$"var",string x} each 1+til N)!multinormalVariateCorr;
};

.quantQ.simul.getNormalVariate:{[nRows]
    // nRows -- dimension of the array
    // the Box-Muller procedure; we have to call it (nRows%2) round-up times
    nCall:ceiling nRows%2;
    // nCall calls of .quantQ.simul.genBoxMuller and then taking the first nRows outcomes
    :(raze {.quantQ.simul.genBoxMuller[]} each til nCall)[til nRows];
 };

.quantQ.simul.genBoxMuller:{[x]
    // x -- two-dimensional array to transform
    // generate 2-dimensional array with radius within (0,1) using convergence
    x1x2:{[x] (2?2f)-1.0}/[{(((x[0]*x[0])+(x[1]*x[1]))>=1) or (((x[0]*x[0])+(x[1]*x[1]))=0.0)};(0.0;0.0)];
    // use it as an input into .quantQ.simul.BoxMuller
    :.quantQ.simul.BoxMuller[x1x2];
 };

.quantQ.simul.BoxMuller:{[x1x2]
    // x1x2 -- 2-dimensional array of uniform numbers
    // the radius of point in 2D space
    rad:(x1x2[0]*x1x2[0])+x1x2[1]*x1x2[1];
    // the output is a 2-dimensional array
    z1z2:(0 0f);
    // the first normal variable
    z1z2[0]:x1x2[0]*sqrt[(neg 2*log[rad])%rad];
    // the second normal variable
    z1z2[1]:x1x2[1]*sqrt[(neg 2*log[rad])%rad];
    :z1z2;
 };

.quantQ.simul.expRandomVariate:{[lambda]
    // lambda -- parameter of exponential distribution
    :neg[log[first 1?1f]]%lambda;
 };


















