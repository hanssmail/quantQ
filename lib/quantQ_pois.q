// Poisson regression: y in integer, using Poisson processes
// Requires .quantQ.stats and .quantQ.opt packages

// Poisson distribution  - probability
.quantQ.pois.poissonPDF:{[lambda;n]
    // lambda -- Poisson process intrensity
    // n -- number of events to occur
    // re-cast the n into integer due to factorial
    n:"j"$n;
    :exp[neg lambda]*xexp[lambda;n]%.quantQ.stats.factorial[n];
 };

// Poisson distribution --- PDF
.quantQ.pois.poissonPDFTab:{[lambda;bucket]
    // lambda -- Poisson process intrensity
    // bucket -- parameters
    // define default 
    bucket:(enlist[`maxN]!enlist[20]),bucket;
    // return table with p(N=n|lambda)
    :{([] n:x;prob:y)}[z;] .quantQ.pois.poissonPDF[lambda;] each z:til bucket[`maxN];
 };

// Poisson distribution --- CDF
.quantQ.pois.poissonCDFTab:{[lambda;bucket]
    // lambda -- Poisson process intrensity
    // bucket -- parameters
    // define default 
    bucket:(enlist[`maxN]!enlist[20]),bucket;
    // return table with p(N=n|lambda)
    :update sums prob from .quantQ.pois.poissonPDFTab[lambda;bucket];
 };

// Negative value of the log-Likelihood for Poisson regression 
.quantQ.pois.logLPoissonNeg:{[tab;theta]  
    // tab -- yVar, xVar columns
    // theta -- vector of parameters
    // log-Likelihood excluding term not depending on theta 
    :neg exec sum (yVar*sum[theta*flip xVar])-exp[sum[theta*flip xVar]] from tab;        
 };

// Negative value of the log-Likelihood for Poisson regression with regularisation
.quantQ.pois.logLPoissonNegReg:{[tab;lambda;bucket;theta]  
    // tab -- yVar, xVar columns
    // theta -- vector of parameters
    // lambda -- L2 regularisation parameter
    // log-Likelihood excluding term not depending on theta 
    logL:neg exec sum (yVar*sum[theta*flip xVar])-exp[sum[theta*flip xVar]] from tab;        
    // add regularisation term -- do not regularise constant
    $[bucket[`addConstant]=1b;
        logL: logL - count[tab]*lambda*wsum[1_theta;1_theta];
        logL: logL - count[tab]*lambda*wsum[theta;theta]
    ];
    // regularised logLikelihood     
    :logL;
 };


// Poisson regression
.quantQ.pois.regPoisson:{[tab;bucket]
    // tab -- yVar, xVar columns
    // bucket -- dictionary with parameters
    // dimensionalioty of the problem
    xDim: count flip exec xVar from tab;
    // defaults
    bucket:((`precision`initPoint)!(0.0001;{[x;y] neg[1]+x?2.0}[xDim;] each til 1+xDim)),bucket;
    // define function to minimise (maximum likelihood)
    sol:.quantQ.opt.amoeba[bucket[`precision]; bucket[`initPoint];.quantQ.pois.logLPoissonNeg[tab;]]; 
    // add extected intensity and sorted list of expected events
    tab:update nExpected: enlist each {[x] exec n from `prob xdesc .quantQ.pois.poissonPDFTab[x;()!()]} each intensity 
        from update intensity: sum sol[`pmin]*flip xVar from tab;
    logLikelihood: sum exec (yVar*intensity)-exp[intensity] from tab;
    // return
    :(`theta`logLikelihood`output)!(sol[`pmin];logLikelihood;select yVar, xVar, intensity, yVarPRED: first each nExpected, yVarLIST: nExpected from tab);    
 };


// Predict the values using the regressed model
.quantQ.pois.predPoisson:{[tab;model]
    // tab -- yVar, xVar columns
    // model -- dictionary with the solution of the model
    // evaluate model on the provided data
    out: update logL: sum  (yVar*intensity)-exp[intensity],yVarPRED: yVarLIST[;0] from update yVarLIST: {[x] exec n from `prob xdesc .quantQ.pois.poissonPDFTab[x;()!()]} each intensity 
        from update intensity: sum model[`theta]*flip xVar from select yVar, xVar from tab;
    // adding correlation    
    out: update corr: cor[yVar;yVarPRED] from out;   
    // return the full table    
    :out;
 };

// Poisson regression, with L2 regularisation
.quantQ.pois.regPoissonReg:{[tab;lambda;bucket]
    // tab -- yVar, xVar columns
    // bucket -- dictionary with parameters
    // lambda -- regularisation parameter
    // defaults
    bucket:((`precision`addConstant)!(0.0001;1b)),bucket; 
    // dimensionalioty of the problem, constant
    xDim: bucket[`addConstant]+count flip exec xVar from tab;
    // add default for init point
    bucket:(enlist[`initPoint]!enlist[{[x;y] neg[1]+x?2.0}[xDim;] each til 1+xDim]),bucket;
    // add constant
    if[bucket[`addConstant]=1b;
        tab:update xVar: (unit,'xVar) from update unit:1.0 from tab;    
    ];
    // define function to minimise (maximum likelihood)
    sol:.quantQ.opt.amoeba[bucket[`precision]; bucket[`initPoint];.quantQ.pois.logLPoissonNegReg[tab;lambda;bucket;]]; 
    // add extected intensity and sorted list of expected events
    tab:update nExpected: enlist each {[x] exec n from `prob xdesc .quantQ.pois.poissonPDFTab[x;()!()]} each intensity 
        from update intensity: sum sol[`pmin]*flip xVar from tab;
    logLikelihood: sum exec (yVar*intensity)-exp[intensity] from tab;
    // return
    :(`theta`logLikelihood`output`lambda`addConstant)!(sol[`pmin];logLikelihood;
        select yVar, xVar, intensity, yVarPRED: first each nExpected, yVarLIST: nExpected from tab;lambda;bucket[`addConstant]);    
 };

// predict the model
.quantQ.pois.predPoissonReg:{[tab;model]
    // tab -- yVar, xVar columns
    // model -- dictionary with the solution of the model
    if[model[`addConstant]=1b;
        tab:update xVar: (unit,'xVar) from update unit:1.0 from tab;    
    ];
    // evaluate model on the provided data
    out: update logL: sum  (yVar*intensity)-exp[intensity],yVarPRED: yVarLIST[;0] from update yVarLIST: {[x] exec n from `prob xdesc .quantQ.pois.poissonPDFTab[x;()!()]} each intensity 
        from update intensity: sum model[`theta]*flip xVar from select yVar, xVar from tab;
    // adding correlation    
    out: update corr: cor[yVar;yVarPRED] from out;   
    // return the full table    
    :out;
 };

// CV fit where logL is being optimised
.quantQ.pois.regPoissonCV:{[tab;bucket;nCV;lambda]
    // tab -- yVar, xVar columns
    // bucket -- dictionary of parameters
    // nCV -- number of folds
    // lambda -- regularisation parameter
    // create cv flags
    foldFlag:not (floor til[count tab]%count[tab]%nCV)=/:til nCV;
    // cv: fits, returns the logL 
    :{[tab;bucket;lambda;cv]
        tabIn: tab where cv;
        tabOut: tab where not cv;
        // fit in-sample model
        modelIn: .quantQ.pois.regPoissonReg[tabOut;lambda;bucket];
        // return soft-margin on the out-sof-sample
        :first exec logL from .quantQ.pois.predPoissonReg[tabIn;modelIn]
    }[tab;bucket;lambda;] each foldFlag;
 };

////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////
// Examples

// First, we explore the Poisson distribution
// Poisson distribution -- probability for a given lambda and n
/ .quantQ.pois.poissonPDF[0.3;3]
// Table with PDF
/ .quantQ.pois.poissonPDFTab[0.3;()!()]
// Table with CDF
/ .quantQ.pois.poissonCDFTab[0.3;()!()]

// Poisson regression
// synthetic data -- Poisson process, three explanatory variables 
/ nTab:500;
/ betaTrue: (-1.0;0.5;2.0);

/ tab: update xVar: (x1,'x2,'x3) from ([] x1:neg[1]+nTab?2.0; x2:nTab?2.0;x3:nTab?3.0);
/ tab: update yVar: {[x] exec first n from .quantQ.pois.poissonCDFTab[x;()!()] where prob>=first 1?1.0 } each sum betaTrue*flip xVar from tab;

// Poisson regression
// Empty dictionary for defaults to be used
/ bucket:()!();
/ model: .quantQ.pois.regPoisson[tab;bucket];
/ model

// Predict the values of the model using the regressed model
/ .quantQ.pois.predPoisson[tab;model]

// Change default values -- decrease precision
/ bucket: enlist[`precision]!enlist[0.01];
/ model: .quantQ.pois.regPoisson[tab;bucket];
/ model

// Regression with regularisation
/ bucket:()!();
/ lambda:0.2;
/ modelReg:.quantQ.pois.regPoissonReg[tab;lambda;bucket];
/ modelReg

// Predict the regularised model
/ .quantQ.pois.predPoissonReg[tab;modelReg]

// Regression model without added constant
/ bucket:enlist[`addConstant]!enlist[0b];
/ lambda:0.2;
/ modelRegNC:.quantQ.pois.regPoissonReg[tab;lambda;bucket];
/ modelRegNC

// Predict the regularised model without added constant
/ .quantQ.pois.predPoissonReg[tab;modelRegNC]

// CV the regularisation, the log-likelihood is presented
/ lambdaArray:0.0 0.001 0.01 0.1 0.2 0.3;
/ bucket:()!();
/ nCV:6;
/ ([] lambda: lambdaArray; logLNeg: avg each .quantQ.pois.regPoissonCV[tab;bucket;nCV;] each lambdaArray)
