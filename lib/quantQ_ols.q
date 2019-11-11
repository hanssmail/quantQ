.quantQ.ols.signLevels:{[ts]
    // ts -- list of t-statistics
    convertPDict:((0j;1j;2j;3j;4j)!("-";"* (90% c.l.)";"** (95% c.l.)";"*** (99% c.l.)";
    "**** (99.9% c.l.)"));
    :convertPDict each binr[(1.282;1.645;2.326;3.091);] abs ts;
 };

.quantQ.ols.rSquared:{[y;yHat]
    // y -- the underlying data
    // yHat -- the model 
    // Total Sum of Squares
    TSS:sum t*t:y-avg[y];
    // Sum Square of Residuals
    RSS: sum t*t:y-yHat;
    // R squared
    :(TSS-RSS)%TSS;
 };

.quantQ.ols.rSquaredAdj:{[y;yHat;p]
    // y -- the underlying data
    // yHat -- the model 
    // p -- the number of explanatory variables excluding constant
    // adjusted Total Sum of Squares
    adjTSS:(1.0%(count[y]-1))*sum t*t:y-avg[y];
    // adjusted Sum Square of Residuals
    adjRSS: (1.0%(count[y]-1-p))*sum t*t:y-yHat;
    // adjusted R squared
    :(adjTSS-adjRSS)%adjTSS;
 };

.quantQ.ols.logL:{[n;RSS]
    // n -- size of the sample
    // RSS -- Residual sum of squares
    :neg[log[2*acos-1]*n%2.0]+neg[log[sqrt[RSS%n]]*n]+neg[(1.0%(2.0%n))];
 };

.quantQ.ols.AIC:{[p;logL]
    // p -- number of parameters
    // logL -- log-likelihood
    :(2.0*p)-2.0*logL;
 };

.quantQ.ols.AICC:{[n;sample;logL]
    // n -- number of parameters
    // sample -- sample size
    // logL -- log-likelihood
    :.quantQ.ols.AIC[n;logL]+("f"$2*(n*n)+n)%("f"$sample-1-n);
 };

.quantQ.ols.fit:{[y;x]
    // y -- array of dependent variables
    // x -- array of arrays of independent variables, constants has to be included
    whereMissing:raze (where y=0n; raze {where x=0n} each x);
    // remove the missing values
    y:y[(til count y) except whereMissing];
    x:{[whereMissing;x] x[(til count x) except whereMissing]}[whereMissing;] each x;
        // y^T X
        ytx: enlist y mmu flip[x];
        // X^T X
        xtx: x mmu flip[x];
        // solve the OLS equation
        solution: ytx lsq xtx;
        // beta
        beta:first solution;
        // var cov matrix
        varcov:(inv flip[xtx])*sum (1.0%(count[y]-count[beta]))*t*t:y-sum[beta*x];
        // errors for every estimation (sample standard deviation)
        errors:{sqrt x[y;y]}[varcov;] each til count varcov;
        // t-stat for every beta
        tStats:beta%errors;
        // TSS  or Total Sum of Squares
        TSS: sum t*t:y-avg[y];
        // RSS or Residual Sum of Squars 
        RSS: sum t*t:y-sum beta*x;
        // R-squared
        Rsquared:1.0-RSS%TSS;
        // adjusted R-squared
        adjRsquared:1-(1-Rsquared)*(count[y]-1.0)%(count[y]-(count[beta]-1)-1.0);
        // significance level in normal approximation (econometric)
        convertPDict:((0j;1j;2j;3j;4j)!("-";"* (90% c.l.)";"** (95% c.l.)";"*** (99% c.l.)";
            "**** (99.9% c.l.)"));
        signLevels: convertPDict each binr[(1.282;1.645;2.326;3.091);abs tStats];
        // log-likelihood        
        logL:neg[log[2*acos -1]*count[y]%2.0]+neg[log[sqrt[RSS%count[y]]]*count[y]]+
            neg[(1.0%(2.0%count[y]))];
        // AIC from log-likelihood
        aic:(2.0*count[beta])-2.0*logL;
        // AIC with small sample correction - depends upon heavy assumptions
        aicc:aic+("f"$2*(k*k)+k)%("f"$count[y]-1-k:count[beta]);
    // output bucket
    bucket:(`beta`errors`varcov`nObservations`dfn`tStats`TSS`RSS`Rsquared`adjRsquared`signLevels`logL`aic`aicc!
           (beta;errors;varcov;count[y];count[beta];tStats;TSS;RSS;Rsquared;adjRsquared;signLevels;logL;aic;aicc));
    :bucket;
 };

.quantQ.olsTab:{[tab]
    // tab -- table with data
    // names
    yName:first cols tab;
    xNames:1_cols tab;
    // running regressions
    reg:.quantQ.ols.fit[tab[yName];tab[xNames]];
    // regression part of table
    tabReg:([] name:xNames; coeff:reg[`beta];error:reg[`errors];tStat:reg[`tStats];
        significance:reg[`signLevels]);
    // stats part of table
    tabStats:([] name:`STATISTICS`nObservations`nParameters`TSS`RSS`Rsquared`RsquaredAdjusted`logLikelihood`AIC`AICC; coeff:(0nf;"f"$count[y];"f"$reg[`dfn];
        "f"$reg[`TSS];"f"$reg[`RSS];reg[`Rsquared];reg[`adjRsquared];reg[`logL];reg[`aic];reg[`aicc]));
    // return the table with results
    :tabReg uj tabStats;
 };

.quantQ.ols.fStatistics:{[model1;model2]
    // model1 -- the regression outcome of the nested model
    // model2 -- the regression outcome of the basis model
    // F-statistics
    :((model2[`nObservations]-model2[`dfn])%(model2[`dfn]-model1[`dfn]))*
        (model1[`RSS]-model2[`RSS])%model2[`RSS];
 };

.quantQ.ols.predict:{[model;tabFeatures]
    // model -- outcome of the .quantQ.ols.fit function
    // tabFeatures -- table with features used to obtain model without y
    // create column in tabFeatures with Prediction
    :([] yPredicted:sum model[`beta]*tabFeatures[cols tabFeatures]),'tabFeatures;
 };

.quantQ.ols.RMSE:{[tabModel]
    // tabModel -- table where first column is true dependent variable and second is the estimated one
    :sqrt t wavg t:tabModel[cols[tabModel][0]]-tabModel[cols[tabModel][1]];
 };

.quantQ.ols.inOut:{[y;x;inOutFlag]
    // y -- array of dependent variables
    // x -- array of arrays of independent variables, constants has to be included
    // inOutFlag -- atom/array to form in/out sample, 1=out
    whereMissing:raze (where y=0n; raze {where x=0n} each x);
    // remove the missing values
    y:y[(til count y) except whereMissing];
    x:{[whereMissing;x] x[(til count x) except whereMissing]}[whereMissing;] each x;
    // inOutFlag is atom, in sample is randomly chosen inOutFlag fraction of sample
    if[(1=count[inOutFlag]);
        // split sample: OUT
        xOUT:flip ceiling[inOutFlag*count y] _ flip x;
        yOUT:ceiling[inOutFlag*count y] _ y;
        // split sample: IN
        x:flip ceiling[inOutFlag*count y]#flip x;
        y:ceiling[inOutFlag*count y]#y
    ];
    // inOutFlag is Boolean array with 1b being In/0b Out
    if[(1<count[inOutFlag]);
        inOutFlag:inOutFlag[(til count inOutFlag) except whereMissing];
        // split sample: OUT
        yOUT:y where inOutFlag;
        xOUT:flip flip[x] where inOutFlag;
        // split sample: IN
        y:y where not inOutFlag;
        x:flip flip[x] where not inOutFlag
    ];
    // IN sample regression  
    regIN:.quantQ.ols.fit[y;x];
    // OUT sample prediction
    yPRED:sum regIN[`beta]*xOUT;
    RMSE:sqrt avg t*t:yOUT-yPRED;
    regOUT:(`RMSE`nObservationsOUT)!(RMSE;count[yOUT]);
    // output bucket
    bucket:regIN,regOUT;
    :bucket;
 };

.quantQ.ols.cv:{[y;x;nCV]
    // y -- array of dependent variables
    // x -- array of arrays of independent variables, constants has to be included
    // nCV -- number of Cross-validation folds
    whereMissing:raze (where y=0n; raze {where x=0n} each x);
    // remove the missing values
    y:y[(til count y) except whereMissing];
    x:{[whereMissing;x] x[(til count x) except whereMissing]}[whereMissing;] each x;
    // create all the Flags for every fold
    foldFlag: (floor til[count y]%count[y]%nCV)=/:til nCV;
    // return list of coresponding RMSEs
    :(.quantQ.ols.inOut[y;x;] each foldFlag)[`RMSE];
 };




