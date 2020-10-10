.quantQ.reg.ridge:{[y;x;reg]
    // y -- array of dependent variables
    // x -- array of arrays of independent variables, constants has to be included
    // reg -- regularization parameter for Ridge regression
    whereMissing:raze (where y=0n; raze {where x=0n} each x);
    // remove the missing values
    y:y[(til count y) except whereMissing];
    x:{[whereMissing;x] x[(til count x) except whereMissing]}[whereMissing;] each x;
    // y^T X
    ytx: enlist y mmu flip[x];
    // X^T X
    xtx: x mmu flip[x];
    // regularization
    xtx: (xtx+reg*.quantQ.mat.diagMatrix[count xtx]);
    // solve the regularised OLS equation
    solution: ytx lsq xtx;
    // beta
    beta:first solution;
    // output bucket
    bucket:(`beta`nObservations`dfn`regParameter!(beta;count[y];count[beta];reg));
    :bucket;
 };

.quantQ.reg.ridgeInOut:{[y;x;reg;inOutFlag]
    // y -- array of dependent variables
    // x -- array of arrays of independent variables, constants has to be included
    // reg -- regularization parameter for Ridge regression
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
    // IN sample Ridge regression  
    regIN:.quantQ.reg.ridge[y;x;reg];
    // OUT sample prediction
    yPRED:sum regIN[`beta]*xOUT;
    RMSE:sqrt avg t*t:yOUT-yPRED;
    regOUT:(`RMSE`nObservationsOUT)!(RMSE;count[yOUT]);
    // output bucket
    bucket:regIN,regOUT;
    :bucket;
 };

.quantQ.reg.ridgeCV:{[y;x;reg;nCV]
    // y -- array of dependent variables
    // x -- array of arrays of independent variables, constants has to be included
    // reg -- regularization parameter for Ridge regression
    // nCV -- number of Cross-validation folds
    whereMissing:raze (where y=0n; raze {where x=0n} each x);
    // remove the missing values
    y:y[(til count y) except whereMissing];
    x:{[whereMissing;x] x[(til count x) except whereMissing]}[whereMissing;] each x;
    // create all the Flags for every fold
    foldFlag: (floor til[count y]%count[y]%nCV)=/:til nCV;
    // return list of coresponding RMSEs
    :(.quantQ.reg.ridgeInOut[y;x;reg;] each foldFlag)[`RMSE];
 };

.quantQ.reg.oneStepLAR:{[bucket]
    // bucket -- all the variables
    beta:bucket[`beta];
    betaUsed:bucket[`betaUsed];
    betaDir:bucket[`betaDir];
    y:bucket[`y];
    xS:bucket[`xS];
    r:bucket[`r];
    path:bucket[`path];
    continue:bucket[`continue];
    step:bucket[`step];
    // set the significance limit
    corThr:1.0%sqrt "f"$count y;
    // step 1, find first variable
    $[1=count path;[
        // find variable with largest correlation
        betaUsed:z=max z:abs cor[r;] each xS;
        // set the direction
        betaDir:"f"$betaUsed*signum[cor[r;] each xS];
        // update r
        r-:sum step*betaDir*xS;
        // update beta
        beta+:step*betaDir;
        // update path
        path,:([] step:1+last path[`step];beta: enlist beta)
    ];];
    // step 2++, there is space for improvement
    $[(1<count path) and (0=prd betaUsed);
        [
        // correlation with existing direction
        corDir: abs cor[r;sum betaDir*xS];
        // best correlation with the candidate directions
        corNew: max z:not[betaUsed]* abs cor[r;] each (),xS;
        corNewWhere: z=max z;
        // decide whether to continue
        $[corThr>max[(corDir;corNew)];
            // stop the algorithm
            continue:0;
            // continue algorithm, decide the direction
            [
            $[corDir>corNew;
            // old direction
                [
                // update r
                r-:sum step*betaDir*xS;
                // update beta
                beta+:step*betaDir;
                // update path
                path,:([] step:1+last path[`step];beta: enlist beta)
                ];
            // new direction
                [
                // new set of variables
                betaUsed:betaUsed or corNewWhere;
                // new direction
                newDir:.quantQ.ols.fit[y;xS where betaUsed][`beta];
                // normalize by L1 to make steps comparable
                newDir:newDir%(sum abs newDir);
                // new beta direction
                betaDir:(count[betaUsed]#0.0);
                betaDir[where 1=betaUsed]:newDir;
                // update r
                r-:sum step*betaDir*xS;
                // update beta
                beta+:step*betaDir;
                // update path
                path,:([] step:1+last path[`step];beta: enlist beta)
                ]
            ]
            ]
        ];
    ];
    ];
    // step 2++, there is no space for improvement
    $[(1<count path) and (1=prd betaUsed);
        [
        // correlation with existing direction
        corDir: abs cor[r;sum betaDir*xS];
        $[corThr>corDir;
            // stop the algorithm
            continue:0;
            // continue algorithm
            [
            // update r
            r-:sum step*betaDir*xS;
            // update beta
            beta+:step*betaDir;
            // update path
            path,:([] step:1+last path[`step];beta: enlist beta)
            ]
        ]
        ];
    ];
    // update bucket
    bucket[`beta]:beta;
    bucket[`betaUsed]:betaUsed;
    bucket[`betaDir]:betaDir;
    bucket[`r]:r;
    bucket[`path]:path;
    bucket[`continue]:continue;
    // return bucket
    :bucket;
 };

.quantQ.reg.LAR:{[tab;step]
    // tab -- table with data with (y,'xS)
    // step -- the step of the LAR algorithm
    // dependent variable - first column of tab
    y:tab first cols tab;
    // independent variables - remaining columns
    xS:  (),tab 1_cols tab;
    // create empty vector with beta
    beta: ((count flip tab)-1)#0.0;
    // vector with beta in a given iteration
    betaDir: ((count flip tab)-1)#0.0;
    // used independent variables;
    betaUsed:((count flip tab)-1)#0b;
    // demeaning y, return meanY
    y:y-meanY:avg y;
    // define residual variable
    r:y;
    // table with entire path
    path:([] step:0;beta:enlist beta);
    // define bucket
    bucket: ((`beta`betaUsed`betaDir`y`xS`r`path`continue`step)!(beta;betaUsed;betaDir;y;xS;r;path;1;step));
    // run the LAR algorithm 
    bucketResult:(.quantQ.reg.oneStepLAR/)[{x[`continue]};bucket];
    // result is the re-arranged path table
    pathTab: bucketResult[`path];
    // create ArcLength, add meanY, and individual betas and return table: the L1 arc-length of piecewise-linear LAR coefficient profile amounts to summing the L1 norms of the changes in coefficients from step to step
    :(select arcLength: sums sum each abs deltas beta, beta0:meanY from pathTab),'(eval 
        (?;pathTab;();0b;({`$"beta",string x} each 1+til count[xS])!{(`beta;::;x)} each til count[xS]));
 };

.quantQ.reg.LARInOut:{[tab;step;inOutFlag]
    // tab -- table with data with (y,'xS)
    // step -- the step of the LAR algorithm
    // inOutFlag -- atom/array to form in/out sample, 1=out
    // dependent variable - first column of tab
    y:tab first cols tab;
    // independent variables - remaining columns
    xS:  (),tab 1_cols tab;
    whereMissing:raze (where y=0n; raze {where x=0n} each xS);
    // remove the missing values
    y:y[(til count y) except whereMissing];
    xS:{[whereMissing;xS] xS[(til count xS) except whereMissing]}[whereMissing;] each xS;
    // inOutFlag is atom, in sample is randomly chosen inOutFlag fraction of sample
    if[(1=count[inOutFlag]);
        // split sample: OUT
        xSOUT:flip ceiling[inOutFlag*count y] _ flip xS;
        yOUT:ceiling[inOutFlag*count y] _ y;
        // split sample: IN
        xS:flip ceiling[inOutFlag*count y]#flip xS;
        y:ceiling[inOutFlag*count y]#y
    ];
    // inOutFlag is Boolean array with 1b being In/0b Out
    if[(1<count[inOutFlag]);
        inOutFlag:inOutFlag[(til count inOutFlag) except whereMissing];
        // split sample: OUT
        yOUT:y where inOutFlag;
        xSOUT:flip flip[xS] where inOutFlag;
        // split sample: IN
        y:y where not inOutFlag;
        xS:flip flip[xS] where not inOutFlag
    ];
    // run LAR algorithm on in-sample
    inLAR: .quantQ.reg.LAR[([] y),'flip ({`$"x",string[x]} each til count[xS])!xS;step];
    // evaluate out -sample for every row of inLAR and calculate RMSE
   :update RMSE:{[yOUT;xSOUT;betas] :sqrt avg t*t:yOUT- betas[0] +sum (1_betas)*xSOUT;}
       [yOUT;xSOUT;] each flip 1_value flip inLAR from inLAR;
 };

// LARS regression with CV
.quantQ.reg.LARCV:{[tab;step;nCV]
    // tab -- table with data with (y,'xS)
    // step -- the step of the LAR algorithm
    // nCV -- number of Cross-validation folds
    // dependent variable - first column of tab
    y:tab first cols tab;
    // independent variables - remaining columns
    xS:  (),tab 1_cols tab;
    whereMissing:raze (where y=0n; raze {where x=0n} each xS);
    // remove the missing values
    y:y[(til count y) except whereMissing];
    xS:{[whereMissing;xS] xS[(til count xS) except whereMissing]}[whereMissing;] each xS;
    // create all the Flags for every fold
    foldFlag: (floor til[count y]%count[y]%nCV)=/:til nCV;
    // run evalaution for every fold
    outCV: (,/){[f;nflag] show first nflag; eval (!;f[last nflag];();0b;(enlist `fold)!enlist first nflag) }[.quantQ.reg.LARInOut[([] y),'flip ({`$"x",string[x]} each til count[xS])!xS;step;];] each flip ((til count foldFlag);foldFlag);
    // create frame (avoid rounding due to arcLength being float) using step and nfolds
    frame: flip (`arcLength`fold)!flip (exec arcLength*step from select distinct "j"$arcLength%step from outCV) cross (til count foldFlag); 
    // aj outCV with frame for every fold separately
    fullCV: (,/){[frame;outCV;n] aj[`arcLength;`arcLength xasc select from frame where fold=n;select arcLength, RMSE from outCV where fold=n]  }[frame;outCV;] peach (til count foldFlag);
    // return RMSE as a function of arcLength averaged across folds
    :select avg RMSE by arcLength from fullCV;
 };
