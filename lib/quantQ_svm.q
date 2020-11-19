// SVM
/////////////////////////////////////////////
//
/////////////////////////////////////////////
// functions


// single observation, soft-margin
.quantQ.svm.softMargin:{[yVar;xVar;b;w;lambda]
    // yVar -- float
    // xVar -- array of x variables
    // b -- intercept
    // w -- array, normal vector of the hyperplane
    // lambda -- regularisation paramter (normalised by sample size)
    :(lambda*wsum[w;w])+{max[(0.0;x)]} 1-yVar*(neg[b]+sum w*xVar);    
 };

// single observation, soft-margin, gradients
.quantQ.svm.softMarginGradient:{[yVar;xVar;b;w;lambda]
    // yVar -- float
    // xVar -- array of x variables
    // b -- intercept
    // w -- array, normal vector of the hyperplane
    // lambda -- regularisation paramter (normalised by sample size)
    // ind -- index of the derivative: 0 b, 1+dim[w]      
    out:raze (0.0;2.0*lambda*w);
    if[(0<1-yVar*(neg[b]+sum w*xVar));
        out:out+raze (yVar;{[yVar;xVar;w;x] (neg yVar*xVar[x])}[yVar;xVar;w;] 
            each til count w)        
    ];
    :out;                    
 };

.quantQ.svm.fitSoftMargin:{[data;bucket]
    // data -- table with yVar (-1/1) and xVar (array)
    // outSVM:.quantQ.svm.fitSoftMargin[data;bucket]
    // default bucket
    bucketDefault:(`w`b)!(neg[1]+count[flip[data[`xVar]]]?2.0;0.0);
    // learning
    bucketDefault:bucketDefault,(`learningRateFunc`learningRateParams)!(
        .quantQ.nn.learningFuncs;(`method`learningRate!(`const;0.01)));
    // lambda (regularisation parameter), initialise counter
    bucketDefault:bucketDefault,(`lambda`counter)!(0.0;0);
    // continue and maximum number iof iterations      
    bucketDefault:bucketDefault,(`continue`stopReason`maxIters)!(1;`none;10000);
    // soft margin objective function initialisation
    bucketDefault:bucketDefault,enlist[`softMargin]!enlist[0wf];
    // decay of soft-margin -- slow and fast, difference below certain threshold
    bucketDefault:bucketDefault,(`sMSlow`sMFast`sMAlpSlow`sMAlpFast)!(0nf;0nf;0.8;0.3);
    // decay error margin, -0wf means never check
    bucketDefault:bucketDefault,enlist[`softMarginDecayError]!enlist[0.0001];
    // monitor convergence default
    bucketDefault:bucketDefault,enlist[`monitorConvergence]!(enlist[0]);
    // merge input with default
    bucket:bucketDefault,bucket; 
    // re-define learning rate function to be function of counter only
    bucket[`learningRateFunc]:(bucket[`learningRateFunc])[bucket[`learningRateParams];];
    // add data
    bucket:bucket,enlist[`data]!enlist[data];
    // terminal condition
    terminus:{ x[`continue]=1};
    // estimate and return
    $[bucket[`monitorConvergence]=1;
        :(.quantQ.svm.oneStep\)[terminus;bucket];
        :(.quantQ.svm.oneStep/)[terminus;bucket]
    ];
 };

// using sgd one step
.quantQ.svm.oneStep:{[bucket]
    // bucket -- dictionary of parameters
    // choose one point -- using SGD with one random observation
    dataPoint: exec from 1?bucket[`data];
    // learning rate
    learningRate:bucket[`learningRateFunc][bucket[`counter]];
    // derivative
    der:.quantQ.svm.softMarginGradient[dataPoint[`yVar];dataPoint[`xVar];bucket[`b];bucket[`w];bucket[`lambda]];
    updatedVec:raze[(bucket[`b];bucket[`w])]-learningRate*der;
    bucket[`b]:first updatedVec;
    bucket[`w]:1_updatedVec;
    // update soft margin objective function
    bucket[`softMargin]: avg .quantQ.svm.softMargin[;;bucket[`b];bucket[`w];bucket[`lambda]]'
        [exec yVar from data;exec xVar from data];
    // ewma of soft margin at two scales
    aSlow: 1%bucket[`sMAlpSlow]*count[bucket[`data]];
    aFast: 1%bucket[`sMAlpFast]*count[bucket[`data]];
    // initiate
    if[null bucket[`sMSlow];
        bucket[`sMSlow]:bucket[`softMargin];
        bucket[`sMFast]:bucket[`softMargin]
    ];
    // update two scales
    bucket[`sMSlow]:((1-aSlow)*bucket[`sMSlow])+aSlow*bucket[`softMargin];
    bucket[`sMFast]:((1-aFast)*bucket[`sMFast])+aFast*bucket[`softMargin];
    // check if continue
    if[bucket[`counter]>=bucket[`maxIters];
        bucket[`continue]:0;
        bucket[`stopReason]:`nIterExceeded
    ];
    // check if continue using decaying soft margin
    if[bucket[`counter]>=min[(bucket[`sMAlpSlow];bucket[`sMAlpFast])]*count bucket[`data];
        if[bucket[`softMarginDecayError]>(bucket[`sMSlow]-bucket[`sMFast]);
            bucket[`continue]:0;
            bucket[`stopReason]:`softMarginConvergence    
        ]  
    ];
    // update counter
    bucket[`counter]:1+bucket[`counter];
    :bucket; 
 };

// cross-validation of lambda, n-folds
.quantQ.svm.cvFit:{[data;bucket;nCV;lambda]
    // data -- input data with xVar and yVar
    // bucket -- dictionary of parameters
    // nCV -- number of folds
    // lambda -- regularisation parameter
    // set the provided lambda
    bucket[`lambda]:lambda; 
    // create cv flags
    foldFlag:not (floor til[count data]%count[data]%nCV)=/:til nCV;
    // cv: fits
    :{[data;bucket;cv]
        dataIn: data where cv;
        dataOut: data where not cv;
        // fit in-sample model
        bucketIn: .quantQ.svm.fitSoftMargin[dataIn;bucket];
        // return soft-margin on the out-sof-sample
        :avg .quantQ.svm.softMargin[;;bucketIn[`b];bucketIn[`w];bucketIn[`lambda]]'
        [exec yVar from dataOut;exec xVar from dataOut]
    }[data;bucket;] each foldFlag;
 };

// Evaluate SVM model on a dataset
.quantQ.svm.evalModel:{[bucket;data]
    // bucket -- dictionary from estimation
    // data -- data set to evaluate model
    // example: dataModel:.quantQ.svm.evalModel[bucket;data];
    :update yVarPred: ?[yVarPredRaw>=0;1;-1] from update yVarPredRaw:(sum bucket[`w]*flip xVar)-bucket[`b] from data;
 };

// Return detail statistics of the contingency table
.quantQ.svm.evalStats:{[dataModel]
    // dataModel -- data with yVar, yVarPred (-1/1 variables)   
    contingencyTable: select state, cnt from ([] yVar:(-1 -1 1 1);yVarPred:-1 1 -1 1; cnt:4#0; state:(`tn`fp`fn`tp) )  
        pj select cnt:count i by "j"$yVar, "j"$yVarPred from dataModel;
    cDict: {x[`state]!x[`cnt]} contingencyTable;
    // output dictionary
    out:(`contingencyTable`truePositive`trueNegative`falsePositive`falseNegative)!
        (contingencyTable;cDict[`tp];cDict[`tn];cDict[`fp];cDict[`fn]);
    // typeIError/typeIIError
    out:out,(`typeIError`typeIIError)!(cDict[`fp];cDict[`fn]);
    // accuracy 
    out:out,enlist[`accuracy]!enlist[(cDict[`tp]+cDict[`tn])%sum value cDict];
    // sensitivity
    out:out,enlist[`sensitivity]!enlist[cDict[`tp]%cDict[`tp]+cDict[`fn]];
    // specificity
    out:out,enlist[`specificity]!enlist[cDict[`tn]%cDict[`tn]+cDict[`fp]];
    // precision
    out:out,enlist[`precision]!enlist[cDict[`tp]%cDict[`tp]+cDict[`fp]];
    // negPredictiveValue
    out:out,enlist[`negPredictiveValue]!enlist[cDict[`tn]%cDict[`tn]+cDict[`fn]];
    // missRate
    out:out,enlist[`missRate]!enlist[cDict[`fn]%cDict[`fn]+cDict[`tp]];
    // fallout
    out:out,enlist[`fallout]!enlist[cDict[`fp]%cDict[`fp]+cDict[`tn]];
    // falseDiscoveryRate
    out:out,enlist[`falseDiscoveryRate]!enlist[cDict[`tp]%cDict[`fp]+cDict[`tp]];
    // falseOmissionRate
    out:out,enlist[`falseOmissionRate]!enlist[cDict[`tn]%cDict[`fn]+cDict[`tn]];
    // prevalenceThreshold
    out:out,enlist[`prevalenceThreshold]!enlist[
        (
        neg[1]+(cDict[`tn]%cDict[`tn]+cDict[`fp])+sqrt[
        (cDict[`tp]%cDict[`tp]+cDict[`fn])*1-(cDict[`tn]%cDict[`tn]+cDict[`fp])]
        )%
        ((cDict[`tp]%cDict[`tp]+cDict[`fn])+(cDict[`tn]%cDict[`tn]+cDict[`fp])-1)
        ];
    // threatScore
    out:out,enlist[`threatScore]!enlist[cDict[`tp]%cDict[`tp]+cDict[`fn]+cDict[`fp]];
    // balancedAccuracy
    out:out,enlist[`balancedAccuracy]!enlist[0.5*(cDict[`tn]%cDict[`tn]+cDict[`fp])
        +(cDict[`tp]%cDict[`tp]+cDict[`fn])];
    // F1
    out:out,enlist[`F1]!enlist[2*cDict[`tp]%(cDict[`fn]+cDict[`fp]+2*cDict[`tp])];
    // MatthewsCorrelation
    out:out,enlist[`MatthewsCorrelation]!enlist[
        ((cDict[`tp]*cDict[`tn])-(cDict[`fp]*cDict[`fn]))%
        sqrt[(cDict[`tp]+cDict[`fp])*(cDict[`tp]+cDict[`fn])*
        (cDict[`tn]+cDict[`fp])*(cDict[`tn]+cDict[`fn])]];
    // FowlkesMallowsIndex
    out:out,enlist[`FowlkesMallowsIndex]!enlist[sqrt 
        (cDict[`tp]%(cDict[`tp]+cDict[`fp]))*(cDict[`tp]%(cDict[`tp]+cDict[`fn]))];
    // informedness
    out:out,enlist[`informedness]!enlist[neg[1]+(cDict[`tp]%cDict[`tp]+cDict[`fn])
        +(cDict[`tn]%cDict[`tn]+cDict[`fp])];
    // markedness
    out:out,enlist[`markedness]!enlist[neg[1]+(cDict[`tp]%cDict[`tp]+cDict[`fp])
        +(cDict[`tn]%cDict[`tn]+cDict[`fn])];
    // return the full statistics
    :out;  
 };    


/////////////////////////////////////////////
//
/////////////////////////////////////////////
// examples




// // example -- create synthetic data
// data:update y:1+2*{-1*x=0} (x1+x2)>1 from ([] x1:100?1.0; x2:100?1.0);
// data: select yVar:y,xVar:(x1,'x2) from data;
// // run default model
// bucket:()!();
// bucketFit:.quantQ.svm.fitSoftMargin[data;bucket];
// // evaluate fitted model
// dataModel:.quantQ.svm.evalModel[bucketFit;data];
// // evaluate statistics
// .quantQ.svm.evalStats[dataModel]
// // CV the regularisation parameter
// lambdaArray:0.0001 0.01 0.1 1.0;
// nCV:6;
// bucket:()!();
// {([] lambda:x; softMarginCV:y)}[lambdaArray;] avg each .quantQ.svm.cvFit[data;bucket;nCV;] each lambdaArray
// // change the parameters
// bucket:()!();
// // power-law decay
// bucket:bucket,(`learningRateFunc`learningRateParams)!(
//         .quantQ.nn.learningFuncs;(`method`learningRate`power!(`powerLaw;0.01;-0.5)));
// // regularisation lambda
// bucket:bucket,enlist[`lambda]!enlist[0.2];        
// // maximum iterations
// bucket:bucket,enlist[`maxIters]!enlist[30000];        
// // slow ewma propto 0.9 dataset size, fast ewma propto 0.1 dataset size
// bucket:bucket,(`sMAlpSlow`sMAlpFast)!(0.9;0.1);
// // do not use soft margin
// bucket:bucket,enlist[`softMarginDecayError]!enlist[-0wf];
// // monitor convergence
// bucket:bucket,enlist[`monitorConvergence]!(enlist[1]);
// bucketFit:.quantQ.svm.fitSoftMargin[data;bucket];
// // visualise convergence, take every 20-th point
// select i, softMargin from bucketFit where (i%20)=floor (i%20), i>0
