.quantQ.apriori.makeDictionary:{[x]
    / x -- list of discrete unique variables
    :x!x=/:x
 };

.quantQ.apriori.createDummies:{[largeTab;i]
    // largeTab -- table to convert
    // i -- index of variable to convert
    // get all variables except the first one
    listCols: 1_cols largeTab;
    // get unique values of selected variable
    distTMP: asc distinct largeTab[listCols[i]];
    // create dictionary
    dictTMP: .quantQ.apriori.makeDictionary[distTMP];
    // vector of values to convert 
    tabTMP: largeTab[listCols[i]];
    // dummy variables
    :?[largeTab;();0b;(`$(raze string listCols[i],"_"),/:string til count distTMP)!flip (dictTMP tabTMP)];
 };

.quantQ.apriori.isIn:{[xOUT;xTEST]
    // xOUT -- unique list of values 
    // xTEST -- unique list of values to be compared
    :count[xOUT]=sum sum xOUT=\:xTEST;
 };

.quantQ.apriori.aprioriOneRun:{[bParams]
    // bParams -- dictionary with all data and parameters
    // increase counter
    bParams[`step]+:1;    
    // create set of all combinations -- distinguish first step
    $[bParams[`step]=1;zz:enlist each 1_cols bParams[`largeTabDummies]; zz:asc each t where
        bParams[`step]=count each t: distinct each distinct asc bParams[`colDummiesStepPrev] cross 
        bParams[`colDummiesStep1]];    
    // if step>1: test whether zz does not have subset from colDummiesStepALLX
    $[bParams[`step]>1;zz:first flip t where 1b=last each t:{[bParams;x] (x; not max 
        .quantQ.apriori.isIn[;x] each bParams[`colDummiesStepALLX])}[bParams;] each zz;];
    // numerical criteria
    yy: flip {[largeTabDummies;whr]
        :(whr;sum prd largeTabDummies[whr])
        }[bParams[`largeTabDummies];] each zz;
    // temporary column variables -- split zz into two sub-sets
    colDummiesStepTMP: asc each zz where bParams[`thresh]<(value first[yy]!last[yy])%
        count bParams[`largeTabDummies];  
    colDummiesStepTMPX: asc each zz where not zz in colDummiesStepTMP;
    // define colDummiesStep1 in step 1
    $[bParams[`step]=1;bParams[`colDummiesStep1]:colDummiesStepTMP;];
    // populate bucket
    bParams[`colDummiesStepPrev]:distinct colDummiesStepTMP;
    bParams[`colDummiesStepALL]:distinct bParams[`colDummiesStepALL],distinct colDummiesStepTMP;
    bParams[`colDummiesStepPrevX]:distinct colDummiesStepTMPX;
    bParams[`colDummiesStepALLX]:distinct bParams[`colDummiesStepALLX], distinct colDummiesStepTMPX;
    // check if further improvement can be done
    $[0=count colDummiesStepTMP;bParams[`canImprove]:0b;bParams[`canImprove]:1b];
    // return bucket
    :bParams;
 };
