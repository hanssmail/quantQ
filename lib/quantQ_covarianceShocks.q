// range of eigenvalues of Dth 
.quantQ.rmt.rangeDTh:{[Nin;Tin;Nout;Tout]
    // Nin,Tin -- dimension of the in-data
    // Nout,Tout -- dimension of the out-data
    qIn:Nin%Tin;
    qOut:Nout%Tout;

    lamMax:((1+qIn+qOut*(1-qIn))+2*sqrt[qIn+qOut*(1-qIn)])%z*z:(1-qIn);
    lamMin:((1+qIn+qOut*(1-qIn))-2*sqrt[qIn+qOut*(1-qIn)])%z*z:(1-qIn);

    blur:xexp[Nin*2.7*1e-3;-2%3];
    // return dict 

    // fraction of zero eigenvalues
    fracZero:max[(0.0;1-1.0%qOut)];
    // return data
    :(`lamMax`lamMin`blur`fracZero)!(lamMax;lamMin;blur;fracZero);    
 };

// assess change of covariance matrices of two data sets
.quantQ.rmt.covChange:{[dataIn;dataOut]
    // dataIn -- table with features, in-sample
    // dataOut -- table with features, out-sample
    // dimensions
    Tin: count (flip value flip dataIn);
    Sin: count (value flip dataIn);
    Tout: count (flip value flip dataOut);
    Sout: Sin;
    // empirical cov matrices
    Ein:.quantQ.mat.covEmpTab[dataIn];
    Eout:.quantQ.mat.covEmpTab[dataOut];
    // inv half matrix
    EinNegHalf1: inv .quantQ.mat.powerSeriesSq[()!();Ein];
    EinNegHalf2: inv .quantQ.mat.DenmanBeavers[()!();Ein];
    // EinNegHalf:EinNegHalf2^EinNegHalf1;
    EinNegHalf:EinNegHalf2;
    // D matrix
    Dmat: EinNegHalf mmu Eout mmu EinNegHalf;
    // eigensystem of Dmat
    eigensystem:.quantQ.mat.eigenSystem[Dmat];
    // thresholds
    thresholds:.quantQ.rmt.rangeDTh[Sin;Tin;Sout;Tout];
    // PCA vectors , row by row
    Omat: exec weights from {x[`tab]} .quantQ.mat.pca Ein;
    // Dmat with in-sample projections
    DO:Omat mmu Dmat mmu flip Omat;
    // calculating first fleeting mode
    eigDO:.quantQ.mat.eigenSystem DO;
    fleetingModes:eigDO[`eigenvectors];
    firstFleetingMode:first[fleetingModes];
    // overlap array
    overlapArray: sums z*z:firstFleetingMode mmu flip Omat;
    // output
    :eigenvalues`lambdaMax`lambdaMaxNCor`overlap`inSampleVecs`fleetingVecs)!(eigensystem[`eigenvalues];thresholds[`lamMax];thresholds[`blur];overlapArray;Omat;fleetingModes);
 };

// wrapper to compare empirical cov matrix step by step
.quantQ.rmt.covTabChange:{[bucket;tab]
    // bucket -- parameters of the algorithm
    // tab -- table with data, first column time index, the 2nd to last are features
    bucket:((`Tin`Tout)!(10;10)),bucket;
    // dimensions of the problem
    iStart:bucket[`Tin]-1;
    iEnd:count[tab]-bucket[`Tout]+1;
    iRange: iStart+til 1+iEnd-iStart;
    // assuming enough data, otherwise add check for data size here
    :delete index from (update index:i from tab) lj 1!raze {[tabF;bucket;iP]      
        // tabF -- input table with features
        // bucket -- parameters
        // iP -- pivot index
        ![0N;iP];
        dataIn: select from tabF where i<=iP, i>iP-bucket[`Tin]; 
        dataOut: select from tabF where i>iP, i<=iP+bucket[`Tout];  
        :([] index: enlist iP),'flip {key[x]!(enlist each value[x])} .quantQ.rmt.covChange[dataIn;dataOut];
    }[?[tab;();0b;(1_cols[tab])!(1_cols[tab])];bucket;] peach iRange;
 };

