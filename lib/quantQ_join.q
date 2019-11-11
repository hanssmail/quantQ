.quantQ.join.asymmetry:{[x]
    :sum[ x = 1 ] % count x;
 };


.quantQ.join.tStatsCorr:{[rho;N]
    // rho -- correlation
    // N -- size of sample
    :rho*sqrt[(N-2)%1-rho*rho];
 };