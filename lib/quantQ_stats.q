.quantQ.stats.getHistogram:{[data;rule]
    // data -- array of data
    // rule -- function which takes data and returns bins
    // calculate bins
    bins: rule data;
    // return histogram table
    :update 0^histogram from
        ([bins:asc til count[bins]] x:bins) uj
        (select histogram: count bins by bins from ([] bins:bins binr data));
 };

.quantQ.stats.expma1:{[lambda;vector]
    / lambda -- memory
    / vector -- data    
    {[x;y;z] (x*y)+z}\[ first vector; 1 - lambda; vector * lambda]
 };

.quantQ.stats.expma:{[lambda;vector]
    // lambda -- memory
    // vector -- data
    :(first vector)(1-lambda)\vector*lambda;
 };

.quantQ.stats.histGrid:{[bins;data]
    // bins -- explicit bins to use
    // data -- array of data (not used)
    // return bins
    :bins;
 };

.quantQ.stats.histWidth:{[h;data]
    // h -- width of bin
    // data -- array of data
    // return bins
    :min[data]+h*til 1+ceiling[(max[data]-min[data])%h];
 };

.quantQ.stats.histBin:{[n;data]
    // n -- number of bins
    // data -- array of data
    // return bins
    :min[data]+((max[data]-min[data])%n)*til 1+n;
 };

.quantQ.stats.histSquareRoot:{[data]
    // data -- array of data
    // number of bins
    n:ceiling sqrt count data;
    // return bins
    :min[data]+((max[data]-min[data])%n)*til 1+n;
 };

.quantQ.stats.histSturge:{[data]
    // data -- array of data
    // number of bins
    n:1+ceiling xlog[2;count data];
    // return bins
    :min[data]+((max[data]-min[data])%n)*til 1+n;
 };

.quantQ.stats.histRice:{[data]
    // data -- array of data
    // number of bins
    n:ceiling 2*xexp[count data;1.0%3.0];
    // return bins
    :min[data]+((max[data]-min[data])%n)*til 1+n;
 };

.quantQ.stats.histScott:{[data]
    // data -- array of data
    // width of bins
    h:ceiling (3.5*dev data)%xexp[count data;1.0%3.0];
    // return bins
    :min[data]+h*til 1+ceiling[(max[data]-min[data])%h];
 };

.quantQ.stats.overviewStats:{[data]
    // data -- array of numerical values
    // number of observations
    K: count data;
    // filter out missing data
    data: data where not null data;
    // number of non-missing data
    KminusN:K-count data;
    // sample average
    mu: avg data;
    // sample standard deviation
    sigma: dev data;
    // sample skewness
    tmp: (data-mu);
    N: count data;
    S: (xexp[N*N-1;0.5]%N-2)*((1.0%"f"$N)*sum tmp*tmp*tmp)%xexp[((1.0%"f"$N)*sum tmp*tmp);1.5];
    // sample kurtosis
    Kurt: (((1.0%"f"$N)*sum tmp*tmp*tmp*tmp)%xexp[((1.0%"f"$N)*sum tmp*tmp);2])-3.0;
    // sample min
    mn: min data;
    // sample max
    mx: max data;
    :([] statistics: ("Sample mean";"Sample standard deviation";"Sample skewness";"Sample kurtosis";"Sample min"; "Sample max";"Number of observations";"Number of missing observations");values: (mu; sigma; S;K; mn; mx; Kurt; KminusN));
 };

.quantQ.stats.pValueTwoSided:{[x]
    // x -- value of test statistics following N(0,1)
    // return corresponding p-value
    :2*(reverse first flip .quantQ.stats.normTable) neg[1]+
        (reverse last flip .quantQ.stats.normTable) binr x;
 };

.quantQ.stats.pValueLeftSided:{[x]
    // x -- value of test statistics following N(0,1)
    // return corresponding p-value for one-sided left hand side test
    :$[x>=0;1-(reverse first flip .quantQ.stats.normTable) neg[1]+
        (reverse last flip .quantQ.stats.normTable) binr x;(reverse first flip .quantQ.stats.normTable) 
        neg[1]+(reverse last flip .quantQ.stats.normTable) binr abs x];
 };

.quantQ.stats.pValueRightSided:{[x]
    // x -- value of test statistics following N(0,1)
    // return corresponding p-value for one-sided right hand side test
    :$[x>=0;(reverse first flip .quantQ.stats.normTable) neg[1]+(reverse last flip .quantQ.stats.normTable) 
        binr x;1-(reverse first flip .quantQ.stats.normTable) neg[1]+(reverse last flip .quantQ.stats.normTable)          
        binr abs x];
 };

.quantQ.stats.tStatsCorr:{[rho;N]
    // rho -- estimated correlation
    // N -- size of sample
    :rho*sqrt[(N-2)%1-rho*rho];
 };

.quantQ.stats.tTestOneSample:{[x1;mean]
    // x1 -- array to be tested
    // mean -- mean
    // return t-statistics
    :(avg[x1]-mean)%(dev[x1]% sqrt count x1);
 };

.quantQ.stats.tTestTwoSample:{[x1;x2]
    // x1,x2 -- two arrays to be tested, unequal length
    // return Welch's statistics
    :(avg[x1]-avg[x2])%sqrt ((s*s:dev[x1])%count[x1])+(z*z:dev[x2])%count[x2];
 };

.quantQ.stats.characteristicsBinomial:{[n;N;p]
    // N -- population size
    // n -- number of "1"
    // p -- probability of "1" in every draw
    // first four moments: mean
    mean: N*p;
    // first four moments: variance
    variance: N*p*1-p;
    // sigma -- convenient
    sigma:sqrt variance;
    // first four moments: skewness
    skewness:(1-2.0*p)%sqrt[N*p*1-p];
    // first four moments: kurtosis
    kurtosis:(1-6.0*p*1-p)%N*p*1-p;
    // the output object
     :(`mean`variance`sigma`skewness`kurtosis)!(mean;variance;sigma;skewness;kurtosis);
 };

.quantQ.stats.pdfBinomial:{[n;N;p]
    // N -- population size
    // n -- number of "1"
    // p -- probability of "1" in every draw
    :.quantQ.stats.coeffBinomial[N;n]*xexp[p;n]*xexp[1-p;N-n]
 };

.quantQ.stats.coeffBinomial:{[a;b]
    // a, b -- integers
    // a choose b calculation
    :.quantQ.stats.factorial[a]%(.quantQ.stats.factorial[b]*.quantQ.stats.factorial[a-b])
 };

.quantQ.stats.factorial:{[n]
    // n  -- the integer input
    // for the purpose of calculation, the numbers are recast into float
    :prd "f"$1 + til n;
 };

.quantQ.stats.pValSignTest:{[n;N]
    // n -- number of positive instances
    // N -- size of sample
    // calculate sets of extreme cases
    leftSample: z where n>=z:til N+1;
    rightSample: z where n<=z:til N+1;
    twoSample:asc distinct $[n<=N%2;(z where n>=z:til N+1),(z where (N-n)<=z:til N+1); 
        (z where (N-n)>=z:til N+1),(z where n<=z:til N+1)];
    // return p-values
    :(`pValueLeft`pValueRight`pValueTwoSided)!(sum .quantQ.stats.pdfBinomial[;N;0.5] 
        each leftSample;sum .quantQ.stats.pdfBinomial[;N;0.5] each rightSample;sum 
        .quantQ.stats.pdfBinomial[;N;0.5] each twoSample);
 };

.quantQ.stats.wilcoxonTest:{[x1;x2]
    // x1,x2 -- arrays of sample 1 and sample 2, respectively, of the same length
    // the test is interpreted as a table
    intTab:([] x1;x2);
    // abs diff and the sign of the difference
    intTab: update absDif: abs[x1-x2], sign: signum[x1-x2] from intTab;
    // remove ties
    intTab: select from intTab where sign <>0;
    // sort
    intTab: `absDif xasc intTab;
    // add order of sorted table
    intTab: update rnk:i from intTab;
    // assign average rank to ties
    intTab: update avgRnk: avg rnk by absDif from intTab;
    // stat element of the stats
    w: exec sum avgRnk*sign from intTab;
    // elements of the test, number of observations
    n: count intTab;
    // elements of the test: variance of the test
    sigma: sqrt (1.0%6.0)*n*(n+1)*(1+2*n);
    // return dictionary with results
    :(`w`n`sigma`wNormalised)!(w;n;sigma;w%sigma);
 };

.quantQ.stats.concordanceRoutine:{[row1;row2]
    // rowJ -- pair of observations (xJ;yJ) with J=1,2
    // explicitly extract x's and y's
    x1: first row1;
    y1: last row1;
    x2: first row2;
    y2: last row2;
    // concordance
    concordance: ((x1>x2) and (y1>y2)) or ((x1<x2) and (y1<y2));
    // discordance
    discordance: ((x1>x2) and (y1<y2)) or ((x1<x2) and (y1>y2));
    // output is triplet
    :(concordance;discordance;not concordance or discordance);
 };

.quantQ.stats.kendallTauRank:{[xS;yS]
    // xS, yS -- arrays of values to compare the rank
    // aggregate concordance statistics
    stats: sum raze {.quantQ.stats.concordanceRoutine/:[y;(1+x?y)_x]}[t] each t: flip(xS;yS);
    // return Kendall's Tau Rank
    :(stats[0]-stats[1])%0.5*count[xS]*count[xS]-1;
 };

.quantQ.stats.somersD:{[yS;xS]
    // xS, yS -- arrays of values to calculate D
    // calculate Somers' D using Kendall's Tau Rank
    :.quantQ.stats.kendallTauRank[xS;yS]%.quantQ.stats.kendallTauRank[xS;xS];
 };

.quantQ.stats.somersDBinaryX:{[yS;xS]
    // xS, yS -- arrays of values to calculate D, xS being binary
    // return modified Somers' D
    :0.5*1.0+.quantQ.stats.somersD[yS;xS];
 };

.quantQ.stats.concordanceRoutineBinaryXIntegerY:{[row1;row2]
    // rowJ -- pair of observations (xJ;yJ) with J=1,2
    // XJ -- being binary
    // YJ -- being integers
    // explicitly x's and y's
    x1: first row1;
    y1: last row1;
    x2: first row2;
    y2: last row2;
    // concordance
    concordance: ((x1>x2) and (y1>y2)) or ((x1<x2) and (y1<y2));
    // discordance
    discordance: ((x1>x2) and (y1<y2)) or ((x1<x2) and (y1>y2));
    // tie
    tie: (x1<>x2) and (y1=y2);
    // output is triplet
    :(concordance;discordance;tie);
 };

.quantQ.stats.somersDBinaryXIntegerY:{[yS;xS]
    // xS, yS -- arrays of values to calculate D
    // xS -- being binary
    // yS being integer
    stats: sum raze {.quantQ.stats.concordanceRoutineBinaryXIntegerY/:[y;(1+x?y)_x]}[t]
        each t:flip (xS;yS);
    // Somers' D
    :("f"$stats[0]-stats[1])%("f"$stats[0]+stats[1]+stats[2]);
 };

.quantQ.stats.kruskallGamma:{[xS;yS]
    // xS, yS -- arrays of values to compare the rank
    // aggregate concordance statistics
    stats: sum raze {.quantQ.stats.concordanceRoutine/:[y;(1+x?y)_x]}[t] each t: flip(xS;yS);
    // return gamma Rank coefficient
    :(stats[0]-stats[1])%stats[0]+stats[1];
 };

.quantQ.stats.tKruskall:{[xS;yS]
    // xS, yS -- arrays of values to compare the rank
    // calculate Kruskall's gamma
    gamma: .quantQ.stats.kruskallGamma[xS;yS];
    // count number of combinations
    stats: sum raze {.quantQ.stats.concordanceRoutine/:[y;(1+x?y)_x]}[t] each t: flip(xS;yS);
    // return t_gamma
    :gamma*sqrt (stats[0]+stats[1])%(count[xS]*count[xS]*1-gamma*gamma);
 };

.quantQ.stats.pBonferroni:{[alpha;m]
    // alpha -- family-wise alpha
    // m -- size of family
    // return Bonferroni alpha
    :alpha%m;
 };

.quantQ.stats.pSidak:{[alpha;m]
    // alpha -- family-wise alpha
    // m -- size of family
    // return Sidak alpha
    :1-xexp[(1-alpha);1.0%m];
 };

.quantQ.stats.pHolm:{[pArray;alpha]
    // alpha -- family-wise alpha
    // pArray -- m-dimensional array of p-values
    // return indices of hypothesis which are rejected
    :exec index from select from (update alphaHolm: alpha%(count[pArray]+1-i) 
        from `pIndividual xasc ([] index: 1+til count[pArray]; pIndividual: pArray)) 
        where pIndividual<=alphaHolm;
 };

.quant.stats.runs:{[x] 0 x\x}"f"$;