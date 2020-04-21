
// Fisher Exact Test
.quantQ.stats.FisherExactTest:{[tabData]
    // tabData: ([] grp: `a`a`b`b; feature: `A`B`A`B  ;val: 3 4 5 6 )
    // tabData:([] grp: `a`a`b`b; feature: `A`B`A`B  ;val: 3 0 0 3 );
    // .quantQ.stats.FisherExactTest[tabData]         
    f:.quantQ.stats.factorial;
    
    dd:(`a`c`b`d)!first ((2#cols tabData) xasc tabData)[-1#cols tabData];
    // Fisher p-value 
    pValue: f[dd[`a]+dd[`b]]*f[dd[`c]+dd[`d]]*f[dd[`a]+dd[`c]]*f[dd[`b]+dd[`d]]%(f[dd[`a]+dd[`b]+dd[`c]+dd[`d]]*f[dd[`a]]*f[dd[`b]]*f[dd[`c]]*f[dd[`d]]);
    :pValue;
 };
 
 
// comb utility function
// array m x n
// start populating offdiagonal region
.quantQ.stats.offDiagRegion:{[m;n;prop]
    // m:3;n:3;prop:0.01;
    // nsteps top explore
    nSteps:2*max[(m;n)];
    coordinates: distinct raze { til[1+first[x]] cross til[1+last[x]]} each distinct {[m;n;prop;x] (min[((m-1);floor m*(prop-x))];min[((n-1);floor n*x)])}[m;n;prop;] each prop*(1%nSteps-1)*til nSteps-1;

    set1:flip (m;0)+(-1 1)* flip coordinates;
    set2:flip (0;n)+(1 -1)* flip coordinates;
    // union of sets, no diagonal is removed
    :distinct raze[(set1;set2)];
 };
 
 // Barnard test weight
 .quantQ.stats.BarnardTestWeight:{[mn;ab;p]
    // wrap factorial
    f:.quantQ.stats.factorial;
    m: first mn;
    n: last mn;
    a: first ab;
    b: last ab;
    :(f[m]*f[n]%f[a]*f[m-a]*f[b]*f[n-b])*xexp[p;a+b]*xexp[1.0-p;m+n-(a+b)];
  }; 
 
 // calculating the region where we can reject the null that both grps shows the 
 // feature at the same probability using Barnard test 
 .quantQ.stats.BarnardRejectionRegion:{[bucket] 
    // bucket -- input parameters
    // example: bucket:(`m`n`xi`dxi`cl`calcCl`pGrid`run)!(3;3;0.0;0.005;0.1;0.0;200;1b);
    // .quantQ.stats.BarnardRejectionRegion bucket
    // the test is taken from https://www.nature.com/articles/156177a0.pdf
    bucketOut:{[d]
        // increase xi
        xi:min[(1;d[`xi]+d[`dxi])];
        // calculate max of weights for given R over all p's
        clCalc: max 
            {[m;n;xi;p]
                sum .quantQ.stats.BarnardTestWeight[(m;n);;p] each .quantQ.stats.offDiagRegion[m;n;xi]
            }[d[`m];d[`n];d[`xi];] each (1.0%d[`pGrid])*til 1+d[`pGrid];    
        // check if cl was exceeded, or xi reached 1
        $[(xi>=1) or clCalc>=d[`cl];d[`run]:0b;[d[`xi]:xi;d[`calcCl]:clCalc]];
        // return back bucket
        :d;
    }/[{x[`run]};bucket];
    :bucketOut,enlist[`rejectionRegion]!enlist[.quantQ.stats.offDiagRegion[bucket[`m];bucket[`n];bucket[`xi]]];    
 };
 
 
 
 