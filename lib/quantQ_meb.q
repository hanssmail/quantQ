// Maximum Entropy Block Bootstrap for Non-stationary Processes
// Based on https://econpapers.repec.org/article/risactuec/0115.htm


// The original formulation of the Maximum entropy bootstrap
.quantQ.meb.BSorig:{[xT]
    // xT -- array of time series (ordered array)
    // order index
    ixT: iasc xT;
    // order statistics
    xTBar: xT@ixT;

    // trimming
    dtrm: avg abs 1_deltas xTBar;

    // trimmed ends
    xTBarIni:first[xTBar]-dtrm;
    xTBarFin:last[xTBar]+dtrm;
    xTBar:raze (xTBarIni;xTBar;xTBarFin);

    // random draws
    ps:count[xT]?1.0;
    // zT
    zT:(-1_xTBar)+0.5*1_deltas xTBar;
    rT: floor ps*count[xT];
    // new stats
    xTS:{[arr;n] arr[n]+(first 1?1.0)*(arr[n+1]-arr[n])    }[zT;] each rT;
    // order index
    ixTS: iasc xTS;
    // re-order xTS based on xT
    xTNew:asc[xTS]@ixT;
    // return BS array
    :xTNew;
 };
// example  .quantQ.meb.BSorig[(0.23 0.24 0.12 0.09 0.02 0.45)]
 
// The new version of the Maximum Entropy bootstrap, exp tails
.quantQ.meb.BS:{[bucket;xT]
    // xT -- array of time series (ordered array)
    // default parameters, using exponential tails
    bucket:((`lambdaCap`zero)!(100.0;1e-6)),bucket;    
    // order index
    ixT: iasc xT;
    // order statistics
    xTBar: xT@ixT;
    // cap for x1=x2 or xT-1=xT, respectively
    $[((xTBar[1]-xTBar[0])%4.0)<(1.0%bucket[`lambdaCap]);
        lambda1:bucket[`lambdaCap];
        lambda1:4.0%(xTBar[1]-xTBar[0])
    ];
    $[((xTBar[count[xTBar]-1]-xTBar[count[xTBar]-2])%4.0)<(1.0%bucket[`lambdaCap]);
        lambdaT:bucket[`lambdaCap];
        lambdaT:4.0%(xTBar[count[xTBar]-1]-xTBar[count[xTBar]-2])
    ];
    bucket:bucket,((`lambda1`lambdaT)!(lambda1;lambdaT));
    // boundaries set +-infty
    xTBar:raze (-0wf;xTBar;0wf);
    // random draws
    ps:count[xT]?1.0;
    // zT
    zT:(-0wf)^(-1_xTBar)+0.5*1_deltas xTBar;
    rT: floor ps*count[xT];
    // new stats, split case 1, case T, and all other cases
    xTS:{[bucket;arr;n]
        // arr:zT;n:5 
            out:0.0;        
            // case - left hand side
            if[n=0;
            out:(log[bucket[`zero]|first 1?1.0]%bucket[`lambda1])+arr[1]];
            // case - right hand side
            if[n=count[arr]-2;  
            out: (arr[count[arr]-2])- log[bucket[`zero]|first 1?1.0]%bucket[`lambdaT]];
            // case other
            if[not (n=0) or (n=count[arr]-2); 
            out:arr[n]+(first 1?1.0)*(arr[n+1]-arr[n])];    
            :out;
        }[bucket;zT;] each rT;
    // order index
    ixTS: iasc xTS;
    // re-order xTS based on xT
    xTNew:asc[xTS]@ixT;
    // return BS array
    :xTNew;
 };
 
// example: .quantQ.meb.BS[()!();(0.23 0.24 0.12 0.09 0.02 0.45)]

// plain block shuffle
.quantQ.meb.blockShuffle:{[bucket;data]
    // bucket -- dictionary with parameters
    // data -- array of data to shuffle
    bucket:(enlist[`block]!enlist[6]),bucket;
    // decide number of blocks
    nBlock:ceiling count[xT]%bucket[`block];
    // starting points
    nStartingPoints:nBlock?1+count[data]-bucket[`block];
    // blocks
    nBlocks:flip nStartingPoints +\(bucket[`block]#1);
    // return bootstrap
    :count[data]#raze data @ nBlocks;
 };

// example: .quantQ.meb.blockShuffle[()!();100?1.0]

// Plain block shuffle with uniform noise, sewing option
.quantQ.meb.blockShuffleUni:{[bucket;data]
    // bucket -- dictionary with parameters
    // data -- array of data to shuffle
    bucket:((`block`uniSpread`sew)!(6;0.1;1b)),bucket;
    // decide number of blocks
    nBlock:ceiling count[data]%bucket[`block];
    // starting points
    nStartingPoints:nBlock?1+count[data]-bucket[`block];
    // blocks
    nBlocks:flip nStartingPoints +\(bucket[`block]#1);
    // last of first block and the first of the second block
    sewingIndices: flip ((last each nBlocks);next[first each nBlocks]);
    // add noise
    dataBlock:{[b;x] x+neg[0.5*b[`uniSpread]]+count[x]?b[`uniSpread] }[bucket;] each (data @ nBlocks);
    if[bucket[`sew];
        // remove difference
        dataBlock:dataBlock+count[dataBlock]#sums neg raze[(0.0;(last each deltas each flip ((last each dataBlock);next[first each dataBlock])))];
        // add difference
        dataBlock:dataBlock+count[dataBlock]#sums raze (0.0;(last each deltas each data @ sewingIndices));
    ];
    :count[data]#raze dataBlock;
 };

// example: .quantQ.meb.blockShuffleUni[()!();100?1.0]

// Block Maximum entropy bootstrap
.quantQ.meb.blockMeb:{[bucket;data]
    // bucket -- dictionary with parameters
    // data -- array of data to shuffle;data:100?1.0
    bucket:((`block`uniSpread`sew)!(6;0.1;1b)),bucket;
    // decide number of blocks
    nBlock:ceiling count[data]%bucket[`block];
    // starting points
    nStartingPoints:nBlock?1+count[data]-bucket[`block];
    // blocks
    nBlocks:flip nStartingPoints +\(bucket[`block]#1);
    // last of first block and the first of the second block
    sewingIndices: flip ((last each nBlocks);next[first each nBlocks]);
    // add noise
    dataBlock:.quantQ.meb.BS[()!();] each (data @ nBlocks);
    // sewing
    if[bucket[`sew];
        // remove difference
        dataBlock:dataBlock+count[dataBlock]#sums neg raze[(0.0;(last each deltas each flip ((last each dataBlock);next[first each dataBlock])))];
        // add difference
        dataBlock:dataBlock+count[dataBlock]#sums raze (0.0;(last each deltas each data @ sewingIndices));
    ];
    :count[data]#raze dataBlock;
 };

// example: .quantQ.meb.blockMeb[()!();100?1.0]


