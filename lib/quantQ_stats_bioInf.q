// Needleman–Wunsch algorithm
// one-step of NW algorithm
.quantQ.stats.NeedlemanWunschOneStep:{[bucket]
    // decompose bucket
    tab:bucket[`tab];
    iter:bucket[`iter];
    seq1:bucket[`seq1];
    seq2:bucket[`seq2];
    // index: seq2-row/seq1-column (iRow;iCol)
    iRow: floor iter%count flip tab;
    iCol: mod[iter;count flip tab];
    // initialise first bucket
    if[(iRow=0)and(iCol=0);
        tab[iRow;iCol;0]:0;
    ];
    if[(iRow=0)and(iCol>0);
        tab[iRow;iCol;0]:neg iCol;
        tab[iRow;iCol;1]:(tab[iRow;iCol;1],`left);
    ];
    if[(iRow>0)and(iCol=0);
        tab[iRow;iCol;0]:neg iRow;
        tab[iRow;iCol;1]:(tab[iRow;iCol;1],`up);
    ];
    if[(iRow>0)and(iCol>0);
        // up
        up:tab[iRow-1;iCol;0]-1;
        // leftup
        leftup:tab[iRow-1;iCol-1;0]+$[seq1[iCol-1]=seq2[iRow-1];1;neg[1]];
        // left
        left:tab[iRow;iCol-1;0]-1;
        // decide
        tab[iRow;iCol;1]:(`up`leftup`left) where t=m:max t:(up;leftup;left);
        tab[iRow;iCol;0]:m;   
    ];
    // compose back bucket and return it
    :bucket,(`tab`iter)!(tab;iter+1);
 };

// extract paths from the NW grid
.quantQ.stats.NeedlemanWunschPathsOneRun:{[tab;x] 
    // tab -- NW solution
    // x -- current solution to be extended
    : raze {[tab;yIn]
        y: last yIn;
        outt:enlist yIn;
        if[not (first[y]=0) and (last[y]=0);
            availableDirection: tab[first y;last y;1];
            // (),enlist; enlist availableDirection:`left`up
            newEndPoints: {:$[1=count x; y;y];}[availableDirection;]   {[x;dir] 
                if[`leftup=dir;out:(first[x]-1;last[x]-1)];
                if[`up=dir;out:(first[x]-1;last[x])];
                if[`left=dir;out:(first[x];last[x]-1)];
                :enlist out;}[y;] each availableDirection;
            // raze each {(x;y)}[yIn;] each newEndPoints
            // cross new end points with In
            outt:raze each {(x;y)}[yIn;] each newEndPoints
        ];
        :outt;
    }[tab;] each x;
 };

.quantQ.stats.NeedlemanWunschVisualisePaths:{[seq1;seq2;direction]
    // seq1, seq2 -- two arrays of the same length to be compared
    // directions -- set of direction composing one path
    seqIn1: last flip direction;                
    seqIn2: first flip direction;                
    seq1sym: seq1@ neg[1]+ 1_ reverse[seqIn1]*{$[x<1;0N;x]  } each deltas reverse seqIn1;                                           
    seq2sym: seq2@ neg[1]+ 1_ reverse[seqIn2]*{$[x<1;0N;x]  } each deltas reverse seqIn2;        
    // return the two sequences, missing is pass
    :(seq1sym;seq2sym);
 };

// Needleman-Wunsch algorithm
.quantQ.stats.NeedlemanWunsch:{[seq1;seq2]
    // seq1, seq2 -- two arrays of the same length to be compared
    Nseq1: count seq1;
    Nseq2: count seq2;
    // empty grid
    NWgrid:((1+Nseq2);(1+Nseq1))#((1+Nseq1)*(1+Nseq2))#enlist (0nf;`$());
    // create bucket to iterate through NW
    bucket:((`tab`iter`seq1`seq2)!(NWgrid;0;seq1;seq2));
    // iterate NW and find solution
    bucketSolution:(.quantQ.stats.NeedlemanWunschOneStep)/[(1+Nseq1)*(1+Nseq2);((`tab`iter`seq1`seq2)!(NWgrid;0;seq1;seq2))];
    // initiate array of paths
    directionsList:(enlist (enlist (Nseq1;Nseq2)));
    // get all the paths from the solution
    directionsList:(.quantQ.stats.NeedlemanWunschPathsOneRun[bucketSolution[`tab];])/[directionsList] ;
    // visualise matches                  
    paths:.quantQ.stats.NeedlemanWunschVisualisePaths[seq1;seq2;] each directionsList;          
    // return solutions
    :(`NWgrid`matches)!(bucketSolution[`tab];paths);              
 };

// example 1
// two input sequences
// seq1:`G`C`A`T`G`C`U; seq2:`G`A`T`T`A`C`A;
// .quantQ.stats.NeedlemanWunsch[seq1;seq2]

// example 1
// two input sequences
// seq1:`A`B`C; seq2:`A`B`C;
// .quantQ.stats.NeedlemanWunsch[seq1;seq2]

// example 1
// two input sequences
// seq1:enlist `A; seq2:enlist `A;
// .quantQ.stats.NeedlemanWunsch[seq1;seq2]

