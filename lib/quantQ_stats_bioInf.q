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
    bucketSolution:(.quantQ.stats.NeedlemanWunschOneStep)/[(1+Nseq1)*(1+Nseq2);bucket];
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

// Smith–Waterman algorithm -- local alignment
// one-step of NW algorithm
.quantQ.stats.SmithWatermanOneStep:{[bucket]
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
        tab[iRow;iCol;0]:0.0;
    ];
    if[(iRow=0)and(iCol>0);
        tab[iRow;iCol;0]:0.0|neg iCol*bucket[`indel];
        tab[iRow;iCol;1]:(tab[iRow;iCol;1],`left);
    ];
    if[(iRow>0)and(iCol=0);
        tab[iRow;iCol;0]:0.0|neg iRow*bucket[`indel];
        tab[iRow;iCol;1]:(tab[iRow;iCol;1],`up);
    ];
    if[(iRow>0)and(iCol>0);
        // up
        up:0.0|tab[iRow-1;iCol;0]-bucket[`indel];
        // leftup
        leftup:0.0|tab[iRow-1;iCol-1;0]+$[seq1[iCol-1]=seq2[iRow-1];bucket[`match];neg[bucket[`match]]];
        // left
        left:0.0|tab[iRow;iCol-1;0]-bucket[`indel];
        // decide
        tab[iRow;iCol;1]:(`up`leftup`left) where t=m:max t:(up;leftup;left);
        tab[iRow;iCol;0]:m;   
    ];
    // compose back bucket and return it
    :bucket,(`tab`iter)!(tab;iter+1);
 };
 
// extract paths from the NW grid
.quantQ.stats.SmithWatermanPathsOneRun:{[tab;x] 
    // tab -- SW solution
    // x -- current solution to be extended
    : raze {[tab;yIn]
        y: last yIn;
        outt:enlist yIn;
        if[not[(first[y]=0) and (last[y]=0)] and 0<tab[first y;last y;0];
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

.quantQ.stats.SmithWatermanVisualisePaths:{[seq1;seq2;direction]
    // seq1, seq2 -- two arrays of the same length to be compared
    // directions -- set of direction composing one path
    
    seqIn1: last flip direction;                
    seqIn2: first flip direction;                
    seq1sym: seq1@ neg[1]+ 1_ reverse[seqIn1]*{$[x<1;0N;x]  } each deltas reverse seqIn1;                                           
    seq2sym: seq2@ neg[1]+ 1_ reverse[seqIn2]*{$[x<1;0N;x]  } each deltas reverse seqIn2;        
    // return the two sequences, missing is pass
    :(seq1sym;seq2sym);
 };

// Smith-Waterman algorithm
.quantQ.stats.SmithWaterman:{[params;seq1;seq2]
    // seq1, seq2 -- two arrays of the same length to be compared
    // params -- parameters
    params:((`indel`match)!(1.0;1.0)),params;
    Nseq1: count seq1;
    Nseq2: count seq2;
    // empty grid
    NWgrid:((1+Nseq2);(1+Nseq1))#((1+Nseq1)*(1+Nseq2))#enlist (0nf;`$());
    // create bucket to iterate through SW
    bucket:((`tab`iter`seq1`seq2)!(NWgrid;0;seq1;seq2)),params;
    // iterate SW and find solution
    bucketSolution:(.quantQ.stats.SmithWatermanOneStep)/[(1+Nseq1)*(1+Nseq2);bucket];
    // points with maximum metrics
    directionsList: enlist each raze z where { not 0=count x} each z:{cross[first x;last x]} each flip {(x;y)}[(til (1+Nseq2));] where each ((1+Nseq2);(1+Nseq1))#((t>0)and max[t]=t:first each raze bucketSolution[`tab]);
    // get all the paths from the solution starting at max and ending at zero
    // directionsList: (enlist (0 1);enlist (0 2))
    
    // directionsList:{-1_x} each (.quantQ.stats.SmithWatermanPathsOneRun[bucketSolution[`tab];])/[directionsList];
    directionsList: (.quantQ.stats.SmithWatermanPathsOneRun[bucketSolution[`tab];])/[directionsList];
    
    // visualise matches                  
    paths:.quantQ.stats.SmithWatermanVisualisePaths[seq1;seq2;] each directionsList;          
    // return solutions
    :(`NWgrid`matches)!(bucketSolution[`tab];paths);              
 };

// examples of edge cases
// seq1:`a`t`g`t`t`a`c`g`g;
// seq2:`z`g`g`t`t`g`a`c`t`a;
// .quantQ.stats.SmithWaterman[()!();seq1;seq2]

// seq1:`a`b`c;
// seq2:`d`e`f;
// .quantQ.stats.SmithWaterman[()!();seq1;seq2]

// seq1:`a`a`s`d`d;
// seq2:`b`b`s`v`v;
// .quantQ.stats.SmithWaterman[()!();seq1;seq2]
// .quantQ.stats.SmithWaterman[(`indel`match)!(3;1);seq1;seq2]

// Levenshtein Distance using Wagner-Fischer algorithm
.quantQ.stats.LevenshteinDistance:{[params;seq1;seq2]
    // params -- dictionary with parameters
    // seq1 -- the source sequence which is compared against the reference sequence
    // seq2 -- the reference sequence
    
    // set default parameters
    params:((`delete`insert`substitute)!(1.0;1.0;1.0)),params;
    // Wagner-Fischer algorithm
    Nseq1: count seq1;    
    Nseq2: count seq2;  
    distanceMatrix: ((1+Nseq2);(1+Nseq1))#((1+Nseq1)*(1+Nseq2))#0.0;
    // delete elements of the source 
    distanceMatrix[0;]:params[`delete]* til Nseq1+1;
    // insert elements to the reference 
    distanceMatrix[;0]:params[`delete]* til Nseq2+1;
    // populate the matrix
    iSeq2:1; 
    while[iSeq2<Nseq2+1;
        iSeq1:1;
        while[iSeq1<Nseq1+1;
            cost:params[`substitute];
            if[seq1[iSeq1-1]=seq2[iSeq2-1];
            cost:0
            ];
            distanceMatrix[iSeq2;iSeq1]:min[(
                distanceMatrix[iSeq2-1;iSeq1]+params[`delete];
                distanceMatrix[iSeq2;iSeq1-1]+params[`insert];
                distanceMatrix[iSeq2-1;iSeq1-1]+cost
                                     )];
            iSeq1+:1;    
        ];
        iSeq2+:1;
    ];
    // return the final distanceMatrix and value for comparison of two sequences
    :(`distance`distanceMatrix)!(distanceMatrix[iSeq2-1;iSeq1-1];distanceMatrix);
 };

// example
// seq1:`H`E`L`L`O; seq2:`W`O`R`L`D;
// .quantQ.stats.LevenshteinDistance[()!();seq1;seq2]



