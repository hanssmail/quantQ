// TDA library

// upper triangular sorted triplets
.quantQ.tda.upperTriangTri:{[loop]
    // loop -- array of numbers
    indices: raze {[x;y] {(y,x)}[;y] each (last[y]+1+til x-1+last[y] )}[count[loop];] 
        each raze {(x,'(x+1+til y-x+1))}[;t] each til neg[1]+t:count[loop];           
    :@[loop;indices];
 };

// distance between pairs, pramatetrised by metric
.quantQ.tda.distance:{[f;x1;x2]
    // f -- distance function
    // x1,x2 -- points (arrays of the same length)
    : f[x1;x2];
 };

// LP norm
.quantQ.tda.Lp:{[p;x1;x2]
    // p -- paratemer  
    // x1,x2 -- points (arrays of the same length)  
    out:0nf;
    if[p=1;out:sum abs x1-x2];
    if[p=0wf;out: max abs x1-x2];
    if[(p>1) and p<0wf; out: xexp[;1.0%p] sum xexp[;p] each abs x1-x2];
    // return norm
    :out;    
 };
// exa:  p:1.2;x1: (1 2 3f); x2:(3 -2 -10f);
// .quantQ.tda.Lp[0wf;x1;x2]

// Get all edges with distances (this table can be large!)
.quantQ.tda.getAllEdges:{[tab;metric]
    // tab -- table with i, point columns
    // metric -- metric function of x1 and x2
    // distance
    dist: raze {[metric;tab;n] {[metric;tab;n;j] metric[tab[`point][n];tab[`point][j]]}[metric;tab;n;]   
        each (n+1)+til count[tab]-n+1   }[metric;tab;] peach til count tab;
    // indices
    ind: raze {[tab;n] {[tab;n;j] (tab[`index][n];tab[`index][j])}[tab;n;]   each (n+1)+til count[tab]-n+1   }[tab;] each til count tab;
    // return tabular form
    :([] edge: ind; distance: dist);
 };
// exa: metric: .quantQ.tda.distance[.quantQ.tda.Lp[2]]
// tabEdges:.quantQ.tda.getAllEdges[tab;metric]


// find all triangles with distance<=thr
.quantQ.tda.getTriangles:{[tabEdges;thr]
    // tabEdges -- outcome of .quantQ.tda.getAllEdges
    // thr -- threshold on edge
    // filter acceptable egdes
    tabEdges: select from tabEdges where distance<=thr;
    
    // find and remove isolated edges (without connection)
    isolatedEdges: {exec t from x where cnt=1} 
        select cnt:count t by t from ([] t:asc raze exec edge from tabEdges);
    tabEdges: delete from tabEdges where {[x;y] max x in y  }[isolatedEdges;] each edge; 
    // triangles
    triangles: {x[`triangles]} ({[bucket]
        // get first edge
        edgeX: first exec edge from bucket[`tab]; 
        // delete first edge and update
        bucket[`tab]:delete from bucket[`tab] where i=0;
        // find all edges (a;.)/(x;.) and (b;.)/(.;b)
        aX: t where not first[edgeX]=t:raze exec edge from bucket[`tab] where {[x;y] max x in y  }[first[edgeX];] each edge;
        bX: t where not last[edgeX]=t:raze exec edge from bucket[`tab] where {[x;y] max x in y  }[last[edgeX];] each edge;
        // intersect of aX/bX will contribute to triangles
        addT:();
        if[0<count aX inter bX; addT: {[x;y] (first[x];last[x];y) }[edgeX;] each  aX inter bX];  
        // update triangles
        newTriangles:bucket[`triangles],addT; 
        :(`tab`cnt`triangles)!(bucket[`tab];bucket[`cnt]+1;newTriangles);}/)  
    [count[tabEdges]-1;(`tab`cnt`triangles)!(tabEdges;0;())];  
    // return all the triangles
    :triangles;
 };
// exa: arrayTriangles:.quantQ.tda.getTriangles[tabEdges;thr]
   
    
// get all VR complexes out of the provided   
.quantQ.tda.getAllVR:{[arrayTriangles]
    // arrayTriangles -- nom-empty array of all triangles
    // arrayEdges -- array of all edges (with given threshold)    
    // infer arrayEdges from triangles
    arrayEdges: distinct raze { @[x;(0 1;0 2;1 2)]} each arrayTriangles;
    // return set of all VR complexes
    :{x[`finalComplexes]}({[bucket]
        // candidate simplex 
        complex: first bucket[`arrayComplexes];    
        // potential nodes
        nodes: {[x;t] t where not t in x }[complex;]  distinct raze  bucket[`arrayEdges];
        // find all nodes, which extend the complex, take first one        
        nodesComplex: first nodes where {[ed;co;n] count[co]=sum (asc each (n,'co)) in ed}
            [bucket[`arrayEdges];complex;] each nodes;
        // update complex
        complexNew: enlist t where not null t:complex,nodesComplex;
        // if empty, remove complex from candidates, add to final list
        if[null nodesComplex;
            bucket[`arrayComplexes]:1_bucket[`arrayComplexes];
            bucket[`finalComplexes]:bucket[`finalComplexes],complexNew
        ];
        // if not empty, add new complex to candidates
        if[not null nodesComplex;
            bucket[`arrayComplexes]:(complexNew,1_bucket[`arrayComplexes]);
            // remove the triangles corresponding to newly added node
            bucket[`arrayComplexes]:bucket[`arrayComplexes] where not bucket[`arrayComplexes] in 
                (@[complex; raze {(x,'(x+1+til y-x+1))}[;t] each til neg[1]+t:count[complex]],'nodesComplex)
        ];
        // decide to continue
        if[0=count[bucket[`arrayComplexes]];
            bucket[`continue]:0;
        ];
        // update bucket
        :(`arrayComplexes`finalComplexes`arrayEdges`continue)!(bucket[`arrayComplexes];bucket[`finalComplexes];bucket[`arrayEdges];bucket[`continue]);
    }/)[{x[`continue]};(`arrayComplexes`finalComplexes`arrayEdges`continue)!(arrayTriangles;();arrayEdges;1)];     
 };      
       
// depends on the order of arrayTriangles, greedy-type algorithm
.quantQ.tda.getAllVRUnique:{[arrayTriangles]
    // arrayTriangles -- nom-empty array of all triangles
    // arrayEdges -- array of all edges (with given threshold)    
    // infer arrayEdges from triangles
    arrayEdges: distinct raze { @[x;(0 1;0 2;1 2)]} each arrayTriangles;
    // return set of all VR complexes
    :{x[`finalComplexes]}({[bucket]
        // candidate simplex 
        complex: first bucket[`arrayComplexes];    
        // potential nodes
        nodes: {[x;t] t where not t in x }[complex;]  distinct raze  bucket[`arrayEdges];
        // find all nodes, which extend the complex, take first one        
        nodesComplex: first u where
            {[complex;aC;x] 
                count[t]=sum (t:@[complex; raze {(x,'(x+1+til y-x+1))}[;t] 
                each til neg[1]+t:count[complex]],'x) in aC}
                    [complex;bucket[`arrayComplexes];] each
                        u:(nodes where {[ed;co;n] count[co]=sum (asc each (n,'co)) in ed}
                        [bucket[`arrayEdges];complex;] each nodes);
        // update complex
        complexNew: enlist t where not null t:complex,nodesComplex;
        // if empty, remove complex from candidates, add to final list
        if[null nodesComplex;
            bucket[`arrayComplexes]:1_bucket[`arrayComplexes];
            bucket[`finalComplexes]:bucket[`finalComplexes],complexNew;
            // remove all the triangles corresponding to newly added node
            bucket[`arrayComplexes]:bucket[`arrayComplexes] where not bucket[`arrayComplexes] in
                @[complex; raze {[x;y] {(y,x)}[;y] each (last[y]+1+til x-1+last[y] )}[count[complex];] each raze {(x,'(x+1+til y-x+1))}[;t] each til neg[1]+t:count[complex]]    
                // .quantQ.tda.upperTriangTri[complex]
        ];
        // if not empty, add new complex to candidates
        if[not null nodesComplex;
            bucket[`arrayComplexes]:(complexNew,1_bucket[`arrayComplexes]);
        ];
        // decide to continue
        if[0=count[bucket[`arrayComplexes]];
            bucket[`continue]:0;
        ];
        // update bucket
        :(`arrayComplexes`finalComplexes`arrayEdges`continue)!(bucket[`arrayComplexes];bucket[`finalComplexes];bucket[`arrayEdges];bucket[`continue]);
    }/)[{x[`continue]};(`arrayComplexes`finalComplexes`arrayEdges`continue)!(arrayTriangles;();arrayEdges;1)];
 };
 
// get all the closed loops from a provided range of edges, this can be computationally very expensive
.quantQ.tda.getAllLoops:{[arrayEdges]
    // arrayEdges -- array with edges
    // discard all isolated edges and then problem above disappears 
    singleNodes: {[x] {exec nodes from x where cnt=1} select cnt:count nodes by nodes from ([] nodes: x)} 
        raze arrayEdges;
    arrayEdges: arrayEdges where not {2=sum y in x}[singleNodes;] each arrayEdges;
    // find all loops
    outt:{x[`loopsFinal]}({[bucketG]
        // get reference edge
        edgeRef: first bucketG[`arrayEdges];
        // remove reference edge from the list
        bucketG[`arrayEdges]:1_bucketG[`arrayEdges];
        // run iteration and get all loops with reference edge        
        loopsWithEdge:{x[`loopsFinal]}({[bucket]
            bucket[`loopsTilde]:raze {[arrayEdges;edgeTilde]
                // get bX
                bX: last edgeTilde;
                // find all edges with bX
                edges2Add: arrayEdges where   {x in y}[bX;] each arrayEdges;
                // check if the edge2Add is closing smaller loop    
                edges2Add: edges2Add where not {[eT;e2A] 2=sum e2A in (1_eT)}[edgeTilde;] each edges2Add;
                // append to edgeTilde
                if[0=count edges2Add;edgeTilde: enlist edgeTilde];
                if[0<count edges2Add;
                    edgeTilde:{(x,y)}[edgeTilde;] each { y where not x=y   }[bX;] each edges2Add        
                ];
                // return the updated edgeTilde
                :edgeTilde;
            }[arrayEdges;] peach  bucket[`loopsTilde];
            // move loops into final; edgeTilde: first bucket[`loopsTilde]
            indLoop: where {first[x]=last[x]} each bucket[`loopsTilde];
            bucket[`loopsFinal]:bucket[`loopsFinal],bucket[`loopsTilde][indLoop];
            bucket[`loopsTilde]: bucket[`loopsTilde] where not til[count[bucket[`loopsTilde]]] in indLoop;      
            // check and remove elements not changing    bucket[`loopsTilde] in bucket[`loopsTildeLast]
            bucket[`loopsTilde]: bucket[`loopsTilde] where not bucket[`loopsTilde] in bucket[`loopsTildeLast];       
            // update loopsTildeLast here
            bucket[`loopsTildeLast]: bucket[`loopsTilde];    
            // update bucket
            :(`edgeRef`arrayEdges`loopsTilde`loopsTildeLast`loopsFinal)!(bucket[`edgeRef];bucket[`arrayEdges];bucket[`loopsTilde];bucket[`loopsTildeLast];bucket[`loopsFinal])  
        }/)[{0<count x[`loopsTilde]};
            (`edgeRef`arrayEdges`loopsTilde`loopsTildeLast`loopsFinal)!(edgeRef;bucketG[`arrayEdges];enlist edgeRef;();())];
        // update registry with loops
        :(`arrayEdges`loopsFinal)!(bucketG[`arrayEdges];bucketG[`loopsFinal],loopsWithEdge);
    }/)[count[arrayEdges]-1;(`arrayEdges`loopsFinal)!(arrayEdges;())];
    // return all the loops
    :outt;
 };
// remark: edgeRef not needed, remove later

// more efficient version keeping only pure loops
.quantQ.tda.getAllPureLoops:{[tabEdges;thr]
    // tabEdges -- outcome of .quantQ.tda.getAllEdges
    // thr -- threshold on edge
    
    // get all triangles
    arrayTriangles:.quantQ.tda.getTriangles[tabEdges;thr]; 
    // get all single loops
    arrayEdges: exec edge from tabEdges where distance<=thr;
    
    // discard all isolated edges and then problem above disappears 
    singleNodes: {[x] {exec nodes from x where cnt=1} select cnt:count nodes by nodes from ([] nodes: x)} 
        raze arrayEdges;
    arrayEdges: arrayEdges where not {2=sum y in x}[singleNodes;] each arrayEdges;
    // find all loops
    outt:{x[`loopsFinal]}({[bucketG]
        // get reference edge
        edgeRef: first bucketG[`arrayEdges];
        // remove reference edge from the list
        bucketG[`arrayEdges]:1_bucketG[`arrayEdges];
        // run iteration and get all loops with reference edge        
        loopsWithEdge:{x[`loopsFinal]}({[bucket]
            bucket[`loopsTilde]:raze {[arrayEdges;edgeTilde]
                // get bX
                bX: last edgeTilde;
                // find all edges with bX
                edges2Add: arrayEdges where   {x in y}[bX;] each arrayEdges;
                // check if the edge2Add is closing smaller loop    
                edges2Add: edges2Add where not {[eT;e2A] 2=sum e2A in (1_eT)}[edgeTilde;] each edges2Add;
                // append to edgeTilde
                if[0=count edges2Add;edgeTilde: enlist edgeTilde];
                if[0<count edges2Add;
                    edgeTilde:{(x,y)}[edgeTilde;] each { y where not x=y   }[bX;] each edges2Add        
                ];
                // return the updated edgeTilde
                :edgeTilde;
            }[arrayEdges;] peach  bucket[`loopsTilde];
            // end of iteration -- check for paths with triangle
            bucket[`loopsTilde]:{x[`allLoops]}({[bucket]
                // get the first triangle
                theTriangle: first bucket[`arrayTriangles];
                // remove the one from the list
                bucket[`arrayTriangles]:1_bucket[`arrayTriangles];
                // remove loops which contain all edges of at least one triangle
                bucket[`allLoops]: bucket[`allLoops] where not {[tH;loop]
                    :$[3>count distinct loop;0b ;1=max .quantQ.tda.upperTriangTri[asc distinct loop] in enlist tH];
                    }[theTriangle;] peach bucket[`allLoops];
                    :bucket;
                }/)[count bucket[`arrayTriangles];(`arrayTriangles`allLoops)!(bucket[`arrayTriangles];bucket[`loopsTilde])];
            // move loops into final; edgeTilde: first bucket[`loopsTilde]
            indLoop: (where {first[x]=last[x]} each bucket[`loopsTilde]);
            indLoopLength: where (4<count each bucket[`loopsTilde]);
            // closed and at least 5 of length (closed and short are de facto discarded)
            bucket[`loopsFinal]:bucket[`loopsFinal],bucket[`loopsTilde][asc indLoop inter indLoopLength];
            // keep all non-closed paths
            bucket[`loopsTilde]: bucket[`loopsTilde] where not til[count[bucket[`loopsTilde]]] in indLoop; 
            // check and remove elements not changing    bucket[`loopsTilde] in bucket[`loopsTildeLast]
            bucket[`loopsTilde]: bucket[`loopsTilde] where not bucket[`loopsTilde] in bucket[`loopsTildeLast];       
            // update loopsTildeLast here
            bucket[`loopsTildeLast]: bucket[`loopsTilde];    
            // update bucket
            :(`edgeRef`arrayEdges`loopsTilde`loopsTildeLast`loopsFinal`arrayTriangles)!(bucket[`edgeRef];bucket[`arrayEdges];bucket[`loopsTilde];bucket[`loopsTildeLast];bucket[`loopsFinal];bucket[`arrayTriangles])  
        }/)[{0<count x[`loopsTilde]};
            (`edgeRef`arrayEdges`loopsTilde`loopsTildeLast`loopsFinal`arrayTriangles)!(edgeRef;bucketG[`arrayEdges];enlist edgeRef;();();bucketG[`arrayTriangles])];
        // update registry with loops
        :(`arrayEdges`loopsFinal`arrayTriangles)!(bucketG[`arrayEdges];bucketG[`loopsFinal],loopsWithEdge;bucketG[`arrayTriangles]);
    }/)[count[arrayEdges]-1;(`arrayEdges`loopsFinal`arrayTriangles)!(arrayEdges;();arrayTriangles)];
    // return all the loops
    :outt;
 };
// exa: \t .quantQ.tda.getAllPureLoops[tabEdges;0.57]

// get all loops, which do not contain any complex inside, brute force approach
.quantQ.tda.getAllPureLoops_NaiveFilter:{[tabEdges;thr]
    // tabEdges -- outcome of .quantQ.tda.getAllEdges
    // thr -- threshold on edge
    // get all triangles
    arrayTriangles:.quantQ.tda.getTriangles[tabEdges;thr]; 
    // get all single loops
    arrayEdges: exec edge from tabEdges where distance<=thr;
    // get all loops
    allLoops:.quantQ.tda.getAllLoops[arrayEdges];
    // remove 3-node loops (a->b->a) and 4-node loops (3-complex)
    allLoops:allLoops where 4<count each allLoops;
    // remove loops which contain at least one triangle    
    pureLoops:{x[`allLoops]}({[bucket]
        // get the first triangle
        theTriangle: first bucket[`arrayTriangles];
        // remove the one from the list
        bucket[`arrayTriangles]:1_bucket[`arrayTriangles];
        // remove loops which contain all edges of at least one triangle
        bucket[`allLoops]: bucket[`allLoops] where not {[tH;loop]
            :1=max .quantQ.tda.upperTriangTri[asc -1_loop] in enlist tH;
        }[theTriangle;] peach bucket[`allLoops];
        :bucket;
    }/)[count arrayTriangles;(`arrayTriangles`allLoops)!(arrayTriangles;allLoops)];
    // return all the pure loops
    :pureLoops;
 };
// exa \t .quantQ.tda.getAllPureLoops_NaiveFilter[tabEdges;0.57]    
  
// Analyse number of triangles as a function of the threshold
.quantQ.tda.nComplexesHistThr:{[tabEdges;thrList]
    // tabEdges -- outcome of .quantQ.tda.getAllEdges
    // thrList -- list of distance thresholds
    // return table with histograms and threshold
    :raze {[t;thr]
        triangles: .quantQ.tda.getTriangles[t;thr];
        outt:([] orderVR: enlist 0; numberVR: enlist 0; thr: enlist thr);
        if[0<count triangles;  
            outt:0!update thr:thr from {select numberVR: count orderVR by orderVR from ([] orderVR:x )}  count each .quantQ.tda.getAllVR triangles;
            if[0=count outt;
                outt:([]orderVR: enlist 0; numberVR: enlist 0; thr: enlist thr);
            ];
        ];
        :outt;
    }[tabEdges; ] each thrList;
 };     
    
// Persistent diagram Number of points which are single and how many are members of at least one edge 
.quantQ.tda.persistentGraphPoint:{[tabEdges;thrList]
    // tabEdges -- outcome of .quantQ.tda.getAllEdges
    // thrList -- list of distance thresholds
    // number of nodes, which are part of edge (dead single points)
    :{([] threshold:x; deathPoint:y )}[thrList;]  
    {[tE;thr] count distinct raze exec edge from tE where distance<=thr}[tabEdges;] each thrList;
 }; 
     
// Analyse numnber of loops as a function of the threshold
.quantQ.tda.nLoopsThr:{[tabEdges;thrList]
    // tabEdges -- outcome of .quantQ.tda.getAllEdges
    // thrList -- list of distance thresholds
    // return table of threshold vs n
    :{([] threshold: x;loopCount:y)}[thrList;] {[t;thr] count .quantQ.tda.getAllLoops exec edge from t where distance<=thr }[tabEdges;] each thrList;
 };       
     
// Persistence diagram (born and death) of loops (1-homology group)
.quantQ.tda.bornAndDeathLoops:{[tabEdges;thrList]
    // tabEdges -- outcome of .quantQ.tda.getAllEdges
    // thrList -- list of distance thresholds
    // get all pure loops
    allLoops: raze {[tE;thr]
     {[x;y] ([] threshold:count[y]#x; loop: y ) }[thr;] .quantQ.tda.getAllPureLoops[tabEdges;thr] 
    }[tabEdges] each thrList;
    // analyse allLoops
    :update isBornBefore: born<=min[thrList], isDeathAfter:death>=max[thrList] from select born: min threshold, death: max threshold, sizeLoop: neg[1]+count first loop by loop from allLoops;
 };     

// TODO: Cech complex
//////////////////////////////////////////////
//    
//////////////////////////////////////////////    
// Examples 
    
// sample data set
// dataSet: update index: i from ([] point:(20 3)#neg[0.5]+60?1.0);   
    
// metric: .quantQ.tda.distance[.quantQ.tda.Lp[2]];
// tabEdges:.quantQ.tda.getAllEdges[dataSet;metric];

// thr:0.4;  
// arrayTriangles:.quantQ.tda.getTriangles[tabEdges;thr];    
// arrayTriangles

// get all complexes  
// .quantQ.tda.getAllVR[arrayTriangles]   
    
// get all unique complexes  
// .quantQ.tda.getAllVRUnique[arrayTriangles]     
    
// Analyse numnber of loops as a function of threshold 
// .quantQ.tda.nLoopsThr[tabEdges;]  0.2+0.005*til 160

// Focus on size=3
// {select thr, numberVR from x where orderVR=3}
// .quantQ.tda.nComplexesHistThr[tabEdges;] (0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6)

// Peristent diagram -- related to 0-homology (points/edges/triangles...)
//  .quantQ.tda.persistentGraphPoint[tabEdges;] 0.1+0.02*til 20

////////////////////////////////////
// loops

// arrayEdges: exec edge from tabEdges where distance<=thr;
// get all loops
// .quantQ.tda.getAllLoops[arrayEdges]  
        
// get all pure loops 
//  .quantQ.tda.getAllPureLoops[tabEdges;0.45]                               

// naive implementation                                                                       
//  .quantQ.tda.getAllPureLoops_NaiveFilter[tabEdges;0.45]  

// Loops analytics
// thrList: (0.25 0.26 0.27 0.28 0.29 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.50 0.51 0.52);
// loopAnalytics:.quantQ.tda.bornAndDeathLoops[tabEdges;thrList]

// select nBorned: count sizeLoop by born from loopAnalytics

// born and death scatter -- persistence diagram
// select born, death from loopAnalytics

// update sums nBorned from select nBorned: count sizeLoop by born from loopAnalytics where sizeLoop=4

