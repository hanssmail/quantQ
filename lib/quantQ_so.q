// Stochastic Optimisation -- random search

// spherical to Euclidean coordinates
.quantQ.so.sph2Eucl:{[sphCoordinates]
    // sphCoordinates -- (r;vecOfPhi)
    // example: .quantQ.so.sph2Eucl[(1.0;(0;0.2;1.3;2.0))] 
    // example: .quantQ.so.sph2Eucl[(1.0;enlist 0)] 
    // last vecOfPhi <0,2*Pi), others <0,Pi) 
    r: first sphCoordinates;
    vecOfPhi: last sphCoordinates;
    // dimensionality of the problem
    nDim: 1+count[vecOfPhi];
    // transform coordinates
    newCoordinates:
        r*{[nDim;vecOfPhi;n] 
            outt:0.0;
            // create the sin/cos multipliers
            $[n=nDim;
                outt:prd sin vecOfPhi;
                    $[n-1>0; 
                    outt:cos[vecOfPhi[n-1]]*prd sin vecOfPhi til n-1;
                    outt:cos[vecOfPhi[n-1]]
                    ]
            ];
        :outt;     
    }[nDim;vecOfPhi;] each 1+til nDim;
    // return Euclidean coordinates,from origin
    :newCoordinates;
 };

// get random point on n-sphere with radius r
.quantQ.so.randNSphere:{[r;nDim]
    // r -- radius of n-sphere
    // nDim -- dimension of the n-sphere
    // example: {wsum[x;x]} .quantQ.so.randNSphere[1.0;3]
    // example: 1={wsum[x;x]} .quantQ.so.randNSphere[1.0;3]
    $[nDim>2;
        coord:.quantQ.so.sph2Eucl[(r;raze[(nDim-2)?2*acos[0];1?4*acos[0]])];
        coord:.quantQ.so.sph2Eucl[(r;1?4*acos[0])]
    ];
    :coord;
 };
 
// find random set of points around the pivot with certain radius 
.quantQ.so.nSphereNeighbours:{[pivotPoint;r;nN]
    // pivotPoint -- centre of the n-sphere
    // r -- radius on n-sphere neighbourhood
    // nN -- size of neighbourhood
    // example: pivotPoint:(2.0;2.0;2.0);r:1.0;nN:5;
    // .quantQ.so.nSphereNeighbours[pivotPoint;r;nN]
    :pivotPoint+/:{[r;pivotPoint;n] .quantQ.so.randNSphere[r;count[pivotPoint]]}[r;pivotPoint;] each til nN;
 };

.quantQ.so.minimise:{[bucket]
    // bucket -- dictionary with parameters    

    // learning function -- power-decay default
    bucket:(enlist[`learningParams]!enlist[(`method`learningRate`power)!(`powerLaw;1.0;neg 0.5)]),bucket;  
    bucket:(enlist[`learningFunc]!enlist[.quantQ.nn.learningFuncs[bucket[`learningParams];]]),bucket;
    // exploring outside of the region    
    bucket:((`pJump`rMult)!(0.05;100)),bucket;
    // convergence parameters
    bucket:((`continue`deltaL2`nL2`nL2counter`maxEpochs`nN`monitorConvergence`counter)!(1;0.01;10;0;1000;50j;1;0)),bucket;
    
    // run the aslgorithm
    $[bucket[`monitorConvergence]=1;
        :(.quantQ.so.oneStep\)[{x[`continue]};bucket];
        :(.quantQ.so.oneStep/)[{x[`continue]};bucket]
    ];
 };

// one random search
.quantQ.so.oneStep:{[bucket]
    // bucket -- dictionary with parameters

    // decide learnign rate (diameter of the n-sphere) -- here add also the jump of random dimension
    learningRate: bucket[`learningFunc][bucket[`counter]];
    if[bucket[`pJump]>=first 1?1.0;learningRate:learningRate*first 1?bucket[`rMult]];    
                
    // create n-neighbours, add current point -- dim>1
    if[1<count bucket[`x0];
    setOfPoints: raze (enlist bucket[`x0];.quantQ.so.nSphereNeighbours[bucket[`x0];learningRate;bucket[`nN]])
    ];
    
    // create n-neighbours, add current point -- dim=1, nN -> 2
    if[1=count bucket[`x0];    
    setOfPoints: raze (enlist bucket[`x0];bucket[`x0] +/: (neg learningRate; learningRate))
    ];
    
    setOfYs: bucket[`func] each setOfPoints;
   
    // the candidate solution, if more, choose randomly
    candidateSolution:first 1?setOfPoints where setOfYs=min[setOfYs];
    // the improvement
    improvement:abs bucket[`func][bucket[`x0]]-bucket[`func][candidateSolution]; 

    // check improvement
    $[improvement<=bucket[`deltaL2];bucket[`nL2counter]+:1;bucket[`nL2counter]:0];

    // stopping rules
    if[bucket[`counter]>=bucket[`maxEpochs];bucket[`continue]:0];
    if[bucket[`nL2counter]>=bucket[`nL2];bucket[`continue]:0];
    // update counter
    bucket[`counter]:bucket[`counter]+1;
    
    // update the solution
    bucket[`x0]:candidateSolution;
  
    // return updated bucket
    :bucket;    
 };



////////////////////////////////////////////////////////////////
// Examples
// bucket explanation
// `func -- real valued function to minimise, function of single argument 
// `x0 -- starting point, vector
// `learningFunc -- learning parameter function[params;counter]
// `learnignParams -- dictionary of parameters of learning function 
// `pJump -- probability of jumping out of regular learning path
// `rMult -- learning rate multiplier, jumping out of the learning path
// `continue -- binary to keep iteration going
// `deltaL2 -- resolution of the function changes between iterations
// `nL2 -- maximum tries without improvement
// `nL2counter -- counter of tries without improvement
// `maxEpochs -- maximum number of epochs
// `nN -- number of points exploring neighbourhood per epoch
// `monitorConvergence -- binary switch to monitor convergence
// `counter -- counter of iterations    

// Example 1
// bucket:(`func`x0)!({ wsum[x;x]};(1.0;2.0;3.0));
// zz:.quantQ.so.minimise[bucket];
// select bucket[`func] each x0 from zz
// select x0 from zz

// Example 2
// bucket:(`func`x0)!({ wsum[x;x]};enlist 45.0);
// zz:.quantQ.so.minimise[bucket];
// select bucket[`func] each x0 from zz
// select x0 from zz



