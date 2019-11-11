.quantQ.knn.distanceLX:{[point;trainingSampleX;LX]
    // point -- N-dimensional array
    // trainingSampleX -- matrix of training sample points
    // LX -- parameter of the distance measure
    :xexp[;1%LX] sum xexp[;LX] abs point-trainingSampleX;
 };

.quantQ.knn.kNNneighbours:{[point;tableFunction;LX;k]
    // point -- N-dimensional array
    // tableFunction -- table with training data, first N are features, N+1-st is dependent
    // LX -- parameter of the distance measure
    // k -- number of nearest neighbours to be selected
    // returns dictionary with dependent variable (`f) and distance (`dist)
    :flip select[k] from `dist xasc
     ([] f: tableFunction[last cols tableFunction];
      dist: .quantQ.knn.distanceLX[point;tableFunction[cols[tableFunction] except `f];LX]);
 };

.quantQ.knn.kNNclassifyMajorityVote:{[kNN]
    // kNN -- output of kNNneighbours function
    // sort classes by number of votes obtained
    votesTab: 0!`cntF xdesc select cntF:count f by f from flip kNN;
    // return class which was favoured by votes (select randomly if more than one)
    :first 1?exec f from votesTab where cntF=max votesTab[`cntF];
 };

.quantQ.knn.kNNclassifyDistanceWeighted:{[kNN]
    // kNN -- output of kNNneighbours function
    // normalize distances to sum up to one
    normDist:exec sum dist from flip kNN;
    // sort classes by total weight of votes
    votesTab: 0!`cntF xdesc select cntF:sum weight by f from update weight:1-dist%normDist from flip kNN;
    // return class which was favoured by weighted votes (select randomly if more than one)
    :first 1?exec f from votesTab where cntF=max votesTab[`cntF];
 };

.quantQ.knn.prototypeKNN:{[table;k;R;iterations]
    // table -- the training data with array of features x and the classifier (dependent variable) f
    // k -- the size of the neighbourhood
    // R -- the number of centres per each class
    // iterations -- the number of iterations
    // all distinct values
    valuesF: exec distinct f from table;
    // "parallel loop" over all distinct values of dependent variable
    tableGlobalOut:raze {[table;k;R;iterations;fRealization]
        // create a sub-sample
        tableSub: select centre:x from table where f=fRealization;
        // output table, contains initial draws of centres
        tableOut: update centreID:1+i from (neg R)?(select centre:x,f from table where f=fRealization);
        // loop over each centre
        centresNew:{[tableSub;k;iterations;centre]
            // for each centre, iterate search using neighbours
            {[tableSub;k;centreX]
                tableX:update distVect:.quantQ.knn.distanceLX[centreX;;2] flip centre from tableSub;
                // sort table
                tableX:`distVect xasc tableX;
                // choose the neighbourhood, include centre itself
                tableX: tableX[til k+1];
                // calculate the mean of the neighbourhood
                :exec avg centre from tableX;
            }[tableSub;k]/[iterations;centre]
        }[tableSub;k;iterations;] each exec centre from tableOut;
        // update the tableOut by setting the new centres
        update centre:centresNew from tableOut
    }[table;k;R;iterations;] peach valuesF;
    :tableGlobalOut;
 };

.quantQ.knn.getPrototypeKNN:{[tablesCentre;point]
    // tablesCentre -- the table with trained centres, contains variables centre and f
    // point -- the array to classify
    // table with proximity of a point to the centres
    tablesCentreX:`distVect xasc update distVect:.quantQ.knn.distanceLX[point;;2] 
        flip centre from tablesCentre;
    // class of nearest centre
    :first exec f from  tablesCentreX;
 };

.quantQ.knn.weightedL2norm:{[x;y;weight]
    // x, y -- n-dimensional arrrays (two points in feature space)
    // weight -- n-dimensional array of weights 
    :t wsum t:weight*x-y;
 };

.quantQ.knn.kNNMeanFunction:{[coordX;coordY;valueY;weights;beta]
    // coordX -- n-dimensional array (point to make regression)
    // coordY -- kNN x n array, k observations of n-dimensional arrays
    // valueY -- kNN-dimensional array of dependent variables    
    // weights -- n-dimensional array of weights 
    // beta -- parameter of the algorithm
    // returns MSE-based score
    :(1%(sum exp neg (1%beta)*.quantQ.knn.weightedL2norm[coordX;;weights] each coordY))*sum
         valueY*exp neg (1%beta) *.quantQ.knn.weightedL2norm[coordX;;weights] each coordY;
 };

.quantQ.knn.errorLocal:{[coordX;coordY;valueY;weights;beta]
    // coordX -- n-dimensional array (point to make regression)
    // coordY -- kNN x n array, k observations of n-dimensional arrays
    // valueY -- kNN-dimensional array of dependent variables 
    // weights -- n-dimensional array of weights 
    // beta -- parameter of the algorithm
    // returns error at local point
    :neg 0.5*t*t:(.quantQ.knn.kNNMeanFunction[coordX;coordY;valueY;weights;beta]-
        .quantQ.knn.kNNMeanFunction[coordX;coordY;valueY;(count weight)#1;beta]);
 };

.quantQ.knn.aFunc:{[coordX;coordXb;coordXbb;beta;weights]
    // coordX, coordXb, and coordXbb -- n-dimensional arrays (three points)
    // weights -- n-dimensional array of weights 
    // beta -- parameter of the algorithm
    // returns aFunc
    :exp neg ((.quantQ.knn.weightedL2norm[coordX;coordXb;weights]+
        .quantQ.knn.weightedL2norm[coordX;coordXbb;weights])%beta);
 };

.quantQ.knn.uFunc:{[coordX;coordXb;coordXbb;weights]
    // coordX, coordXb, and coordXbb -- n-dimensional arrays (three points)
    // weight -- n-dimensional array of weights 
    // returns uFunc
    :weights*(t*t:(coordX-coordXb))+(u*u:(coordX-coordXbb))
 };

.quantQ.knn.gradFunc:{[coordX;coordY;valueY;weights;beta]
    // coordX -- n-dimensional array
    // coordY -- kNN x n array, k observations of n-dimensional arrays
    // valueY -- kNN-dimensional array of dependent variables
    // weights -- n-dimensional array of weights
    // beta -- parameter of the algorithm
    // returns gradient of a function
    :(1%sum {[coordX;coordY;beta;weights;x]
             .quantQ.knn.aFunc[coordX;coordY[x[0]];coordY[x[1]];beta;weights]
            }[coordX;coordY;beta;weights;] each
            (t cross t:til count coordY))*neg (4%beta)*
          sum {[valueY;coordX;coordY;beta;weights;x]
               valueY[x[1]]*
               .quantQ.knn.aFunc[coordX;coordY[x[0]];coordY[x[1]];beta;weights]*
               .quantQ.knn.uFunc[coordX;coordY[x[0]];coordY[x[1]];weights]
               }[valueY;coordX;coordY;beta;weights;] each
               (t cross t:til count coordY);
 };

.quantQ.knn.gradError:{[coordFull;valueFull;weights;beta;fractionS;kNN]
    // coordFull -- k x n array, k observations of n-dimensional arrays
    // valueFull -- k-dimensional array of dependent variables 
    // fractionS -- relative size of the data used in assessing the score
    // kNN -- number of nearest neighbours used for regression around each point
    // weight -- n-dimensional array of weights 
    // beta -- parameter of the algorithm
    dim: count coordFull;
    setS: asc(neg ceiling[fractionS*dim])?dim;
    // returns gradient of error at given weights
    :neg sum{[coordFull;valueFull;weights;beta;i;kNN]
        coordX:coordFull[i];
        valueX:valueFull[i];
        dim: count coordFull;
        vecDist: .quantQ.knn.weightedL2norm[coordFull[i];;weights] each coordFull;
        neighbourhood: where vecDist<=(first (kNN-1)_(asc vecDist));
        x1:coordFull[neighbourhood];
        y1:valueFull[neighbourhood];
        :.quantQ.knn.gradFunc[coordX;x1;y1;weights;beta]*(valueX-
            .quantQ.knn.kNNMeanFunction[coordX;x1;y1;weights;beta]);
    }[coordFull;valueFull;weights;beta;;kNN] each setS;
 };

.quantQ.knn.kNNFeatureSelection:{[bucket]
    // bucket -- dictionary with variables of the model
    // bucket[`coordFull] -- k x n array, k observations of n-dimensional arrays
    // bucket[`valueFull] -- k-dimensional array of dependent variables 
    // bucket[`learningRatio] -- learning ratio
    // bucket[`decay] -- parameter for decaying learning ratio
    // bucket[`kNN] -- number of nearest neighbours used for regression around each point
    // bucket[`beta] -- parameter of the algorithm
    // bucket[`fractionS] -- relative size of the data used in each iteration of the algorithm
    // bucket[`nIterations] -- number of iterations
    nFeatures: count flip bucket[`coordFull];
    // initiate weights with 1's
    weightInit:nFeatures#1.0;
    // one-iteration function, defined internally
    kNNUpdateWeights:{[coordFull;valueFull;learningRatio;decay;kNN;beta;fractionS;bucket]
        :`weight`cnt!(bucket[`weight]+learningRatio*(decay xexp bucket[`cnt])*
            .quantQ.knn.gradError[coordFull;valueFull;bucket[`weight];beta;fractionS;kNN];bucket[`cnt]+1);
    };
    // loop to iterate the algorithm nIterations times, return outcome
    :(kNNUpdateWeights[bucket[`coordFull];bucket[`valueFull];bucket[`learningRatio];
        bucket[`decay];bucket[`kNN];bucket[`beta];bucket[`fractionS];]\)[bucket[`nIterations];
        `weight`cnt!(weightInit;0)];
 };

