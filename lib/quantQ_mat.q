.quantQ.mat.diagMatrix:{[dim]
    // dim -- dimension of diagonal matrix
    :"f"$v=/:v:til dim
 };

.quantQ.mat.diagMatrix2:{[dimM]
    // dimM -- dimension of diagonal matrix
    :{[x;y] t:x#0.0; t[y]:1.0; t}[dimM;] each til dimM
 };

.quantQ.mat.givensRotation:{[dimM;i;j;theta]
    // dimM -- dimension of the matrix
    // i -- dim1 (0 to dimM-1)
    // j -- dim2 (0 to dimM-1)
    // theta -- angle of rotation
    // unit matrix
    givens: .quantQ.mat.diagMatrix[dimM];
    // insert cos
    givens[i;i]:cos[theta];
    givens[j;j]:cos[theta];
    // insert sin
    givens[min[(i;j)];max[(i;j)]]:neg[sin[theta]];
    givens[max[(i;j)];min[(i;j)]]:sin[theta];
    // rotation matrix
    :givens
 };

.quantQ.mat.QRdecomp:{[matt]
    // matt -- input matrix
    // calculate dimension of the matrix
    dimM: count[matt];
    // list of coordinates to iterate over
    coord: raze {[dimM;x]    {[x;y] (y;x)}[x]'[t where x<t:reverse til dimM]}[dimM;] each til dimM;
    // inputmatrix
    inputMat:matt;
    // initialised rotation matrix
    rotMat:.quantQ.mat.diagMatrix[dimM];
    // counter for bucket
    counter:0;
    // bucket with four variables for iterative procedure
    bucket:((`inputMat`rotMat`coord`counter)!(inputMat;rotMat;coord;counter));
    // iteration to find QR decomposition
    bucket:({[bucket]
        // coordinates at the given step
        i:first bucket[`coord][bucket[`counter]];
        j:last bucket[`coord][bucket[`counter]];
        // theta to rotate the matrix
        theta:atan[neg[bucket[`inputMat][i;j]]%bucket[`inputMat][min[(i;j)];min[(i;j)]]];
        // dimension of the problem
        dimM: count bucket[`inputMat];
        // rotation matrix at the given step
        rotTMP:.quantQ.mat.givensRotation[dimM;j;i;theta];
        // update inputMat
        bucket[`inputMat]:rotTMP mmu bucket[`inputMat];
        // update rotMat
        bucket[`rotMat]:rotTMP mmu bucket[`rotMat];  // tohle bylo prehozene
        // increase counter
        bucket[`counter]:bucket[`counter]+1;
        // return bucket
        :bucket
    }/)[count[coord];bucket];
    // output bucket: initial matrix, rotated matrix and rotation matrix
    :((`matrix`matrixRotated`rotationMatrix)!(matt;bucket[`inputMat];bucket[`rotMat]));
 };

.quantQ.mat.iterQR:{[bucket]
    // bucket -- contains matrix and counter
    // QR decomposition
    bucketTMP:.quantQ.mat.QRdecomp[bucket[`matrix]];
    // define Q and R
    Q:flip bucketTMP[`rotationMatrix];
    R: bucketTMP[`matrixRotated];
    // QT:bucket[`rotationMatrix];
    // define new matrix, similar to original one
    bucket[`matrix]:R mmu Q;
    // increase counter
    bucket[`counter]:bucket[`counter]+1;
    // return bucket
    :bucket
 };

.quantQ.mat.algoQR:{[bucket]
    // bucket -- contains matrix, counter, max number of iterations and precision
    // QR decomposition
    bucketTMP:.quantQ.mat.QRdecomp[bucket[`matrix]];
    // define Q and R
    Q:flip bucketTMP[`rotationMatrix];
    R: bucketTMP[`matrixRotated];
    dimM: count bucket[`matrix];
    // define new matrix, similar to original one
    bucket[`matrix]:R mmu Q;
    // increase counter
    bucket[`counter]:bucket[`counter]+1;
    // decide if continue
    // criterion 1: counts
    $[bucket[`counter]>=bucket[`maxCounts];bucket[`continue]:0b;];
    // criterion 2: precision
    $[bucket[`thresholdZero]>max abs {[matt;coord] matt[first[coord];last[coord]]}[bucket[`matrix];] 
        each raze {[dimM;x]    {[x;y] (y;x)}[x]'[t where x<t:reverse til dimM]
        }[dimM;] each til count[bucket[`matrix]];bucket[`continue]:0b;];
    :bucket
 };

.quantQ.mat.extractEigenvalQR:{[bucket]
    // bucket -- outcome from algoQR
    // extract diagonal and sort
    :desc {x[y;y]}[bucket[`matrix];] each til count bucket[`matrix];
 };

.quantQ.mat.normL1:{[matt]
    // matt -- matrix
    // return L1 (column norm)
    :max sum abs matt;
 };

.quantQ.mat.inverseIteration:{[mattInit;eigenValue;toler;maxIter]
    // mattInit -- the matrix
    // eigenValue -- the eigenvalue
    // toler -- numerical tolerance
    // maxIter -- maximum number of iterations
    // dimension of the matrix
    dimM: count mattInit;
    // initial starting vector, random
    vecB:neg[1.0]+dimM?2.0;
    // set the first element of the random vector to be positive
    $[0<=first vecB; ; vecB:neg vecB];
    // prepare the bucket
    bucketI:((`matrix`eigenVal`vecB`tolerance`counter`maxIter`continue)!(mattInit;eigenValue;vecB;toler;0j;maxIter;1b));
    bucketI:({[bucketI]
        // the matrix for inverse iteration
        invMat:(bucketI[`matrix]-bucketI[`eigenVal]*.quantQ.mat.diagMatrix[count[bucketI[`matrix]]]);
        // the inverse iteration method works for non-singular matrices
        $[0<sum{[x] (x=0nf) or (x=0wf) or (x=-0wf)} each raze inv[invMat];
            // introduce random perturbation
            bucketI[`eigenVal]:bucketI[`eigenVal]+(1e-7)*(neg[1.0]+first 1?2.0);
        [
            vecBOld:bucketI[`vecB];
            vecB: (1.0% .quantQ.mat.normL1[t])*t:(inv[invMat] mmu vecBOld);
            // set the first element of the random vector to be positive
            $[0<=first vecB; ; vecB:neg vecB];
            // plug into bucket
            bucketI[`vecB]:vecB;
            // increase counter
            bucketI[`counter]:bucketI[`counter]+1;
            // switch off iteration
            $[((max[abs[vecBOld-vecB]])<bucketI[`tolerance]) or (bucketI[`counter]>bucketI[`maxIter]);bucketI[`continue]:0b; ];
            // in the last step, normalize to L2=1
            $[bucketI[`continue]=0b; bucketI[`vecB]:bucketI[`vecB]*(1.0%sqrt[sum bucketI[`vecB]*bucketI[`vecB]]); ]
        ]
        ];
        // output bucket
        : bucketI
        }/)[{x[`continue]};bucketI];
    // eigenvector
    :bucketI[`vecB];
 };

.quantQ.mat.eigenSystem:{[matt]
    // matt -- input matrix
    // set the bucket, default controls
    bucket:((`matrix`counter`maxCounts`thresholdZero`continue)!(matt;0j;10000j;1e-10;1b));
    // solve for eigenvalues
    bucket:(.quantQ.mat.algoQR/)[{x[`continue]};bucket];
    // extract eigenvalues
    eigenvalues:.quantQ.mat.extractEigenvalQR[bucket];
    // solve for eigenvectors
    eigenvectors: .quantQ.mat.inverseIteration[matt;;1e-10;10000j] each eigenvalues;
    // output bucket
    :((`eigenvalues`eigenvectors)!(eigenvalues;eigenvectors))
 };

.quantQ.mat.covarianceMatrixFast:{[data]
    // data -- array of data
    :(data+flip(not n=\:n)*data:(n#'0.0),'(data$/:'(n:til count data)_\:data)%count first data)
        -a*\:a:avg each data;
 };

.quantQ.mat.pca:{[data]
    // data -- array of input data
    // calculate eigensytem
    pTab:flip .quantQ.mat.eigenSystem[.quantQ.mat.covarianceMatrixFast[data]];
    // normalised data
    dataNorm: data - avg each data;
    // add variance for each eigenvector (eigenvalues)
    pTab: update variance: {var y wsum x}[dataNorm;] each eigenvectors from pTab;
    // caluclate cumulative variance
    pTab: update relVariance: variance%sum variance from pTab;
    // output format
    pTab:select variance, weights:eigenvectors, relVariance from pTab;
    // transformed data
    dataTransformed: wsum[;data] each pTab[`weights];
    // return overview data and transformed data
    :(`tab`dataTransformed)!(pTab;dataTransformed);
 };

.quantQ.mat.minorMatrix:{[mat;i;j]
    // mat -- square matrix
    // i -- row index
    // j -- column index
    :mat[(t where not i=t:til count[mat]);(s where not j=s:til count[mat])];
 };

.quantQ.mat.det:{[mat]
    // mat - square matrix
    // dimension of matrix
    dimM: count mat;
    // split based on the dimension of the problem
    // default output value
    det:0.0;
    // dim=1
    $[dimM=1;det:mat; ];
    // dim=2
    $[dimM=2;det:(mat[0;0]*mat[1;1])-(mat[0;1]*mat[1;0]); ];
    // dim=3, Sarrus rule
    $[dimM=3;
        [
        // extend matrix
        mat:raze (mat ; mat);
        // positive diagonal
        posDiagonal:sum {[mat;dimM;offset]  prd {[mat;offset;i] mat[i+offset;i]}[mat;offset;] each til dimM}[mat;dimM;]  each til dimM;
        // negative diagonal
        negDiagonal:sum {[mat;dimM;offset]  prd {[mat;offset;i] mat[(count[flip mat]-i)+offset;i]}[mat;offset;] each til dimM}[mat;dimM;]  each til dimM;
        // determinant
        det:posDiagonal-negDiagonal
        ];
    ];
    // dim>3, Laplace expansion, always use first row
    $[dimM>3;
        [
        rowPivot:0;
        signFact: {[rowPivot;x]xexp[neg[1.0];(rowPivot+1+x+1)] }[rowPivot;] each til dimM;
        det:sum mat[rowPivot]*signFact*(.z.s each .quantQ.mat.minorMatrix[mat;rowPivot;] each til dimM);
        ];
    ];
    // output determinant
    :det;
 };

.quantQ.mat.DenmanBeavers:{[bucket;mat]
    // mat -- matrix
    // bucket -- parameters
    bucket:(enlist[`maxIter]!enlist[1000]), bucket;
    unit:.quantQ.mat.diagMatrix[count mat];
    // iteration
    res:({[x]     
        x1:0.5*(first[x]+inv last[x]);
        x2:0.5*(last[x]+inv first[x]);
    :(x1;x2);}/)[bucket[`maxIter];(mat;unit)];
    :first[res];
 };

.quantQ.mat.pow:{[mat;n]
    // mat -- matrix
    // pow -- power
    :last ({[x] (first[x];first[x] mmu last[x])  }/)[n;(mat;.quantQ.mat.diagMatrix[count mat])];
 };

.quantQ.mat.powerSeriesSq:{[bucket;mat]
    // mat -- matrix
    // bucket -- parameters
    bucket:(enlist[`maxIter]!enlist[100]), bucket;
    :sum {[mat;n] 
        // mat -- matrix
        // n -- coeff
        :xexp[-1.0;n]*.quantQ.stats.binomialGen[0.5;n]*.quantQ.mat.pow[.quantQ.mat.diagMatrix[count mat]-mat;n];
    }[mat;] each til bucket[`maxIter];
 }; 
 
.quantQ.mat.covEmpTab:{[tab]
    // tab -- table with features
    :{(1.0%count[x])*flip[x] mmu x}  (flip value flip tab);
 }; 

