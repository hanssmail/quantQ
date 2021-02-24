// quantum simulator

// indexing:
// (a|0>+b|1>) tensor (c|0>+d|1>) = ac|00>+ad|01>+bc|10>+bd|11>;
// index 0 -> (c;d), 1 -> (a;b)
// there are no measurement gates, measurement is thoight to be taken at the end

/////////////////////////////////////////////////  
// Helpers and utility functions

// pi
.quantQ.quantum.pi: acos -1;

// Pauli matrices
.quantQ.quantum.PauliX:(((`re`im)!(0.0;0.0);(`re`im)!(1.0;0.0));((`re`im)!(1.0;0.0);(`re`im)!(0.0;0.0)));
.quantQ.quantum.PauliY:(((`re`im)!(0.0;0.0);(`re`im)!(0.0;neg 1.0));((`re`im)!(0.0;1.0);(`re`im)!(0.0;0.0)));
.quantQ.quantum.PauliZ:(((`re`im)!(1.0;0.0);(`re`im)!(0.0;0.0));((`re`im)!(0.0;0.0);(`re`im)!(neg 1.0;0.0)));

// complex identity matrix
.quantQ.quantum.one:(((`re`im)!(1.0;0.0);(`re`im)!(0.0;0.0));((`re`im)!(0.0;0.0);(`re`im)!(1.0;0.0)));

// complex null matrix
.quantQ.quantum.null:(((`re`im)!(0.0;0.0);(`re`im)!(0.0;0.0));((`re`im)!(0.0;0.0);(`re`im)!(0.0;0.0)));

// U2 gate
.quantQ.quantum.U2:{[params]
    // params -- dictionary with parameters
    params:((`phi`lambda)!(0.5*.quantQ.complex.pi;0.5*.quantQ.complex.pi)),params;
    a12:.quantQ.complex.polar2Canonical[(`radius`angle)!(neg 1.0;params[`lambda])];
    a21:.quantQ.complex.polar2Canonical[(`radius`angle)!(1.0;params[`phi])];
    a22:.quantQ.complex.polar2Canonical[(`radius`angle)!(1.0;params[`lambda]+params[`phi])];
    // output
    :(((`re`im)!(1.0;0.0);a12);(a21;a22));
 };

// U1 gate
.quantQ.quantum.U1:{[params]
    // params -- dictionary with parameters
    params:(enlist[`lambda]!enlist[0.5*.quantQ.complex.pi]),params;
    a22:.quantQ.complex.polar2Canonical[(`radius`angle)!(1.0;params[`lambda])];
    // output
    :(((`re`im)!(1.0;0.0);(`re`im)!(0.0;0.0));((`re`im)!(1.0;0.0);a22));
 };

// Hadamard gate
.quantQ.quantum.H:(((`re`im)!(1.0%sqrt[2];0.0);(`re`im)!(1.0%sqrt[2];0.0));((`re`im)!(1.0%sqrt[2];0.0);(`re`im)!(neg 1.0%sqrt[2];0.0)));

// S gate
.quantQ.quantum.S:(((`re`im)!(1.0;0.0);(`re`im)!(0.0;0.0));((`re`im)!(0.0;0.0);(`re`im)!(0.0;1.0)));

// Sconj
.quantQ.quantum.Sconj:(((`re`im)!(1.0;0.0);(`re`im)!(0.0;0.0));((`re`im)!(0.0;0.0);(`re`im)!(0.0;neg 1.0)));

// T gate
.quantQ.quantum.T:(((`re`im)!(1.0;0.0);(`re`im)!(0.0;0.0));((`re`im)!(0.0;0.0);.quantQ.complex.polar2Canonical[(`radius`angle)!(1.0;0.25*.quantQ.quantum.pi)]));

// Tconj
.quantQ.quantum.Tconj:(((`re`im)!(1.0;0.0);(`re`im)!(0.0;0.0));((`re`im)!(0.0;0.0);.quantQ.complex.polar2Canonical[(`radius`angle)!(1.0;neg 0.25*.quantQ.quantum.pi)]));

// Rx
.quantQ.quantum.Rx:{[params]
    // params -- dictionary with parameters
    params:(enlist[`theta]!enlist[0.5*.quantQ.complex.pi]),params;
    // output
    :(((`re`im)!(cos[0.5*params[`theta]];0.0);(`re`im)!(0.0;neg sin[0.5*params[`theta]]));((`re`im)!(0.0;neg sin[0.5*params[`theta]]);(`re`im)!(cos[0.5*params[`theta]];0.0)));
 };

// Ry
.quantQ.quantum.Ry:{[params]
    // params -- dictionary with parameters
    params:(enlist[`theta]!enlist[0.5*.quantQ.complex.pi]),params;
    // output
    :(((`re`im)!(cos[0.5*params[`theta]];0.0);(`re`im)!(neg sin[0.5*params[`theta]];0.0));((`re`im)!(sin[0.5*params[`theta]];0.0);(`re`im)!(cos[0.5*params[`theta]];0.0)));
 };

// Rz
.quantQ.quantum.Rz:{[params]
    // params -- dictionary with parameters
    params:(enlist[`theta]!enlist[0.5*.quantQ.complex.pi]),params;
    a11:.quantQ.complex.polar2Canonical[(`radius`angle)!(1.0;neg 0.5*params[`phi])];
    a22:.quantQ.complex.polar2Canonical[(`radius`angle)!(1.0;0.5*params[`phi])];
    // output
    :((a11;(`re`im)!(0.0;0.0));((`re`im)!(0.0;0.0);a22));
 };

// dictionary with single quibit operations -- all supported elements here
.quantQ.quantum.singleDict:(`PauliX`PauliY`PauliZ`One`Null`U2`U1`H`S`Sconj`T`Tconj`Rx`Ry`Rz)!(
    (.quantQ.quantum.PauliX;0b);
    (.quantQ.quantum.PauliY;0b);
    (.quantQ.quantum.PauliZ;0b);
    (.quantQ.quantum.one;0b);
    (.quantQ.quantum.null;0b);
    (.quantQ.quantum.U2;1b);
    (.quantQ.quantum.U1;1b);
    (.quantQ.quantum.H;0b);
    (.quantQ.quantum.S;0b);
    (.quantQ.quantum.Sconj;0b);
    (.quantQ.quantum.T;0b);
    (.quantQ.quantum.Tconj;0b);
    (.quantQ.quantum.Rx;1b);
    (.quantQ.quantum.Ry;1b);
    (.quantQ.quantum.Rz;1b)
    );

// multiplication of the complex matrix and the complex vector
.quantQ.quantum.multMatVec:{[mat;vec]
    // mat -- complex matrix
    // vec -- complex vector
    :{[x;y] sum .quantQ.complex.mult'[x;y]}[vec;] each mat;    
 };
// example: mat:.quantQ.quantum.PauliY; vec:(((`re`im)!(1.0;0.0));((`re`im)!(neg 1.0;0.0)));
// .quantQ.quantum.multMatVec[mat;vec]

// multiplication of a constant and matrix
.quantQ.quantum.multConstMat:{[const;mat]
    // const -- complex scalar
    // mat -- complex matrix    
    :{[const;vec] .quantQ.complex.mult[const;] each vec}[const;] each mat;
 };

// tensor product of two complex matrices
.quantQ.quantum.tensorProd:{[mat1;mat2]
    // mat1, mat2 -- two complex matrices, mat1 is 2x2
    // create tensor
    tensor: ((.quantQ.quantum.multConstMat[mat1[0;0];mat2];.quantQ.quantum.multConstMat[mat1[0;1];mat2]);
        (.quantQ.quantum.multConstMat[mat1[1;0];mat2];.quantQ.quantum.multConstMat[mat1[1;1];mat2]));
    // translate into bigger matrix
    tensorMat:raze {[x] {[x;y]  raze (x;y)}'[first x;last x]}   each tensor;
    // output
    :tensorMat;
 };
// exa: mat1: .quantQ.quantum.one; mat2: .quantQ.quantum.one;
// .quantQ.quantum.tensorProd[mat1;mat2]

// get all bit combinations
.quantQ.quantum.getCompState:{[Nquibits]
    // Nquibits -- number of quibits
    :(cross/) Nquibits#enlist[(0;1)];
 };    
// exa:  .quantQ.quantum.getCompState[3]
  
// quibits2number 
.quantQ.quantum.q2number:{[state] 
    // state -- array of bits
    :wsum[reverse[state];{xexp[2;x]} til count state]; 
 }; 

// map between bits and numbers 
.quantQ.quantum.q2nMap:{[Nquibits]
    // Nquibits -- number of quibits 
    :{[t;x] t!x}[t;] .quantQ.quantum.q2number each t:.quantQ.quantum.getCompState[Nquibits];
    // {x[(0 1 1)]}
 }; 

/////////////////////////////////////////////////  
// Quantum systems

// initialise the Bench with unit matrices -- the Bench is a list of matrices, which act upon the quantum system
.quantQ.quantum.setBench:{[Nquibits]
    // Nquibits -- number of quibits
    un:.quantQ.quantum.one;
    $[Nquibits=1;
        benchOut:un;
        benchOut:({[x] :.quantQ.quantum.tensorProd[.quantQ.quantum.one;x];  }/)[Nquibits-1;un] 
    ];
    // output view
    benchView: flip raze[(`index;{ `$"quibit",string[x]} each til Nquibits)]!raze[(0;enlist each Nquibits#enlist un)];
    // dictionary
    :(`bench`view)!(enlist benchOut;benchView);
 };
// bench1: .quantQ.quantum.setBench[Nquibits:3];
// bench1[`view]

// initialise quantum system out of |0> -> |000>,|001>,|010>,|011>...
.quantQ.quantum.setSystem:{[Nquibits]
    // Nquibits -- number of quibits
    un: ((`re`im)!(1.0;0.0);(`re`im)!(0.0;0.0));
    $[Nquibits=1;
        systemOut:un;
        systemOut:(.quantQ.complex.mult/) each (cross/)Nquibits#enlist (enlist each un)    
    ];
    // output
    :systemOut;
 };
 
// Apply the Bench on the state    
.quantQ.quantum.evolveState:{[state;bench]
    // state -- quantum state
    // bench -- the Bench
    // Evolve state
    stateOut: {x[`s]} ({[bucket]
        sOut:.quantQ.quantum.multMatVec[bucket[`b][bucket[`cnt]];bucket[`s]];
        // increase counter 
        cOut:bucket[`cnt]+1;
        // return updated bucket
        :(`b`s`cnt)!(bucket[`b];sOut;cOut);
       }/)[count bench[`bench];(`b`s`cnt)!(bench[`bench];state;0)];
    // return state
    :stateOut;
 };   
// exa: state2:.quantQ.quantum.evolveState[state;bench]  

// Perform measurement    
.quantQ.quantum.measurement:{[state]    
    // state -- quantum state
    // get Nquibits from the bench
    Nquibits: "j"$xlog[2;] count state;
    // probabilities
    probs: {x*x} each .quantQ.complex.norm each state;
    // tabular output
    :([] state: .quantQ.quantum.getCompState[Nquibits];probs:probs);
 };
   
// create the Bench matrix for single quibit operation 
.quantQ.quantum.addSingleQuibitMat:{[Nquibits;matName;id;params]
    // Nquibits -- number of quibits
    // matName -- name of the single quibit matrix
    // id -- index of the affected quibit
    // params -- dictionary with parameters (for some matrices)
    // unit matrix
    un:.quantQ.quantum.one;
    // get the matrix (default if wrong name is unit matrix)
    mat:un;
    if[matName in key[.quantQ.quantum.singleDict];
        $[1b=ps:last .quantQ.quantum.singleDict[matName];
            mat:first[.quantQ.quantum.singleDict[matName]][params];            
            mat:first[.quantQ.quantum.singleDict[matName]]
        ]
    ];
    // iterate tensor product of unit matrices and the matrix
    $[Nquibits=1;
        benchOut:mat;
        benchOut:{x[`out]}({[x] 
        
            if[x[`cnt]=0;
                if[x[`id]=0;newMat:.quantQ.quantum.tensorProd[x[`un];x[`mat]]];
                if[x[`id]=1;newMat:.quantQ.quantum.tensorProd[x[`mat];x[`un]]];
                if[x[`id]>1;newMat:.quantQ.quantum.tensorProd[x[`un];x[`un]]]
            ];
            if[x[`cnt]>0;
                $[x[`id]=x[`cnt]+1;   
                    newMat:.quantQ.quantum.tensorProd[x[`mat];x[`out]];
                    newMat:.quantQ.quantum.tensorProd[x[`un];x[`out]]    
                ]           
            ];
            // return the full bucket    
            :(`out`un`mat`id`cnt)!(newMat;x[`un];x[`mat];x[`id];x[`cnt]+1);
        }/)[Nquibits-1;(`out`un`mat`id`cnt)!(un;un;mat;id;0)] 
    ];
    // output enlisted matrix
    :enlist benchOut;
 };
// exa:  Nquibits:3;matName:`PauliX;id:0;params:()!()
// .quantQ.quantum.addSingleQuibitMat[3;`PauliX;0;()!()]
  
// the Bench is a list of matrices, which act on the quantum system
// append single quibit matrix into the existing bench
.quantQ.quantum.appendSingleQuibitMat:{[bench;matName;id;params]
    // bench -- the bench
    // matName -- name of the single quibit matrix
    // id -- index of the affected quibit
    // params -- dictionary with parameters (for some matrices)
    // get Nquibits from the bench
    Nquibits: "j"$xlog[2;] count first bench[`bench];
    // append new matrix to the bench
    benchOut:bench[`bench],.quantQ.quantum.addSingleQuibitMat[Nquibits;matName;id;params];
    // append view
    benchUpdate: flip raze[(`index;{ `$"quibit",string[x]} each til Nquibits)]!raze[((1+max exec index from bench[`view]);enlist each {[n;i;m] t:n#`$"||";t[i]:m;:t; }[Nquibits;id;matName])];
    // return the updated bench
    :(`bench`view)!(benchOut;(bench[`view],benchUpdate));
    // return the updated bench
    //:bench;
 };
 
// bench:bench1;matName:`PauliX;id:1;params:()!()
// bench2:  .quantQ.quantum.appendSingleQuibit[bench;matName;id;params]
// bench2[`view]
// bench2:  .quantQ.quantum.appendSingleQuibit[bench2;matName;0;params]
// bench2[`view]
// count bench2[`bench]

/////////////////////////////////////////////////  
// Multi-quibit
 
// extract indices of affected elements for CNOT
.quantQ.quantum.cnotInd:{[Nquibits;cid;tid] 
    // Nquibits -- number of quibits
    // cid -- id of control bit
    // tid -- id of target bit 
 
    // map q and number
    q2nMap:.quantQ.quantum.q2nMap[Nquibits];
    // targeted cids
    cidq: {[cid;x] x where 1={[cid;x] x[cid]}[cid;] each reverse each x}[cid;]  key[q2nMap]; 
    // get index of those bits (those rows will be changed)
    cidqN: "j"$q2nMap[cidq];    
    // change the target
    tidq:    {[t;x] x:reverse x;x[t]:"j"$not x[t];:reverse x;}[tid] each cidq;   
    // get index of noted bits (the column ids)
    tidqN: "j"$q2nMap[tidq];
    // return the row/columns map
    :(cidqN;tidqN);
 };    
// exa:  Nquibits:3;cid:0;tid:2
// .quantQ.quantum.cnotInd[Nquibits;cid;tid]    

// create the Bencg matrix for CNOT            
.quantQ.quantum.addCNOTMat:{[Nquibits;cid;tid] 
    // Nquibits -- number of quibits
    // cid -- id of control bit
    // tid -- id of target bit     
    // create unit bench
    unitBench:first .quantQ.quantum.setBench[Nquibits][`bench];
    // get indices
    ctid: .quantQ.quantum.cnotInd[Nquibits;cid;tid]; 
    // calculate the Bench matrix
    benchOut:{x[`b]} ({[ctid;b]  
        bench:b[`b];
        // set diagonal term to 0+i0
        bench[first[ctid][b[`cnt]];first[ctid][b[`cnt]]]:(`re`im)!(0.0;0.0);
        // flip the corresponding term
        bench[first[ctid][b[`cnt]];last[ctid][b[`cnt]]]:(`re`im)!(1.0;0.0);
        // update counter  
        b[`cnt]:b[`cnt]+1;  
        // return updated dictionary
        :(`b`cnt)!(bench;b[`cnt])  
     }[ctid;]/)[count flip ctid;(`b`cnt)!(unitBench;0)];
    // return benchmat with cnot
    :enlist benchOut;
 };
// exa: .quantQ.quantum.addCNOTMat[3;0;1]     
    
// append CNOT to the Bench    
.quantQ.quantum.appendCNOTMat:{[bench;cid;tid]     
    // Nquibits -- number of quibits
    // cid -- id of control bit
    // tid -- id of target bit     
    // get Nquibits from the bench
    Nquibits: "j"$xlog[2;] count first bench[`bench];
    // append new matrix to the bench
    benchOut:bench[`bench],.quantQ.quantum.addCNOTMat[Nquibits;cid;tid];
    // append view
    benchUpdate: flip raze[(`index;{ `$"quibit",string[x]} each til Nquibits)]!raze[((1+max exec index from bench[`view]);enlist each {[n;cid;tid] t:n#`$"||";t[cid]:`cCNOT;t[tid]:`tCNOT;:t; }[Nquibits;cid;tid])];
    // return the updated bench
    :(`bench`view)!(benchOut;(bench[`view],benchUpdate));
  };
 // bench3:  .quantQ.quantum.appendCNOTMat[bench2;0;2]
 // bench3[`view]
 // count bench3[`bench]    
   
// dictionary with double quibit operations -- all supported elements here
.quantQ.quantum.doubleDict:(enlist `CNOT)!(enlist (.quantQ.quantum.appendCNOTMat;0b));
    
/////////////////////////////////////////////////  
// Ideas    
// 1) general wrapper .quantQ.quantum.appendMat
// 2) add more multiquibit mats   

    
    