
.quantQ.tse.simAR:{[phi;eps]
    // equivalent to .quantQ.ts.simAR but using scan
    // phi -- list of parameters (phi1;...,phip)
    // eps -- list of n random normal numbers, use .quantQ.simul.getNormalVariate n
    r:{[n;phi;x;eps] (eps+phi $ x),n#x}[cp-1;reverse phi]\[(cp:count phi)#0f;eps];
    :r[;0]
 };
