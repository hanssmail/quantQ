// building a prototype of the RSA algorithm

// check for prime number
.quantQ.rsa.isPrime:{[prime]
    // prime -- candidate prime number; prime:31
    :prime=1+count {[prime;x] not floor[z]=z:prime%1+x}[prime;]{x+1}\1;   
 };
 
// example .quantQ.rsa.isPrime[31]
 
// find the first prime bigger than x
.quantQ.rsa.primeUp:{[thr]
    // thr -- prime number; thr:14
    :thr+neg[1]+last {[thr;x] not .quantQ.rsa.isPrime (thr-1)+x}[thr;]{x+1}\1; 
  };
// example .quantQ.rsa.primeUp[14]

// Euclidean algorithm -- greatest common divisor
.quantQ.rsa.gcd:{[x1;x2]
    // x1, x2 -- two integers to find their common divisior x1:14;x2:14
    :{x[`a]}({[bucket]
        // sort if needed
        if[bucket[`a]<bucket[`b]; tmp:bucket[`a];bucket[`a]:bucket[`b];bucket[`b]:tmp];
        tmp: bucket[`b];
        res:mod[bucket[`a];bucket[`b]];
        bucket[`a]:tmp;
        bucket[`b]:res;
        :bucket;
        }/)[{x[`b]>0};(`a`b)!(x1;x2)];
 };
// example  .quantQ.rsa.gcd[14;23]
 
// test for co-prime
.quantQ.rsa.isCoprime:{[x1;x2]
    // x1, x2 -- two integers to check for being coprime x1:14;x2:14
    :1=.quantQ.rsa.gcd[x1;x2];
 };
// example  .quantQ.rsa.isCoprime[14;28]

// Carmichael function
.quantQ.rsa.CarmichaelFunc:{[n]
    // n -- positive integer; n:14
    allCoprimes:z where .quantQ.rsa.isCoprime[n;] each z:1+til n;
    :last {[cp;n;x] not prd {[cp;n;x]  1=mod["j"$xexp[cp;x];n]}[;n;x] each cp   }[allCoprimes;n;]{x+1}\1;   
 };
// example .quantQ.rsa.CarmichaelFunc[14]

// least common multiplier
.quantQ.rsa.lcm:{[x1;x2]
    // x1,x2 -- two integers
    :$[0=abs[x1]+abs[x2];0;abs[x1*x2]%.quantQ.rsa.gcd[x1;x2]];
 };
// example .quantQ.rsa.lcm[6;14]

// modular inverse
.quantQ.rsa.inverseMod:{[e;n]
    // e -- number to invert; e:17        
    // n -- modulo; n:21798
    sol:({[bucket]
        quot:floor bucket[`r]%bucket[`rn];
        bucketTMP:bucket;
        bucket[`t]:bucketTMP[`tn];
        bucket[`tn]:bucketTMP[`t]-quot*bucketTMP[`tn];
        bucket[`r]:bucketTMP[`rn];
        bucket[`rn]:bucketTMP[`r]-quot*bucketTMP[`rn];
        :bucket; 
    }/)[{x[`rn]<>0};(`t`tn`r`rn)!(0;1;n;e)];                
    sol[`status]:1;
    if[sol[`r]>1;sol[`status]:0];
    if[sol[`t]<0;sol[`t]:sol[`t]+n];   
    // return solution, t is the inverse
    :sol;                       
 };
// example: .quantQ.rsa.inverseMod[17;21798]

// modular multiplicative exponent
.quantQ.rsa.modExp:{[b;e;n]
    // parameters: b^e mod n; b:4;e:13;n:497;
    :{x[`sol]}({[bucket]
        bucket[`eBar]:bucket[`eBar]+1; 
        bucket[`sol]:mod[bucket[`sol]*bucket[`b];bucket[`md]];
        :bucket;
     }/)[{x[`eBar]<x[`e]};(`e`eBar`sol`b`md)!(e;0;1;b;n)];
  };
// example: .quantQ.rsa.modExp[4;13;497]    

// generate RSA keys
.quantQ.rsa.genKeys:{[bucket]
    // bucket -- parameters; bucket:()!()    
    // set default 
    bucket:((`pMin`qMin`eStarting)!(123;345;17)),bucket;
    // calculate p and q
    p:.quantQ.rsa.primeUp[bucket[`pMin]];
    q:.quantQ.rsa.primeUp[bucket[`qMin]];
    // set the RSA object
    rsa:(`p`q`status)!(p;q;1);
    // add modulus
    rsa[`modulus]:rsa[`p]*rsa[`q];
    // calculate Carmichael function of modulus  
    rsa[`lam]:"j"$.quantQ.rsa.lcm[rsa[`p]-1;rsa[`q]-1];
    // find e
    rsa[`e]:{[es;x] es+x-1}[bucket[`eStarting];] last {[ln;es;x] not .quantQ.rsa.isCoprime[ln;es+x-1]}[rsa[`lam];bucket[`eStarting];]{x+1}\1;
    // calculate d
    invM:.quantQ.rsa.inverseMod[rsa[`e];rsa[`lam]];
    rsa[`d]:0;
    rsa[`status]:0;
    if[invM[`status]=1;rsa[`d]:invM[`t];rsa[`status]:1];
    // add the public part
    rsaOut:(`p`q`status`modulus`lam`e`d`publicKey)!(rsa[`p];rsa[`q];rsa[`status];rsa[`modulus];rsa[`lam];rsa[`e];rsa[`d];(`modulus`e)!(rsa[`modulus];rsa[`e]));
    :rsaOut;   
 };
// example .quantQ.rsa.genKeys[()!()]

// encryption
.quantQ.rsa.encrypt:{[rsaPublicKey;m]
    // rsaPublicKey -- public key
    // m -- message to encode
    cipher:0;
    // algorithm can cipher only m<modulus
    if[m<rsaPublicKey[`modulus];
        cipher:.quantQ.rsa.modExp[m;rsaPublicKey[`e];rsaPublicKey[`modulus]];
    ];
    :(`rsaPublicKey`cipher)!(rsaPublicKey;cipher);
 };
// example rsaX:.quantQ.rsa.genKeys[()!()]; .quantQ.rsa.encrypt[rsaX[`publicKey];123]

// decryption
.quantQ.rsa.decrypt:{[message;rsa]
    // message -- the encrypted message including the public key
    // rsa -- full key
    // return decrypted message    
    :.quantQ.rsa.modExp[message[`cipher];rsa[`d];rsa[`modulus]];
 };
// example .quantQ.rsa.decrypt[.quantQ.rsa.encrypt[rsaX[`publicKey];123];.quantQ.rsa.genKeys[()!()]]
