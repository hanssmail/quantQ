.quantQ.pi: acos -1;

.quantQ.complex.mult:{[x;y]
    :`re`im!(1 -1 wsum x*y;x wsum reverse value y);
 };

.quantQ.complex.div:{[num;den]
    // num -- numerator
    // den -- denominator
    :`re`im!(wsum[num;den]%yy;(1 -1 wsum den*reverse value num)%yy:den wsum den);
 };

.quantQ.complex.conjugate:{[x]
    :@[x;`im;neg];
 };

.quantQ.complex.polar:{[cx]
    // cx -- complex number
    // radial part
    r: sqrt cx wsum cx;
    // angular part
    phi:$[0<cx`re;atan cx[`im]%cx[`re];
        0>cx`re;$[0<=cx`im;
            atan[cx[`im]%cx[`re]] + .quantQ.pi;
            atan[cx[`im]%cx[`re]] - .quantQ.pi];
        0=cx`re;$[0<cx`im;  .5 * .quantQ.pi;
            0>cx`im; -.5 * .quantQ.pi;
            0n];
        0n];
    :`radius`angle!(r;phi);
 };

.quantQ.complex.polar2Canonical:{[polar]
    // polar -- complex number in polar coordinates
    :`re`im!polar[`radius]*(cos;sin)@\:polar`angle;
 };

.quantQ.complex.realExp:{[phi]
    // phi -- angular component
    :`re`im!(cos;sin)@\:phi;
 };

.quantQ.complex.exp:{[cx]
    // cx -- complex number
    :exp[cx`re]*.quantQ.complex.realExp cx`im;
 };

.quantQ.complex.log:{[cx;k]
    // cx -- complex number
    // k -- branching cut
    polar:.quantQ.complex.polar cx;
    :`re`im!(log polar`radius;polar[`angle] + 2 * .quantQ.pi * k);
 };

.quantQ.complex.real2Cx:{[real]
    // real -- real number
    :`re`im!(`float$real;$[0>type real;0f;count[real]#0f]);
 };

.quantQ.complex.dft1k:{[cxs;n;tiln;k]
    // cxs -- series of complex numbers
    // n -- length of the series
    // tiln -- index of the series
    // k -- coefficient of the DFT coefficient
    :sum .quantQ.complex.mult'[cxs;flip .quantQ.complex.realExp neg (2 * .quantQ.pi * k * tiln)%n];
 };

.quantQ.complex.dft:{[cxs]
    // cxs -- series of complex numbers to be transformed
    :.quantQ.complex.dft1k[cxs;N;til N] peach til N:count cxs;
 };

.quantQ.complex.dftTs:{[dfts]
    // dfts -- DFT series
    :update i from dfts,'.quantQ.complex.polar each dfts;
 };

.quantQ.complex.dftTab:{[tab]
    // tab -- sorted time-series table with 2 columns (t and xt)
    t:tab`t;
    xt:update im:0f from `re xcol delete t from tab;
    // sampling frequency
    samplingFrequency:(count[t]-1)%(last[t]-first[t]);
    // dft
    dftts: .quantQ.complex.dftTs .quantQ.complex.dft xt;
    // discrete fourier transform decomposition of the time series
    :(([] f:samplingFrequency*(til count t)%neg[1]+count[t]),'dftts);
 };

.quantQ.complex.fractal:{[cx]
    // cx -- complex offset
    // iterate z_n+1= (z_n^2)+ c and exit if the norm converges or hits max number of iterations
    :{[c;list]
       newp: .quantQ.complex.mult[p;p:first list] + c;
       if[(50<i:last list) or 2<.quantQ.complex.norm newp;:list];
       :(newp;1+i)
    }[cx]/[(`re`im!0 0f;0)];
 };

.quantQ.complex.mandelbrot_:{[p;r;i]
    // p -- starting point for iteration
    // r -- real part
    // i -- imaginary part
    :`row`col`res!(r;i;50<last .quantQ.complex.fractal[`re`im!(r;i)])
 };

.quantQ.complex.mandelbrot:{[h;w]
    // h -- screen height
    // w -- screen width
    // iterate over all rows and columns
    // set the constant c by scaling the row and column and centre to 0
    c:{[h;w;row;col] ((row - .5 * w ) * 4.0 % w; (col - .5 * h) * 4.0 % w)}[h;w]. flip t cross t:til w;
    // iterate over all rows and columns
    :.quantQ.complex.mandelbrot_\[0 0;c 0;c 1]
 };


