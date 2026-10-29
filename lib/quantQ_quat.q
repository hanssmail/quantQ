.quantQ.quat.visual:{[x]
    // x -- quaternion
    :raze (enlist[""],("+";"")(-1=signum 1_value x)),'string[value x],'enlist[""],string 1_key x
 };

.quantQ.quat.mult:{[qt1;qt2]
    // qt1,qt2 -- pair of quaternions
    // Hamilton product; every term is parenthesised because q evaluates right to left
    a1:qt1`r; b1:qt1`i; c1:qt1`j; d1:qt1`k;
    a2:qt2`r; b2:qt2`i; c2:qt2`j; d2:qt2`k;
    :`r`i`j`k!(
        (a1*a2)-((b1*b2)+(c1*c2)+(d1*d2));
        ((a1*b2)+(b1*a2)+(c1*d2))-(d1*c2);
        ((a1*c2)+(c1*a2)+(d1*b2))-(b1*d2);
        ((a1*d2)+(b1*c2)+(d1*a2))-(c1*b2));
 };

.quantQ.quat.inverse:{[qt]
    // qt -- quaternion
    :@[qt;`i`j`k;neg]%qt wsum qt;
 };

.quantQ.quat.norm:.quantQ.complex.norm;

.quantQ.quat.conjugate:{[x]
    // x -- quaternion
    :@[x;`i`j`k;neg];
 };

.quantQ.quat.polarDecomp:{[qt]
    // qt -- quaternion
    norm: .quantQ.quat.norm qt;
    // return a norm and unit quaternion
    :(norm; qt%norm);
 };

.quantQ.quat.crossProduct:{[qt1;qt2]
    // qt1, qt2 -- quaternion
    // the real part is ignored and in the output set to zero
    :`r`i`j`k!(0.0;(qt1[`j]*qt2[`k])-(qt1[`k]*qt2[`j]);(qt1[`k]*qt2[`i])-(qt1[`i]*qt2[`k]);
        (qt1[`i]*qt2[`j])-(qt1[`j]*qt2[`i]));
 };

.quantQ.quat.exp:{[qt]
    // qt -- quaternion
    realPart: qt`r;
    imPart: `r _qt;
    newRealPart: enlist[`r]!enlist cos sqrt imPart wsum imPart;
    newImPart: $[0=normImPart:sqrt sum imPart*imPart;imPart;(imPart%normImPart)*
        sin[sqrt sum imPart*imPart]];
    // return exp of quaternion
    :exp[realPart]*newRealPart,newImPart;
 };

.quantQ.quat.log:{[qt]
    // qt -- quaternion
    realPart: qt`r;
    imPart: `r _qt;
    // exp of quaternion
    newRealPart: enlist[`r]!enlist log .quantQ.quat.norm qt ;
    newImPart: $[0=normImPart:sqrt sum imPart*imPart;0.0*imPart;
        (imPart%normImPart)]*acos[realPart % .quantQ.quat.norm qt];
    :newRealPart,newImPart;
 };

