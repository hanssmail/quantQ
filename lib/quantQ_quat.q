.quantQ.quat.visual:{[x]
    // x -- quaternion
    :raze (enlist[""],("+";"")(-1=signum 1_value x)),'string[value x],'enlist[""],string 1_key x
 };

.quantQ.quat.mult:{[qt1;qt2]
    // qt1,qt2 -- pair of quaternions
    :`r`i`j`k!(
        {x[`r]-x[`i]-x[`j]-x[`k]}qt1*qt2;
        (qt1[`r]*qt2[`i])+(qt1[`i]*qt2[`r])+(qt1[`j]*qt2[`k])-qt1[`k]*qt[`j];
        (qt1[`r]*qt2[`j])-(qt1[`i]*qt2[`k])+(qt1[`j]*qt2[`r])+qt1[`k]*qt2[`i];
        {x[`r]+x[`i]-x[`j]+x[`k]}qt1*key[qt2]!reverse value qt2);
 };

.quantQ.quat.inverse:{[qt]
    // qt -- quaternion
    :@[qt;`i`j`k;neg]%qt wsum qt;
 };

.quantQ.quat.norm:{[x]
    // x -- quaternion
    :sqrt x wsum x;
 };

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

