.quantQ.ts.simulate:{[process;params;eps]
    // process -- type of process to simulate. Accepted values: `AR,`MA or `ARMA
    // params -- list of floating numbers containing the parameters of the process
    // eps -- list of random normal numbers
    $[process=`AR;
        $[9h=type params;
            .quantQ.ts.simAR[params;eps];
            '"wrong type for params"      // note the single tick (') which is used to return an exception
	];
        $[process=`MA;
            $[9h=type params;
                .quantQ.ts.simMA[count params;params;eps];
                '"wrong type of params"
            ];
            $[process=`ARMA;
                $[(2=count params) & (9h = type params[0]) & (9h = type params[1]);
                    .quantQ.ts.simARMA[params[0];params[1];eps];
                    '"wrong type of params"
                    ];
                    '"process is `AR, `MA or `ARMA"
            ]
        ]
    ]
 };

.quantQ.ts.simAR:{[phi;eps]
    // phi -- list of parameters (phi1;...,phip)
    // eps -- list of random normal numbers
    p:count phi;
    x:(count eps)#0f;
    t:p;
    while[t<(count eps);
        // here all the error terms are used, can be replace with eps[t]
        x[t]:(x[(t-p)+til p] mmu phi)+eps[t-p];
        t:t+1];
    :x;
 };

.quantQ.ts.simMA:{[q;alpha;eps]
    // q -- order of the MA(q) process to simulate. q is an integer.
    // alpha -- list of parameters of the MA process (alpha1,...,alphaq)
    // eps -- list of random normal numbers
    :sum flip $[q>1;
        (,')[(q _eps,q#0);((q-1) _alpha[q-1]*eps,(q-1)#0)],'.quantQ.ts.simMA[q-1;alpha;eps];
        (,')[(q _eps,0);(q-1) _alpha[q-1]*eps]];
 };

.quantQ.ts.simARMA:{[phi;alpha;eps]
    // phi -- list of parameters of the AR(p) process.
    // alpha -- list of parameters of the MA(q) process
    // eps -- list of random numbers
    p:count phi;
    q:count alpha;
    x:(count eps)#0f;
    t:p;
    while[t<(count eps);
        x[t]:(x[(t-p)+til p] mmu phi)+(eps[(t-q)+til q] mmu alpha);
        t:t+1];
    :x;
 }

.quantQ.ts.estimate:{[data;process;lags]
    // data -- historical data (list of floating numbers) used to estimate the process. 
    // process -- type of process to simulate. Accepted values: `AR,`MA or `ARMA
    // lags -- order of the process. integer for AR and MA process. List of two integers for ARMA
    $[process=`AR;
        $[-7h = type lags;
            .quantQ.ts.DurbinLevinson[data;lags];
            $[7h = type lags;
                .quantQ.ts.DurbinLevinson[data;first lags];
                '"wrong type of the parameter lags"]
         ];
        $[process=`MA;
            $[-7h = type lags;
                .quantQ.ts.innovations[data;lags];
                $[7h = type lags;
                    quantQ.ts.innovations[data;first lags];
                    '"wrong type of the parameter lags"]
            ];
            $[process=`ARMA;
                $[7h = type lags;
                    .quantQ.ts.HannanRissanen[data;lags[0];lags[1]];
                    '"The estimation of ARMA(p,q) processes requires a list of two integers for the lag parameter"
                ];
                '"unknown process type. Accepted values: `AR; `MA and `ARMA"
            ]
        ]
    ]
 };

.quantQ.ts.DurbinLevinson:{[data;p]
    // data -- historical data (list of floating numbers) used to estimate the AR(p) process.
    // p -- order of the AR(p) model.
    p:p+1;
    phi:(p;p)#0f;
    v:p#0f;
    phi[1;1]:.quantQ.ts.gamma[1;data]%.quantQ.ts.gamma[0;data];
    v[0]:.quantQ.ts.gamma[0;data];
    n:1;
    s:0f;
    while[n<p;
        j:1;
        while[j<n;
            s:s+phi[n-1;j]*.quantQ.ts.gamma[n-j;data];
            j:j+1];
        phi[n;n]:(.quantQ.ts.gamma[n;data]-s)%v[n-1];
        s:0f;
        v[n]:v[n-1]*(1f-(phi[n;n] xexp 2));
        j:1;
        while[j<n;
            phi[n;j]:phi[n-1;j]-phi[n;n]*phi[n-1;n-j];
            j:j+1];
        n:n+1];
    :1_phi[p-1;];
 };

.quantQ.ts.gamma:{[h;data] 
    // h -- order of autocorrelation
    // data -- data
    :(neg[h]_data) cov (h _data);
 };

.quantQ.ts.r:{[data;h]
    // h -- order of autocorrelation
    // data -- data
    :.quantQ.ts.gamma[h;data]%.quantQ.ts.gamma[0;data];
 };

.quantQ.ts.toeplitz:{[data;p]
    // data -- historical data (list of floating numbers) 
    // p -- order of the AR(p) process
    t:abs (til p)-/:(til p);
    :(p;p)#{[x;data].quantQ.ts.r[data;x]}[;data] each raze t;
 };

.quantQ.ts.yuleWalker:{[data;p]
    // data -- historical data (list of floating numbers)
    // p -- order of the AR(p) model
    :inv[.quantQ.ts.toeplitz[data;p]] mmu (.quantQ.ts.r[data;] each 1+til p);
 };

.quantQ.ts.estARP:{[data;p]
    // data -- historical data (list of floating numbers) used to estimate the AR(p) model
    // p -- order of the AR(p) process
    n:count data;
    dim:(n-p;p+1);
    x:dim#0f;
    i:0;j:1;
    while[i<(n-p);
        x[i;0]:1f;
        while[j<p+1;
            $[0<=i-j;
                x[i;j]:data[i-j];
            x[i;j]:0f];
         j:j+1];
     j:1;
     i:i+1];
     y:neg[p]_ data;
     :inv[(flip x) mmu x] mmu (flip x) mmu y;
 };

.quantQ.ts.innovations:{[data;q]
    // data -- historical data (list of floating numbers) used to estimate the MA(q) process
    // q -- order of the MA(q) process
    p:10*q; // until convergence
    v:(p+1)#0f;
    alpha:(p+1;p+1)#0f;
    v[0]:.quantQ.ts.gamma[0;data];
    n:1;
    k:0;
    while[n<=p;
        k:0;
        while[k<n;
            j:0; s:0f;
            while[j<k;
                s:s+(alpha[k;k-j]*alpha[n;n-j]*v[j]);
                j:j+1;
            ];
            alpha[n;n-k]:(.quantQ.ts.gamma[n-k;data] - s)%v[k];
            k:k+1;
            s:0f; j:0;
            while[j<n;
                s:s+(alpha[n;n-j] xexp 2)*v[j];
                j:j+1;
            ];
            v[n]:.quantQ.ts.gamma[0;data]-s;
        ];
        n:n+1;
    ];
    :alpha[p-1;1+ til q];
 };

.quantQ.ts.HannanRissanen:{[data;p;q]
    // data -- historical data (list of floating numbers) used to estimate the ARMA(p,q) process
    // p -- order of the AR(p) process
    // q -- order of the MA(q) process
    // Step 1 - Estimate AR(max(p+q)+1)
    lags:1+(p|q);
    arEst:.quantQ.ts.DurbinLevinson[data;lags];
    res:.quantQ.ts.residualsAR[arEst;data];
    res:(lags#0f),res;
    // Step 2 - regress data[t] on (data[t-1],...,data[t-p],res[t-1],...,res[t-q])
    // for t = lag+q, ..., n by OLS
    t:(lags+q) _ til count data;
    x:(reverse each .quantQ.ts.matrix[t;p;data]),'(reverse each .quantQ.ts.matrix[t;q;res]);
    :(inv[(flip x) mmu x] mmu flip x) mmu data[t];
 };

.quantQ.ts.matrix:{[t;q;d]
    // t -- index
    // q -- order of the MA(q) process
    // d -- series
    :$[q>1; d[t-q],'.quantQ.ts.matrix[t;q-1;d];d[t-1]];
 };

.quantQ.ts.forecast:{[data;period;process;lags]
    // data -- historical data (list of floating numbers) used to forecast the process
    // period -- forecast period-step ahead the process where period is a positive integer.
    // process -- type of process to simulate. Accepted values: `AR,`MA or `ARMA
    // lags -- order of the process. integer for AR and MA process. List of two integers for ARMA.
    $[process=`AR;
        $[(-7h=type period) & (period > 0);
            $[-7h = type lags;
            	.quantQ.ts.forecastAR[data;period;lags];
                $[7h = type lags;
                   .quantQ.ts.forecastAR[data;period;first lags];
                   '"wrong type of the parameter lags for AR processes"]
            ];
            '"period needs to be a positive integer"
        ];
        $[process=`MA;
            $[(-7h=type period) & (period > 0);
                $[-7h = type lags;
                    .quantQ.ts.forecastMA[data;period;lags];
                    $[7h = type lags;
                        .quantQ.ts.forecastMA[data;period;first lags];
                        '"wrong type of the parameter lags for MA processes"]
                ];
                '"period needs to be a positive integer"
            ];
            $[process=`ARMA;
                $[(-7h=type period) & (period > 0);
                    $[7h = type lags;
                        .quantQ.ts.forecastARMA[data;period;lags[0];lags[1]];
                        '"lags needs to be a list of two positive integers for ARMA(p,q) process"
                    ];
                    '"period needs to be a positive integer"
                ];
                '"unknown process type. Accepted values `AR; `MA and `ARMA"
            ]
        ]
    ]
 };

.quantQ.ts.forecastAR:{[data;period;lags]
    // data --  historical data (list of floating) used to forecast the AR(lags) model
    // period -- number of steps (positive integer) forecasted  
    // lags -- positive integer containing the order of the autoregressive process
    forecast:((neg lags)#data) mmu .quantQ.ts.DurbinLevinson[data;lags];
    data:data,forecast;
    $[period=1;
        :forecast;
        :forecast,.quantQ.ts.forecastAR[data;period-1;lags]
    ]
 };

.quantQ.ts.forecastMA:{[data;q]
    // data --  historical data (list of floating) used to forecast the AR(lags) model
    // q -- degree of MA
    est:(q+2)#0f;
    v:(q+1)#0f;
    alpha:(q+1;q+1)#0f;
    v[0]:.quantQ.ts.gamma[0;data];
    n:1;
    k:0;
    while[n<=q;
        k:0;
        while[k<n;
            j:0; s:0f;
            while[j<k;
                s:s+(alpha[k;k-j]*alpha[n;n-j]*v[j]);
                j:j+1;
            ];
        alpha[n;n-k]:(.quantQ.ts.gamma[n-k;data] - s)%v[k];
        k:k+1;
        s:0f; j:0;
            while[j<n;
                s:s+(alpha[n;n-j] xexp 2)*v[j];
                j:j+1;
            ];
            v[n]:.quantQ.ts.gamma[0;data]-s;
        ];
        n:n+1;
    ];
    // one-step predictor of MA(q)
    n:1;
    while[n <= q; 
        j:1;s:0f;
        while[j <= n;
            s:s+(alpha[n;j]*(data[j]-est[j]));
            j:j+1];
        est[n+1]:s;
        n:n+1;
        s:0f
    ];
    :1_est;
 };

.quantQ.ts.vec:{[a]
    // a -- matrix	
    i:(count raze a);
    j:1;
    v:(flip a)[0;];
    while[j<i;
        v:v,/(flip a)[j;];
        j:j+1];
    :v;
 };

.quantQ.ts.zt:{[y;t;p]
    // y -- historical data (list of floating numbers) 
    // t -- time point (positive integer)
    // p -- lags (positive integer) corresponding to the order of the VAR(p) process
    :1f,.quantQ.ts.vec[y[;(t-1)-til p]];
 };

.quantQ.ts.Z:{[y;p]
    // y -- historical data (list of floating numbers) 
    // p -- lags (positive integer) corresponding to the order of the VAR(p) process
    :flip {[y;i;p]
        .quantQ.ts.zt[y;i;p]}[y;;p] each p + til (count y[0])-p;
 };

.quantQ.ts.Y:{[y;p]
    // y -- historical data (list of floating numbers) 
    // p -- lags (positive integer) corresponding to the order of the VAR(p) process
    :y[;p+til (count y[0;])-p];
 };

.quantQ.ts.varpest:{[y;p]
    // y -- historical data (list of floating numbers) used to estimate the VAR(p) model.
    // p -- order of the VAR(p) model
    :.quantQ.ts.Y[y;p] mmu (flip .quantQ.ts.Z[y;p]) mmu
        inv[.quantQ.ts.Z[y;p] mmu flip .quantQ.ts.Z[y;p]];
 };

.quantQ.ts.covarianceResidual:{[y;p]
    // y -- historical data (list of floating numbers)
    // p -- order of the VAR(p) process
    T:(count y[0])-p;
    K:count y;
    coef:1f%((T-K*p)-1f);
    :coef*(.quantQ.ts.Y[y;p] mmu flip .quantQ.ts.Y[y;p]) -
        (.quantQ.ts.Y[y;p] mmu flip .quantQ.ts.Z[y;p]) mmu
        inv[.quantQ.ts.Z[y;p] mmu flip .quantQ.ts.Z[y;p]] mmu
        (.quantQ.ts.Z[y;p] mmu flip .quantQ.ts.Y[y;p]);
 };

.quantQ.ts.eye:{[k]
    // k -- rank of matrix
    :`float$(til k)=/:til k;
 };

.quantQ.ts.beta:{[y;p]
    // y -- historical data (list of floating numbers)
    // p -- order of the VAR(p) process
    :(1f%(count y[0])-p) * (.quantQ.ts.Y[y;p] mmu
    (.quantQ.ts.eye[(count y[0])-p] - ((flip .quantQ.ts.Z[y;p]) mmu
        inv[.quantQ.ts.Z[y;p] mmu flip .quantQ.ts.Z[y;p]] mmu
        .quantQ.ts.Z[y;p]))) mmu flip .quantQ.ts.Y[y;p];
 };

.quantQ.ts.varOrder:{[y;m]
    // y -- historical data (list of floating numbers
    // m -- maximum order of the VAR process
    T:count y[0];
    K:count y;
    i:m;
    k:v1:v2:();
    while[i>0;
        lr:T*(log .quantQ.mat.det[.quantQ.ts.beta[y;i-1]])-(log .quantQ.mat.det[.quantQ.ts.beta[y;i]]);
        k:k,(m-i);
        v1:v1,lr;
        v2:v2,.quantQ.ts.pValueChi2[lr;K*K];
        i:i-1];
    :([p:k]lr:v1;pvalue:v2);
 };

.quantQ.ts.pValueChi2:{[chi;nu]
    // chi -- instance of distribution (floating)
    // nu -- degrees of freedom of the Chi2 distribution
    u:chi%2f;
    v:nu%2f;
    p:(u-v)+1f;
    pi:3.141592654f;
    term1:exp[v-u]%(p*sqrt[2*pi]);
    term2:(u % v) xexp v;
    term3:1f-((v-1f)%((p*p)+2f*u));
    term4:(12f*(v xexp 1.5f))%((12f*v)+1);
    :term1*term2*term3*term4;
 };

.quantQ.ts.tstat:{[y;p]
    // y -- historical data (list of floating numbers)
    // p -- order of the VAR(p) process
    :.quantQ.ts.kron[.quantQ.ts.covarianceResidual[y;p]; inv[.quantQ.ts.Z[y;p] mmu
        flip .quantQ.ts.Z[y;p]]];
 };

.quantQ.ts.tRatios:{[y;p]
    // y -- historical data (list of floating numbers)
    // p -- order of the VAR(p) process
    varp:.quantQ.ts.varpest[y;p];
    diag:(count varp;count varp[0])#sqrt[.quantQ.ts.diag[.quantQ.ts.tstat[y;p]]];
    :varp%diag;
 };

.quantQ.ts.kron:{[a;b]
    // a -- matrix of floating numbers
    // b -- matrix of floating numbers
    // calculating Kroenecker product
    nbcols:(count b[0;])*(count a[0;]);
    nbrows:(count b[;0])*(count a[;0]);
    :(nbrows;nbcols) #(raze over) {[i;a;b](,/)(flip (flip a)[;i]*\:b)}[;a;b] each til count a;
 };

.quantQ.ts.diag:{[m]
    // m -- array of numbers to put on diagonal	
    :{[i;m] m[i;i]}[;m] each til count m;
 };

.quantQ.ts.estimate1:{[y;p]
    // y -- historical data (list of floating numbers)
    // p -- order of the VAR(p) process
    b:.quantQ.ts.varpest[y;p];
    v:b[;0];
    k:count b;
    :v+sum {[i;y;b;k]b[;(1+k*i)+til k] mmu y[;(count y[0])-(1+i)]}[;y;b;k] each til p;
 };

.quantQ.ts.estimateN:{[n;y;p]
    // n -- positive integer indicating a n-step ahead forecast
    // y -- historical data (list of floating numbers)
    // p -- order of the VAR(p) process
    output:.quantQ.ts.estimate1[y;p];
    :$[n>1;
        .quantQ.ts.estimateN[n-1;y,'output;p];
        output];
 };

.quantQ.ts.forecastMSE:{[y;p]
    // y -- historical data (list of floating numbers)
    // p -- order of the VAR(p) process
    T:(count y[0])-p;
    K:count y;
    coef:(T+(K*p)+1)%T;
    :coef*.quantQ.ts.covarianceResidual[y;p];
 };

.quantQ.ts.confInterval:{[y;p]
    // y -- historical data (list of floating numbers)
    // p -- order of the VAR(p) process
    mse:1.96*sqrt[.quantQ.ts.diag[.quantQ.ts.forecastMSE[y;p]]];
    forecast:.quantQ.ts.estimate1[d;2];
    :flip ((forecast-mse);(forecast+mse));
 };

.quantQ.ts.zt:{[y0;t;p;flag]
    // y0 -- vectorization of y0 
    // t -- time (positive integer)
    // p -- number of lags (positive integer)
    // flag -- boolean value deciding whether a leading 1 needs to be set.
    :$[flag;
        1f,.quantQ.ts.vec[y0[;(t-1)-til p]];
        .quantQ.ts.vec[y0[;(t-1)-til p]]];
 };

 .quantQ.ts.Z:{[y;p;flag]
    // y -- historical data of the endogenous variable (list of floating numbers)
    // p -- order of the VAR(p) process
    // flag -- boolean value deciding whether a leading 1 needs to be set for zt.
    :flip {[y;i;p;flag].quantQ.ts.zt[y;i;p;flag]}[y;;p;flag] each p + til (count y[0])-p;
 };

.quantQ.ts.ZX:{[x;y;q;p]
    // x -- historical data of the exogenous variable (list of floating numbers)
    // y -- historical data of the endogenous variable (list of floating numbers)
    // q -- order of exogenous process
    // p -- order of the VAR(p) process
    :.quantQ.ts.Z[y;p;1b],.quantQ.ts.Z[x;q;0b],x[;q+til (count x[0])-q];
 };

.quantQ.ts.Y:{[y;p]
    // y -- historical data of the endogenous variable (list of floating numbers)
    // p -- order of the VAR(p) process
    :y[;p+til (count y[;0])+1-p];
 };

.quantQ.varXest:{[x;y;q;p]
    // x -- historical data of the exogenous variable (list of floating numbers)
    // y -- historical data of the endogenous variable (list of floating numbers)
    // q -- order of the exogenous process
    // p -- order of the VAR(p) process
    z0:.quantQ.ts.ZX[x;y;q;p];
    :y[;p+til (count y[0])-p] mmu flip z0 mmu inv (flip z0) mmu z0;
 };

.quantQ.ts.covarianceResidual:{[x;y;q;p]
    // x -- historical data of the exogenous variable (list of floating numbers)
    // y -- historical data of the endogenous variable (list of floating numbers)
    // q -- order of exogenous process
    // p -- order of the VAR(p) process
    yp0:.quantQ.ts.Y[y;p];
    z0:.quantQ.ts.ZX[x;y;q;p];
    T:count y[0];
    K:count y;
    coef:1%(T-(K*p)-1f);
    :coef*(yp0 mmu flip yp0) - yp0 mmu (flip z0) mmu inv[z0 mmu flip z0] mmu z0 mmu flip yp0;
 };

.quantQ.ts.estimateX:{[x;y;q;p]
    // x -- historical data of the exogenous variable (list of floating numbers)
    // y -- historical data of the endogenous variable (list of floating numbers)
    // q -- order of exogenous process
    // p -- order of the VAR(p) process
    M:.quantQ.ts.varXest[x;y;q;p];
    v:M[;0];
    k:count y;
    m:count x;
    :v+sum[{[i]M[;(1+k*i)+til k] mmu y[;(count y[0])-(1+i)]} each til p]+
        sum[{[i]M[;1+(k*p)+(m*i)+til m] mmu x[;(count x[0])-(1+i)]} each til q]+
        M[;1+(k*p)+(m*q)+til m] mmu x[;0];
 };

.quantQ.ts.estimateNX:{[n;x;y;q;p]
    // n -- positive integer indicating the number of steps of the foreast
    // x -- historical data of the exogenous variable (list of floating numbers)
    // y -- historical data of the endogenous variable (list of floating numbers)
    // q -- order of exogenous process
    // p -- order of the VAR(p) process
    output:.quantQ.ts.estimateX[x;y;q;p];
    :$[n>1;
    {
        newx:.quantQ.ts.estimateN[(n-1);q;x];
        .quantQ.ts.estimateNX[n-1;x,'newx;y,'output;q;p];
    };
    output];
 };
