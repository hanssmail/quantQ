// Rosenbrock's valley (banana) function:
.quantQ.func.rosenbrock: {[a;b;x] t1: a-x[0]; t2: x[1]-x[0]*x[0]; (t1*t1) + b*t2*t2};

.quantQ.func.himmelblau: {[x] t1: (x[0]*x[0]) + x[1] - 11f; t2: x[0] + (x[1]*x[1]) - 7f; (t1*t1) + t2*t2};

.quantQ.amoeba_.amotry: {[ihi;fac;func]
    fac1: (1f-fac)%.quantQ.amoeba_.ndim;
    fac2: fac1-fac;
    ptry: (.quantQ.amoeba_.psum*fac1) - (.quantQ.amoeba_.p[ihi]*fac2);
    ytry: func[ptry];
    $[ytry < .quantQ.amoeba_.y[ihi]; [
        .quantQ.amoeba_.y[ihi]: ytry;
        .quantQ.amoeba_.psum: .quantQ.amoeba_.psum + (ptry - .quantQ.amoeba_.p[ihi]);
        .quantQ.amoeba_.p[ihi]: ptry;
    ];];
    ytry
  };

.quantQ.amoeba_.amoeba: {[func]
    ilo: 0;
    $[.quantQ.amoeba_.y[0]>.quantQ.amoeba_.y[1]; [ihi:0; inhi:1]; [ihi:1; inhi:0]];
    i: 0;
    while[i<.quantQ.amoeba_.mpts;
        $[.quantQ.amoeba_.y[i]<=.quantQ.amoeba_.y[ilo]; ilo:i;];
        $[.quantQ.amoeba_.y[i]>.quantQ.amoeba_.y[ihi]; [inhi:ihi; ihi:i]; $[and[.quantQ.amoeba_.y[i]>.quantQ.amoeba_.y[inhi]; i<>ihi]; inhi:i;]];
        i: i+1;
    ];

    rtol: (2f*abs[.quantQ.amoeba_.y[ihi]-.quantQ.amoeba_.y[ilo]]) % (abs[.quantQ.amoeba_.y[ihi]]+abs[.quantQ.amoeba_.y[ilo]]+.quantQ.amoeba_.TINY);

    $[rtol<.quantQ.amoeba_.ftol; [
        temp: .quantQ.amoeba_.y[0]; .quantQ.amoeba_.y[0]: .quantQ.amoeba_.y[ilo]; .quantQ.amoeba_.y[ilo]: temp;
        i: 0;
        while[i<.quantQ.amoeba_.ndim;
            temp: .quantQ.amoeba_.p[0;i]; .quantQ.amoeba_.p[0;i]: .quantQ.amoeba_.p[ilo;i]; .quantQ.amoeba_.p[ilo;i]: temp;
            .quantQ.amoeba_.pmin[i]: .quantQ.amoeba_.p[0;i];
            i: i+1;
        ];
        `.quantQ.amoeba_.fmin set .quantQ.amoeba_.y[0];
        :.quantQ.amoeba_.fmin;
    ];];

    $[.quantQ.amoeba_.nfunc>=.quantQ.amoeba_.NMAX; '`MAX_EVALUATIONS_EXCEEDED;];
    `.quantQ.amoeba_.nfunc set .quantQ.amoeba_.nfunc+2;
    ytry: .quantQ.amoeba_.amotry[ihi;-1f;func];
    $[ytry<=.quantQ.amoeba_.y[ilo]; [
        ytry: .quantQ.amoeba_.amotry[ihi;2f;func];
    ]; [ $[ytry>=.quantQ.amoeba_.y[inhi]; [
        ysave: .quantQ.amoeba_.y[ihi];
        ytry: .quantQ.amoeba_.amotry[ihi;.5;func];
        $[ytry>=ysave; [
            i: 0;
            while[i<.quantQ.amoeba_.mpts;
                $[i <> ilo; [
                    `.quantQ.amoeba_.psum set .5*(.quantQ.amoeba_.p[i] + .quantQ.amoeba_.p[ilo]);
                    .quantQ.amoeba_.p[i]: .quantQ.amoeba_.psum;
                    .quantQ.amoeba_.y[i]: func[.quantQ.amoeba_.psum];
                ];];
                i: i+1;
            ];
            `.quantQ.amoeba_.nfunc set .quantQ.amoeba_.nfunc + .quantQ.amoeba_.ndim;
            `.quantQ.amoeba_.psum set sum .quantQ.amoeba_.p;
        ];];
    ]; [`.quantQ.amoeba_.nfunc set .quantQ.amoeba_.nfunc-1;] ]; ] ];
  };

.quantQ.amoeba: {[ftol;p;func]
    `.quantQ.amoeba_.ftol set ftol;
    `.quantQ.amoeba_.NMAX set 5000;
    `.quantQ.amoeba_.TINY set 1.0e-10;
    `.quantQ.amoeba_.mpts set count p;
    `.quantQ.amoeba_.ndim set count p[0];

    `.quantQ.amoeba_.y set func each p;

    `.quantQ.amoeba_.nfunc set 0;
    `.quantQ.amoeba_.psum set sum p;

    `.quantQ.amoeba_.p set p;
    `.quantQ.amoeba_.pmin set .quantQ.amoeba_.ndim#0f;

    `.quantQ.amoeba_.fmin set 0nf;

    while[null .quantQ.amoeba_.amoeba[func];];

    `pmin`fmin`nfunc!(.quantQ.amoeba_.pmin; .quantQ.amoeba_.fmin; .quantQ.amoeba_.nfunc)
  };

// Examples
// --------

// 762.5:
// .quantQ.func.rosenbrock[1f;100f;(-1.5;-.5)];
// 756.5:
// .quantQ.func.rosenbrock[1f;100f;(1.5;-.5)];
// Global minimum (0):
// .quantQ.func.rosenbrock[1f;100f;(1f;1f)];
// .quantQ.amoeba[0.0001; (-1 1f; -1 0f;.5 0); .quantQ.func.rosenbrock[1f;100f]]
// .quantQ.amoeba[0.0001; (-1 1f; -1 0f;.25 0); .quantQ.func.rosenbrock[2.5f;100f]]
// .quantQ.amoeba[0.0001; (-1 1f; -1 0f;.5 0); {neg .quantQ.func.himmelblau[x]}]
// .quantQ.amoeba[0.0000001; (-1 1f; -1 0f; .5 0); {neg .quantQ.func.himmelblau[x]}]
