// Rosenbrock's valley (banana) function:
.quantQ.opt.func.rosenbrock: {[a;b;x] t1: a-x[0]; t2: x[1]-x[0]*x[0]; (t1*t1) + b*t2*t2};

.quantQ.opt.func.himmelblau: {[x] t1: (x[0]*x[0]) + x[1] - 11f; t2: x[0] + (x[1]*x[1]) - 7f; (t1*t1) + t2*t2};

.quantQ.opt.amoeba_.amotry: {[ihi;fac;func]
    fac1: (1f-fac)%.quantQ.opt.amoeba_.ndim;
    fac2: fac1-fac;
    ptry: (.quantQ.opt.amoeba_.psum*fac1) - (.quantQ.opt.amoeba_.p[ihi]*fac2);
    ytry: func[ptry];
    $[ytry < .quantQ.opt.amoeba_.y[ihi]; [
        .quantQ.opt.amoeba_.y[ihi]: ytry;
        .quantQ.opt.amoeba_.psum: .quantQ.opt.amoeba_.psum + (ptry - .quantQ.opt.amoeba_.p[ihi]);
        .quantQ.opt.amoeba_.p[ihi]: ptry;
    ];];
    ytry
  };

.quantQ.opt.amoeba_.amoeba: {[func]
    ilo: 0;
    $[.quantQ.opt.amoeba_.y[0]>.quantQ.opt.amoeba_.y[1]; [ihi:0; inhi:1]; [ihi:1; inhi:0]];
    i: 0;
    while[i<.quantQ.opt.amoeba_.mpts;
        $[.quantQ.opt.amoeba_.y[i]<=.quantQ.opt.amoeba_.y[ilo]; ilo:i;];
        $[.quantQ.opt.amoeba_.y[i]>.quantQ.opt.amoeba_.y[ihi]; [inhi:ihi; ihi:i]; $[and[.quantQ.opt.amoeba_.y[i]>.quantQ.opt.amoeba_.y[inhi]; i<>ihi]; inhi:i;]];
        i: i+1;
    ];

    rtol: (2f*abs[.quantQ.opt.amoeba_.y[ihi]-.quantQ.opt.amoeba_.y[ilo]]) % (abs[.quantQ.opt.amoeba_.y[ihi]]+abs[.quantQ.opt.amoeba_.y[ilo]]+.quantQ.opt.amoeba_.TINY);

    $[rtol<.quantQ.opt.amoeba_.ftol; [
        temp: .quantQ.opt.amoeba_.y[0]; .quantQ.opt.amoeba_.y[0]: .quantQ.opt.amoeba_.y[ilo]; .quantQ.opt.amoeba_.y[ilo]: temp;
        i: 0;
        while[i<.quantQ.opt.amoeba_.ndim;
            temp: .quantQ.opt.amoeba_.p[0;i]; .quantQ.opt.amoeba_.p[0;i]: .quantQ.opt.amoeba_.p[ilo;i]; .quantQ.opt.amoeba_.p[ilo;i]: temp;
            .quantQ.opt.amoeba_.pmin[i]: .quantQ.opt.amoeba_.p[0;i];
            i: i+1;
        ];
        `.quantQ.opt.amoeba_.fmin set .quantQ.opt.amoeba_.y[0];
        :.quantQ.opt.amoeba_.fmin;
    ];];

    $[.quantQ.opt.amoeba_.nfunc>=.quantQ.opt.amoeba_.NMAX; '`MAX_EVALUATIONS_EXCEEDED;];
    `.quantQ.opt.amoeba_.nfunc set .quantQ.opt.amoeba_.nfunc+2;
    ytry: .quantQ.opt.amoeba_.amotry[ihi;-1f;func];
    $[ytry<=.quantQ.opt.amoeba_.y[ilo]; [
        ytry: .quantQ.opt.amoeba_.amotry[ihi;2f;func];
    ]; [ $[ytry>=.quantQ.opt.amoeba_.y[inhi]; [
        ysave: .quantQ.opt.amoeba_.y[ihi];
        ytry: .quantQ.opt.amoeba_.amotry[ihi;.5;func];
        $[ytry>=ysave; [
            i: 0;
            while[i<.quantQ.opt.amoeba_.mpts;
                $[i <> ilo; [
                    `.quantQ.opt.amoeba_.psum set .5*(.quantQ.opt.amoeba_.p[i] + .quantQ.opt.amoeba_.p[ilo]);
                    .quantQ.opt.amoeba_.p[i]: .quantQ.opt.amoeba_.psum;
                    .quantQ.opt.amoeba_.y[i]: func[.quantQ.opt.amoeba_.psum];
                ];];
                i: i+1;
            ];
            `.quantQ.opt.amoeba_.nfunc set .quantQ.opt.amoeba_.nfunc + .quantQ.opt.amoeba_.ndim;
            `.quantQ.opt.amoeba_.psum set sum .quantQ.opt.amoeba_.p;
        ];];
    ]; [`.quantQ.opt.amoeba_.nfunc set .quantQ.opt.amoeba_.nfunc-1;] ]; ] ];
  };

.quantQ.opt.amoeba: {[ftol;p;func]
    `.quantQ.opt.amoeba_.ftol set ftol;
    `.quantQ.opt.amoeba_.NMAX set 5000;
    `.quantQ.opt.amoeba_.TINY set 1.0e-10;
    `.quantQ.opt.amoeba_.mpts set count p;
    `.quantQ.opt.amoeba_.ndim set count p[0];

    `.quantQ.opt.amoeba_.y set func each p;

    `.quantQ.opt.amoeba_.nfunc set 0;
    `.quantQ.opt.amoeba_.psum set sum p;

    `.quantQ.opt.amoeba_.p set p;
    `.quantQ.opt.amoeba_.pmin set .quantQ.opt.amoeba_.ndim#0f;

    `.quantQ.opt.amoeba_.fmin set 0nf;

    while[null .quantQ.opt.amoeba_.amoeba[func];];

    `pmin`fmin`nfunc!(.quantQ.opt.amoeba_.pmin; .quantQ.opt.amoeba_.fmin; .quantQ.opt.amoeba_.nfunc)
  };

// Examples
// --------

// 762.5:
// .quantQ.opt.func.rosenbrock[1f;100f;(-1.5;-.5)]
// 756.5:
// .quantQ.opt.func.rosenbrock[1f;100f;(1.5;-.5)]
// Global minimum (0):
// .quantQ.opt.func.rosenbrock[1f;100f;(1f;1f)]
// .quantQ.opt.amoeba[0.0001; (-1 1f; -1 0f;.5 0); .quantQ.opt.func.rosenbrock[1f;100f]]
// .quantQ.opt.amoeba[0.0001; (-1 1f; -1 0f;.25 0); .quantQ.opt.func.rosenbrock[2.5f;100f]]
// .quantQ.opt.amoeba[0.0001; (-1 1f; -1 0f;.5 0); {neg .quantQ.opt.func.himmelblau[x]}]
// .quantQ.opt.amoeba[0.0000001; (-1 1f; -1 0f; .5 0); {neg .quantQ.opt.func.himmelblau[x]}]
