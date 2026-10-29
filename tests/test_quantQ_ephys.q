/ Unit tests for lib/quantQ_ephys.q -- econophysics: power laws and the Ising market model (.quantQ.ephys)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`ephys;

/ powerLaw
.qt.approx["powerLaw alpha*x^-beta";{.quantQ.ephys.powerLaw[2;2;4]};0.125;1e-12];
.qt.approx["powerLaw is vectorised";{.quantQ.ephys.powerLaw[1;1;1 2 4f]};1 0.5 0.25;1e-12];

/ getSquareAround -- 8 neighbours with periodic boundaries
.qt.eq["getSquareAround: interior point";{asc .quantQ.ephys.getSquareAround[5;2 2]};asc (1 1;1 2;1 3;2 1;2 3;3 1;3 2;3 3)];
.qt.eq["getSquareAround: corner wraps around";{asc .quantQ.ephys.getSquareAround[5;0 0]};asc (4 4;4 0;4 1;0 4;0 1;1 4;1 0;1 1)];
.qt.true["getSquareAround never returns the point itself";{not any (0 0)~/:.quantQ.ephys.getSquareAround[5;0 0]}];

/ probPositiveSpin -- logistic in 2*beta*energy
.qt.approx["probPositiveSpin is 1/2 at zero energy";{.quantQ.ephys.probPositiveSpin[1;0]};0.5;1e-12];
.qt.approx["probPositiveSpin at beta=1, E=1";{.quantQ.ephys.probPositiveSpin[1;1]};0.8807970779778823;1e-12];
.qt.approx["probPositiveSpin(E) + probPositiveSpin(-E) = 1";{.quantQ.ephys.probPositiveSpin[0.7;1.3]+.quantQ.ephys.probPositiveSpin[0.7;-1.3]};1f;1e-12];

/ entropy term -p log p
.qt.approx["entropy(0) = 0";{.quantQ.ephys.entropy 0};0f;1e-12];
.qt.approx["entropy(1) = 0";{.quantQ.ephys.entropy 1};0f;1e-12];
.qt.approx["entropy(1/2) = log(2)/2";{.quantQ.ephys.entropy 0.5};0.34657359027997264;1e-12];

/ Ising market: one Monte-Carlo update
isingN:4;
isingModel:{[priceType]
    `nGrid`priceIsing`wallStreet`getSquareAround`alpha`beta`probPositiveSpin!(
        isingN;
        ([] t:enlist 0; price:enlist priceType$isingN*isingN);
        ([] coord:(til isingN) cross til isingN; spin:(isingN*isingN)#1);
        .quantQ.ephys.getSquareAround;
        1.0;
        0.5;
        .quantQ.ephys.probPositiveSpin)
 };
.qt.eq["oneRunIsing appends one price";{count (.quantQ.ephys.oneRunIsing isingModel[`long])`priceIsing};2];
.qt.true["oneRunIsing: price is the sum of the spins";{m:.quantQ.ephys.oneRunIsing isingModel[`long]; (last m[`priceIsing]`price)=sum m[`wallStreet]`spin}];
.qt.true["oneRunIsing: at most one spin changes";{m:.quantQ.ephys.oneRunIsing isingModel[`long]; 1>=sum m[`wallStreet][`spin]<>16#1}];
.qt.true["oneRunIsing: spins stay in {-1,1}";{m:(.quantQ.ephys.oneRunIsing/)[50;isingModel[`long]]; all m[`wallStreet][`spin] in -1 1}];
.qt.true["oneRunIsing_v2: price is the average spin";{m:.quantQ.ephys.oneRunIsing_v2 isingModel[`float]; (last m[`priceIsing]`price)=avg m[`wallStreet]`spin}];
.qt.eq["oneRunIsing_v2: time advances by one";{last (.quantQ.ephys.oneRunIsing_v2 isingModel[`float])[`priceIsing]`t};1];

.qt.finish[];
