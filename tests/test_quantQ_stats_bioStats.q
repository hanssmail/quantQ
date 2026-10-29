/ Unit tests for lib/quantQ_stats_bioStats.q -- 2x2 contingency tables: Fisher and Barnard (.quantQ.stats)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`stats_bioStats;

bioTab:{[v] ([] grp:`a`a`b`b; feature:`A`B`A`B; val:v)};

/ Fisher: hypergeometric probability of the observed table (scipy.stats.hypergeom.pmf)
.qt.approx["FisherExactTest on a perfectly separated 3/3 table";{.quantQ.stats.FisherExactTest bioTab 3 0 0 3};0.05;1e-12];
.qt.approx["FisherExactTest on the 3 4 5 6 table";{.quantQ.stats.FisherExactTest bioTab 3 4 5 6};0.3695324283559578;1e-12];
.qt.approx["FisherExactTest does not depend on the row order of the input";{.quantQ.stats.FisherExactTest reverse bioTab 3 4 5 6};0.3695324283559578;1e-12];
.qt.approx["FisherExactTest probabilities over all tables with fixed margins sum to 1";{sum {.quantQ.stats.FisherExactTest bioTab (x;3-x;4-x;1+x)} each til 4};1f;1e-12];

/ Barnard
.qt.eq["offDiagRegion with a tiny proportion is the two extreme tables";{.quantQ.stats.offDiagRegion[3;3;0.01]};(3 0;0 3)];
.qt.true["offDiagRegion stays inside the table";{r:.quantQ.stats.offDiagRegion[4;5;0.3]; all raze (r[;0] within 0 4),r[;1] within 0 5}];
.qt.approx["BarnardTestWeight is a product of binomial probabilities";{.quantQ.stats.BarnardTestWeight[3 4;1 2;0.3]};0.11668859999999996;1e-12];
.qt.approx["BarnardTestWeight sums to 1 over all outcomes";{sum .quantQ.stats.BarnardTestWeight[3 4;;0.3] each (til 4) cross til 5};1f;1e-12];
bioBar:.qt.fixture["bioBar";{.quantQ.stats.BarnardRejectionRegion `m`n`xi`dxi`cl`calcCl`pGrid`run!(3;3;0.0;0.005;0.1;0.0;200;1b)}];
.qt.true["BarnardRejectionRegion: size of the region is below the confidence level";{bioBar[`calcCl]<0.1}];
.qt.true["BarnardRejectionRegion contains the most extreme tables";{r:bioBar`rejectionRegion; all {[r;p] any p~/:r}[r] each (3 0;0 3)}];
.qt.true["BarnardRejectionRegion: the reported region has size below cl for every p";{r:bioBar`rejectionRegion; 0.1>max {[r;p] sum .quantQ.stats.BarnardTestWeight[3 3;;p] each r}[r] each 0.01*til 101}];

.qt.finish[];
