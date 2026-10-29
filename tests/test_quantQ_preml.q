/ Unit tests for lib/quantQ_preml.q -- feature preprocessing (.quantQ.preml)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`preml;

.qt.eq["alpha2Int maps sorted distinct values to 0..n-1";{.quantQ.preml.alpha2Int `c`a`b`a};`a`b`c!0 1 2];
.qt.eq["distanceCategorical is an equality mask";{.quantQ.preml.distanceCategorical[`a;`a`b`a]};101b];
.qt.approx["zScoreNorm uses the sample standard deviation";{.quantQ.preml.zScoreNorm 1 2 3f};-1 0 1f;1e-12];
.qt.approx["zScoreNorm: result has mean 0";{avg .quantQ.preml.zScoreNorm 3 7 1 9 4f};0f;1e-12];
.qt.approx["zScoreNorm: result has sample sd 1";{sdev .quantQ.preml.zScoreNorm 3 7 1 9 4f};1f;1e-12];
.qt.approx["zScoreNormExt uses the supplied mean and sd";{.quantQ.preml.zScoreNormExt[1 2 3f;`mean`std!2 0.5]};-2 0 2f;1e-12];
.qt.approx["rangeScaleNorm maps [min,max] onto [a,b]";{.quantQ.preml.rangeScaleNorm[1 2 3f;0 10f]};0 5 10f;1e-12];
.qt.approx["rangeScaleNorm onto a negative range";{.quantQ.preml.rangeScaleNorm[2 4 6f;-1 1f]};-1 0 1f;1e-12];
.qt.approx["rangeScaleNormExt uses the supplied min and max";{.quantQ.preml.rangeScaleNormExt[1 2 3f;0 1f;`min`max!0 4f]};0.25 0.5 0.75;1e-12];

.qt.finish[];
