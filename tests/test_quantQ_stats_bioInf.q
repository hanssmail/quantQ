/ Unit tests for lib/quantQ_stats_bioInf.q -- sequence alignment and edit distance (.quantQ.stats)
/ Reference values from direct Python implementations (match +1, mismatch -1, gap -1).
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`stats_bioInf;

bioS1:`G`C`A`T`G`C`U;
bioS2:`G`A`T`T`A`C`A;

/ Needleman-Wunsch global alignment
.qt.approx["NeedlemanWunsch: optimal score of GCATGCU vs GATTACA is 0";{last last[.quantQ.stats.NeedlemanWunsch[bioS1;bioS2]`NWgrid] @\: 0};0f;0f];
.qt.eq["NeedlemanWunsch: all three optimal alignments are returned";{count .quantQ.stats.NeedlemanWunsch[bioS1;bioS2]`matches};3];
.qt.approx["NeedlemanWunsch: identical sequences score their length";{last last[.quantQ.stats.NeedlemanWunsch[`A`B`C;`A`B`C]`NWgrid] @\: 0};3f;0f];
.qt.eq["NeedlemanWunsch: identical sequences align to themselves";{.quantQ.stats.NeedlemanWunsch[`A`B`C;`A`B`C]`matches};enlist (`A`B`C;`A`B`C)];
.qt.approx["NeedlemanWunsch: first row and column are gap penalties";{g:.quantQ.stats.NeedlemanWunsch[`A`B;`A`C`D]`NWgrid; ((first g)@\:0;first each g@\:0)};(0 -1 -2;0 -1 -2 -3);0f];
/ unique optimal alignment ABCD- / A-CDE (Python)
.qt.eq["NeedlemanWunsch: gaps appear as nulls in the alignment";{.quantQ.stats.NeedlemanWunsch[`A`B`C`D;`A`C`D`E]`matches};enlist (`A`B`C`D`;`A``C`D`E)];
/ the traceback starts at (Nseq1;Nseq2) but the grid is indexed [row=seq2;col=seq1], so this needs equal lengths
.qt.eq["NeedlemanWunsch handles sequences of different lengths (A-C / ABC)";{.quantQ.stats.NeedlemanWunsch[`A`C;`A`B`C]`matches};enlist (`A``C;`A`B`C)];

/ Smith-Waterman local alignment
.qt.approx["SmithWaterman: best local score";{max raze .quantQ.stats.SmithWaterman[()!();`a`t`g`t`t`a`c`g`g;`z`g`g`t`t`g`a`c`t`a][`NWgrid] @\:\: 0};4f;0f];
.qt.approx["SmithWaterman: scores are never negative";{min raze .quantQ.stats.SmithWaterman[()!();`a`b`c;`d`e`f][`NWgrid] @\:\: 0};0f;0f];
.qt.eq["SmithWaterman: no alignment between unrelated sequences";{count .quantQ.stats.SmithWaterman[()!();`a`b`c;`d`e`f]`matches};0];
.qt.approx["SmithWaterman with a heavy indel penalty";{max raze .quantQ.stats.SmithWaterman[`indel`match!(3f;1f);`a`a`s`d`d;`b`b`s`v`v][`NWgrid] @\:\: 0};1f;0f];
.qt.true["SmithWaterman: the best local alignment of a shared block is found";{m:.quantQ.stats.SmithWaterman[()!();`x`a`b`c`y;`q`a`b`c`r]`matches; any {(x 0)~x 1} each m}];

/ Levenshtein distance (Wagner-Fischer)
.qt.approx["LevenshteinDistance HELLO -> WORLD";{.quantQ.stats.LevenshteinDistance[()!();`H`E`L`L`O;`W`O`R`L`D]`distance};4f;0f];
.qt.approx["LevenshteinDistance kitten -> sitting";{.quantQ.stats.LevenshteinDistance[()!();"kitten";"sitting"]`distance};3f;0f];
.qt.approx["LevenshteinDistance with substitution cost 2";{.quantQ.stats.LevenshteinDistance[enlist[`substitute]!enlist 2f;"kitten";"sitting"]`distance};5f;0f];
.qt.approx["LevenshteinDistance to itself is 0";{.quantQ.stats.LevenshteinDistance[()!();"quant";"quant"]`distance};0f;0f];
.qt.approx["LevenshteinDistance to the empty sequence is the length";{.quantQ.stats.LevenshteinDistance[()!();"abc";""]`distance};3f;0f];
.qt.eq["LevenshteinDistance returns the full matrix";{count each .quantQ.stats.LevenshteinDistance[()!();"ab";"xyz"]`distanceMatrix};3 3 3 3];

.qt.finish[];
