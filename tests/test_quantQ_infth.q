/ Unit tests for lib/quantQ_infth.q -- information theory: word statistics, tf-idf, protein kernels (.quantQ.infth)
if[not `qt in key `; system "l ",{$[count w:where x="/";(last w)#x;"."]}[string .z.f],"/quantQ_test.q"];
.qt.setModule`infth;

infCorpus:("the cat sat";"the dog sat down";"a cat";"the cat ate the fish");

/ word frequencies
.qt.eq["wordFrequency counts words across the corpus";{.quantQ.infth.wordFrequency ("a b";"b c")};(enlist "a";enlist "b";enlist "c")!1 2 1];
.qt.eq["wordFrequency on the test corpus: 'the'";{.quantQ.infth.wordFrequency[infCorpus]"the"};4];
.qt.approx["wordMeasure: in-document count / corpus count";{first .quantQ.infth.wordMeasure[enlist "a b";.quantQ.infth.wordFrequency ("a b";"a a c")]};(enlist "a";enlist "b")!(1%3;1f);1e-12];
.qt.approx["wordAddAnalytics adds AVERAGE and MAX";{d:first .quantQ.infth.wordAddAnalytics enlist (enlist "a";enlist "b")!0.25 0.75; d("AVERAGE";"MAX")};0.5 0.75;1e-12];

/ idf = log(N / number of documents containing the word)
.qt.approx["idf of 'the' (3 of 4 documents)";{.quantQ.infth.idf[infCorpus]"the"};log 4%3;1e-12];
.qt.approx["idf of 'fish' (1 of 4 documents)";{.quantQ.infth.idf[infCorpus]"fish"};log 4f;1e-12];
.qt.approx["idf of a word present everywhere is 0";{.quantQ.infth.idf[("x y";"y z")]"y"};0f;1e-12];
.qt.approx["tfidfSentence = term count * idf";{t:.quantQ.infth.tfidfSentence["the cat ate the fish";.quantQ.infth.idf infCorpus]; t("the";"fish")};(2*log 4%3;log 4f);1e-12];
.qt.eq["findBestSentence picks the sentence with the highest tf-idf";{.quantQ.infth.findBestSentence[("a b";"a a c";"c d");"a"]};"a a c"];

/ char <-> int
.qt.eq["char2int";{.quantQ.infth.char2int "A"};65i];
.qt.eq["int2char";{.quantQ.infth.int2char 66};"B"];
.qt.eq["char2int/int2char round trip";{.quantQ.infth.int2char .quantQ.infth.char2int "quantQ"};"quantQ"];

/ proteins
.qt.eq["getProtein length within bounds (fixed length)";{count .quantQ.infth.getProtein[7;7;4]};7];
.qt.true["getProtein uses the first nAcids letters";{all .quantQ.infth.getProtein[5;30;3] in "ABC"}];
.qt.eq["proteinGetSubstring, N=3";{.quantQ.infth.proteinGetSubstring["ABCDE";3]};("ABC";"BCD";"CDE")];
.qt.eq["proteinGetSubstring honours N (N=2)";{.quantQ.infth.proteinGetSubstring["ABCDE";2]};("AB";"BC";"CD";"DE")];
.qt.eq["proteinGetFeatureVector counts substrings";{.quantQ.infth.proteinGetFeatureVector ("ABC";"BCD";"ABC")};("ABC";"BCD")!2 1];
.qt.approx["proteinInnerProduct of a protein with itself";{.quantQ.infth.proteinInnerProduct["ABCABC";"ABCABC";3]};6f;1e-12];
.qt.approx["proteinInnerProduct of disjoint proteins is 0";{.quantQ.infth.proteinInnerProduct["AAAA";"BBBB";3]};0f;1e-12];
.qt.approx["proteinDistance to itself is 0";{.quantQ.infth.proteinDistance["ABCD";"ABCD";3]};0f;1e-12];
.qt.approx["proteinDistance ABCD vs BCDA";{.quantQ.infth.proteinDistance["ABCD";"BCDA";3]};2f;1e-12];

.qt.finish[];
