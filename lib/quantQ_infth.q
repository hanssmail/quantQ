.quantQ.infth.wordFrequency:{[corpus]
    // corpus -- corpus with cleaned text
    :count each group" "vs" "sv corpus;
 };

.quantQ.infth.wordMeasure:{[document;wordFrequency]
    // document -- text we want to analyse
    // wordFrequency -- dictionary with corpus word frequency
    // calculates word frequency measure
    :distinct {[x;y] tmp:count each group" "vs x;(tmp% y key tmp)}[;wordFrequency]each document;
 };

.quantQ.infth.wordAddAnalytics:{[documentMeasure]
    // documentMeasure -- analysed document with defined measure
    // analyse sentence by sentence and add average and max measure
    :{tmpAvg: avg value x; tmpMax: max value x;  (raze (key x;enlist "AVERAGE";enlist "MAX"))!
        (raze (value x;enlist tmpAvg;enlist tmpMax))} each documentMeasure;
 };

.quantQ.infth.idf:{[corpus]
    // corpus -- corpus with cleaned text
    uniqueWords: distinct " "vs" "sv corpus;
    // count of documents with a given word
    vals: {[corpus;word] sum { x in " "vs y}[word;] each corpus}[corpus;] each uniqueWords;
    // idf
    :log uniqueWords!(count corpus)%vals;
 };

.quantQ.infth.tfidfSentence:{[sentence;idf]
    // sentence -- sentence of words to be analysed   
    // idf -- idf calculated for the corpus 
    // word frequency dictionary   
    wordFrequency: .quantQ.infth.wordFrequency[enlist sentence];
    // tf-idf for every word
    :t!(value wordFrequency)*idf each t:key wordFrequency;
 };

.quantQ.infth.findBestSentence:{[corpus;word]
    // corpus -- corpus with cleaned text
    // word -- word to be searched for
    // subset of corpus with given word
    corpusSubset:corpus where { first x in " "vs y}[word;] each corpus;
    // return the sentence with the maximum tf-idf for a given word, if more then one, select randomly
    :first 1?corpusSubset where max[t]=t:{first x[y]}[;word] each
        .quantQ.infth.tfidfSentence[;.quantQ.infth.idf[corpus]] each corpusSubset;
 };

.quantQ.infth.char2int:{[char] 
    // char -- argument of type char
    // conversion from char to ASCII
    :"i" $ char;
 };

.quantQ.infth.int2char:{[int] 
    // int -- argument of type integer
    // conversion from ASCII to char
    :"c" $ int;
 };

.quantQ.infth.getProtein:{[minLength;maxLength;nAcids]
    // minLength -- min length of the protein
    // maxLength -- max length of the protein
    // nAcids -- number of different amino acids, number should not exceed 26 to have alphabetic representation
    // capital "A" offset for ASCII
    offset:.quantQ.infth.char2int["A"];
    // choose the length of the protein    
    proteinLength: first minLength+1?(1+maxLength-minLength);
    // create ASCII protein
    protein: offset+proteinLength?nAcids;
    // convert to chars and export
    :.quantQ.infth.int2char[protein];
 };

.quantQ.infth.proteinGetSubstring:{[p;N]
    // p -- protein 
    // N -- length of substring
    // return all substrings of length N
    :({[dict] dict[`strings]:dict[`strings],enlist 3#dict[`protein]; dict[`protein]:1_dict[`protein];dict }/)
        [{[x;N] N<=count x[`protein]}[;N];`strings`protein!(();p)][`strings];
 };

.quantQ.infth.proteinGetFeatureVector:{[substrings]
    // substrings -- outcome of .quantQ.infth.proteinGetSubstring function
    // returns the vector of features as dictionary of substrings and its occurrence
    :(#:)'[group substrings];
 };

.quantQ.infth.proteinInnerProduct:{[p1;p2;N]
    // p1 -- protein 1
    // p2 -- protein 2
    // N -- length of substring
    // calculate vector of features for both proteins   
    d1:.quantQ.infth.proteinGetFeatureVector[.quantQ.infth.proteinGetSubstring[p1;3]];
    d2:.quantQ.infth.proteinGetFeatureVector[.quantQ.infth.proteinGetSubstring[p2;3]];
    // return inner product
    :exec (0f^d1)$0f^d2 from ([substring: key d1];value d1) uj ([substring: key d2];value d2);
 };

.quantQ.infth.proteinDistance:{[p1;p2;N]
    // p1 -- protein 1
    // p2 -- protein 2
    // N -- length of substring
    // return measure
    :.quantQ.infth.proteinInnerProduct[p1;p1;N]+.quantQ.infth.proteinInnerProduct[p2;p2;N]-
        2.0*.quantQ.infth.proteinInnerProduct[p1;p2;N];
 };
