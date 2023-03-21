// classical technical analysis


// The functions (user interface) in this repository follow general structure:
// .quantQ.ta.f[sourceColumns;params;tab]
// columns -- string or list of strings, ordered, names of source columns
// params -- dictionary with parameters, ()!() always acceptable producing default setup
// tab -- source table, which contains columns and which is updated

// using .quantQ.stats

//////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////
// List of implemented TIs

// simple moving average: .quantQ.ta.ma

// exponential moving average: .quantQ.ta.ema

// simple moving standard deviation: .quantQ.ta.msd

// the Bollinger Bands: .quantQ.ta.bollingerBands

// momentum indicator .quantQ.ta.momentum  

// Acceleration bands .quantQ.ta.accBands

// Ichimoku indicator .quantQ.ta.Ichimoku 

// MACD .quantQ.ta.MACD

// Aroon Up/Down/Oscillator .quantQ.ta.Aroon

// Relative Strength Indicator .quantQ.ta.RSI 

// Average True Range .quantQ.ta.ATR
    
//////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////
// Functions 

// simple moving average
.quantQ.ta.ma:{[inp;params;tab]        
    // inp -- name of the source column
    // params -- parameters
    // tab -- table
    params:(enlist[`memory]!enlist[10]),params;
    :![tab; (); 0b; 
    enlist[`$ string[inp],"MA",string[params[`memory]]]!enlist[(mavg;params[`memory];inp)]];
 };

// exponential moving average
.quantQ.ta.ema:{[inp;params;tab]        
    // inp -- name of the source column
    // params -- parameters
    // tab -- table
    params:(enlist[`memory]!enlist[10]),params;
    :![tab; (); 0b; 
    enlist[`$ string[inp],"EMA",string[params[`memory]]]!enlist[(.quantQ.stats.expma1[2.0%(params[`memory]+1)];inp)]];
 };

// simple moving standard deviation
.quantQ.ta.msd:{[inp;params;tab]        
    // inp -- name of the source column
    // params -- parameters
    // tab -- table
    params:(enlist[`memory]!enlist[10]),params;
    :![tab; (); 0b; 
    enlist[`$ string[inp],"MSD",string[params[`memory]]]!enlist[(sqrt;(mdev;params[`memory];inp))]];
 };

// the Bollinger Bands
.quantQ.ta.bollingerBands:{[inp;params;tab]        
    // inp -- name of the source column
    // params -- parameters
    // tab -- table
    params:(enlist[`memory]!enlist[10]),params;
    // temporaray tab
    temp: ([] x: tab[inp]);
    header: {`$string[x],string[z],string[y]}[inp;params[`memory];] each `MiddleBound`LowerBound`UpperBound;
    :tab,'header xcol select xM, xM-2*xS, xM+2*xS from `x`xM`xS xcol .quantQ.ta.msd[`x;params;] .quantQ.ta.ma[`x;params;temp];
 };

// momentum indicator
.quantQ.ta.momentum:{[inp;params;tab]        
    // inp -- name of the source column
    // params -- parameters
    // tab -- table
    params:(enlist[`memory]!enlist[10]),params;
    :![tab; (); 0b;enlist[`$ string[inp],"Mom",string[params[`memory]]]!enlist[(-;inp;({ raze(x#y;(neg[x] _ y))  };params[`memory];inp))]];
 };

// Acceleration bands
.quantQ.ta.accBands:{[inp;params;tab]        
    // inp -- ordered names of the OHLC column 
    // params -- parameters
    // tab -- tab
    // name is derived from "close" price column
    params:(enlist[`memory]!enlist[10]),params;   
    :![tab; (); 0b;(`$ string[inp[3]],"UpperBand",string[params[`memory]];
    `$ string[inp[3]],"MiddleBand",string[params[`memory]];
    `$ string[inp[3]],"LowerBand",string[params[`memory]])!
    ((mavg;params[`memory];
    (*;inp[1];
    (+;1.0;(*;4.0;(%;
    (-;inp[1];inp[2]);
    (+;inp[1];inp[2])
    )))));
    (mavg;params[`memory];inp[3]);
    (mavg;params[`memory];
    (*;inp[2];
    (-;1.0;(*;4.0;(%;
    (-;inp[1];inp[2]);
    (+;inp[1];inp[2])
    ))))))];
 };

// Ichimoku indicator
.quantQ.ta.Ichimoku:{[inp;params;tab]        
    // inp -- ordered names of the OHLC column 
    // params -- parameters
    // tab -- table
    params:((`memoryShort`memoryMedium`memoryLong)!(9;26;52)),params;
    tab:![tab; (); 0b;
    (`$ string[inp[3]],"TurningLine",string[params[`memoryShort]],"x",string[params[`memoryMedium]],"x",string[params[`memoryLong]];
    `$ string[inp[3]],"StandardLine",string[params[`memoryShort]],"x",string[params[`memoryMedium]],"x",string[params[`memoryLong]])!
    ((%;(+;(mmax;params[`memoryShort];inp[1]);(mmin;params[`memoryShort];inp[2]));2.0);
    (%;(+;(mmax;params[`memoryMedium];inp[1]);(mmin;params[`memoryMedium];inp[2]));2.0))];
    :![tab; (); 0b;(`$ string[inp[3]],"LeadingSpan1x",string[params[`memoryShort]],"x",string[params[`memoryMedium]],"x",string[params[`memoryLong]];
    `$ string[inp[3]],"LeadingSpan2x",string[params[`memoryShort]],"x",string[params[`memoryMedium]],"x",string[params[`memoryLong]])!
    (({ raze(x#0.0;(neg[x] _ y))  };params[`memoryMedium];
    (%;(+;`$ string[inp[3]],"TurningLine",string[params[`memoryShort]],"x",string[params[`memoryMedium]],"x",string[params[`memoryLong]];
    `$ string[inp[3]],"StandardLine",string[params[`memoryShort]],"x",string[params[`memoryMedium]],"x",string[params[`memoryLong]]
    );2.0));
    ({ raze(x#0.0;(neg[x] _ y))  };params[`memoryMedium];
    (%;(+;(mmax;params[`memoryLong];inp[1]);(mmin;params[`memoryLong];inp[2]));2.0)))];
 };
 
// MACD
.quantQ.ta.MACD:{[inp;params;tab]        
    // inp -- name of the source column 
    // params -- parameters
    // tab -- table
    params:((`memoryFast`memorySlow`memory)!(10;20;15)),params; 
    :![tab; (); 0b;(`$ string[inp],"MACDHistogram",string[params[`memoryFast]],"x"
        ,string[params[`memorySlow]],"x",string[params[`memory]];`$ string[inp],"MACDSignalLine",string[params[`memoryFast]],"x"
        ,string[params[`memorySlow]],"x",string[params[`memory]];
        `$ string[inp],"MACD",string[params[`memoryFast]],"x"
        ,string[params[`memorySlow]]
        )!(
        (-;(-;(mavg;params[`memoryFast];inp);(mavg;params[`memorySlow];inp));
        (mavg;params[`memory];(-;(mavg;params[`memoryFast];inp);
        (mavg;params[`memorySlow];inp))));
        (mavg;params[`memory];(-;(mavg;params[`memoryFast];inp);
        (mavg;params[`memorySlow];inp)));
        (-;(mavg;params[`memoryFast];inp);
        (mavg;params[`memorySlow];inp))
        )];
 };

// Aroon Up/Down/Oscillator
.quantQ.ta.Aroon:{[inp;params;tab]        
    // inp -- name of the source column 
    // params -- parameters
    // tab -- table
    params:(enlist[`memory]!enlist[20]),params; 
    // temp function
    f:({raze [(0nf; neg[1] _ x)]}\)[params[`memory];];
    tab:![tab; (); 0b;
        (`$ string[inp],"AroonUp",string[params[`memory]];
        `$ string[inp],"AroonDown",string[params[`memory]])!
        ((%;(-;(params[`memory]-1);
        (each;{first where 1=x};
        (=;(mmax;params[`memory];inp);
        (flip;(f;inp)))));params[`memory]-1);
        (%;(-;(params[`memory]-1);
        (each;{first where 1=x};
        (=;(mmin;params[`memory];inp);
        (flip;(f;inp)))));params[`memory]-1)
    )];
    // add oscillator
    :![tab;();0b;
        enlist[`$ string[inp],"AroonOscillator",string[params[`memory]]]!enlist[
        (-;`$ string[inp],"AroonUp",string[params[`memory]];
        `$ string[inp],"AroonDown",string[params[`memory]])]
    ];
 };

// Relative Strength Indicator
.quantQ.ta.RSI:{[inp;params;tab]        
    // inp -- ordered names of the open and close columns
    // params -- parameters
    // tab -- table
    params:(enlist[`memory]!enlist[20]),params; 
    // temp function
    f:({raze [(0nf; neg[1] _ x)]}\)[params[`memory];];
    g:{
        p:0^avg[x where x>=0];
        n:0^abs[avg[x where x<0]];
        :$[(p=0) and (n=0);1.0;$[(n=0);0wf;p%n]];
    };
    :![tab;();0b;enlist[`$ string[inp[1]],"RSI",string[params[`memory]]]!enlist[ 
        (-;100.0;(%;100.0;(+;1;(each;g;(flip;(f;(-;inp[1];inp[0])))))))
    ]];
 };

// Average True Range
.quantQ.ta.ATR:{[inp;params;tab]        
    // inp -- ordered names of the open/high/low/close columns
    // params -- parameters
    // tab -- table
    params:(enlist[`memory]!enlist[20]),params; 
    // temp function
    f:{max[(x;y;z)]};
    :![tab;();0b;enlist[`$ string[inp[3]],"ATR",string[params[`memory]]]!enlist[(mavg;params[`memory];(f;(-;inp[1];inp[2]);(abs;(-;inp[1];inp[3]));(abs;(-;inp[2];inp[3]))))]];
 };
 
 
//////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////
// Examples
// tab:([] price:100?1.0);
// we ensure that OHLC is >0
// tabOHLC: update high: {max[(x;y)]}'[high;close],low: {min[(x;y)]}'[low;close] from select open: price, high: price+100?0.3, low: price-100?0.3, close:price+neg[0.1]+100?0.2 from update price:1+price from tab;

// .quantQ.ta.ma[`price;()!()] tab

// .quantQ.ta.ema[`price;()!()] tab

// .quantQ.ta.msd[`price;()!()] tab

// meta .quantQ.ta.bollingerBands[inp;()!()] tab 

// meta .quantQ.ta.momentum[`price;()!()] tab

// meta .quantQ.ta.accBands[`open`high`low`close;()!()] tabOHLC

// meta .quantQ.ta.Ichimoku[`open`high`low`close;()!()] tabOHLC 
  
// meta .quantQ.ta.MACD[`price;()!()] tab  

// meta .quantQ.ta.Aroon[`price;()!()] tab  

// meta .quantQ.ta.RSI[`open`close;()!()] tabOHLC

// meta .quantQ.ta.ATR[`open`high`low`close;()!()] tabOHLC