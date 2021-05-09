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


//////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////
// Examples
// tab:([] ret:100?1.0);

// .quantQ.ta.ma[`ret;()!()] tab

// .quantQ.ta.ema[`ret;()!()] tab


