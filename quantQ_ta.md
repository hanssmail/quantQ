# Technical Analysis

The library ```.quantQ.ta``` comprises many standard technical indicators used in the Technical Analysis. The provided set of functions aims to help users with feature engineering and analyse the data. For the sake of simplicity, it is the user's responsibility to have data sorted by the correct time index. Time index is used only when necessary, e.g., the time distance between observations is needed.

The technical analysis library is using functions from the ```.quantQ``` space.

## General Convention

The functions in the library follow a general convention. Each function, which is supposed to be used by the user, i.e., not utility functions, is defined as:

```
.quantQ.ta.f[sourceColumns;params;tab]
```
where

* ```columns``` is a string or list of strings, ordered, with names of source columns used for the calculation,

* ```params``` is a dictionary with parameters,  ```()!()``` is always acceptable such that the function is using default parameters

* ```tab``` is a source table, which contains columns and which will be  updated by the function with new columns

The function adds one or more columns, depending on the indicator. Whenever possible, the parameters are part of the column name, .i.e., the simple moving average technical indicator with 10 step memory applied on a column named ```price``` update the table with a new column called ```priceMA10```. When the technical indicator is based on the Open/High/Low/Close columns, we derive the name from the ```close``` column.


This approach allows users to add a number of the same technical indicators with different parameters, and the process is down the line (either by visualising the new columns, or, for instance, using a more quantitative approach with the Lasso regression, as implemented in ```.quantQ.reg.LAR```).

The library assumes the input is the price. Some indicators expect a single price per row, while the other indicators expect the price quadruplet of the open, high, low and close prices. The primary purpose of the indicators is either to produce a quantity, which is of the same level as the price or to create a standardised indicator normalised to some predefined scale.   

## Implemented Functions and Examples

We present a primary use of the provided functions.

### Synthetic Data Generation

We use several tables with synthetically generated data to illustrate functionality.

The simple table with a single dummy column called ```price```, which represents the price of an asset:

```
tab:([] price:100?1.0);
```

Another table we use has four columns with Open/High/Low/Close columns. This is frequently used in technical analysis. Every record in the table corresponds to a specific time interval. Then, the Open/High/Low/Close corresponds to first, max, min, last price within the time interval. This is frequently visualised as a candle bar. Our synthetic table is created as:

```
tabOHLC: update high: {max[(x;y)]}'[high;close],low: {min[(x;y)]}'[low;close] from select open: price, high: price+100?0.3, low: price-100?0.3, close:price+neg[0.1]+100?0.2 from update price:1+price from tab;
```
where we have ensured that the OHLC values are ```>0```.


### Simple Moving Average

The simple moving average is provided as

```
.quantQ.ta.ma[`price;enlist[`memory]!enlist[20]] tab
```
with default values of the ```memory``` being set to ```10```.

The table is updated as:

| c	| t	|
|  ----- |  ---- |
| price | f |
| priceMA20	| f |

### Exponential Moving Average

The exponential moving average is provided as

```
.quantQ.ta.ema[`price;enlist[`memory]!enlist[20]] tab
```
with default values of the ```memory``` being set to ```10```. The exponential moving average is updated as:

```
emaNew = column * (2.0%(memory+1))+ emaOld * (1-2.0%(memory+1))
```

The function updates the table as:

| c	| t	|
|  ----- |  ---- |
| price | f |
| priceEMA20	| f |


### Simple Moving Standard Deviation:

The simple moving standard deviation is provided as

```
.quantQ.ta.msd[`price;enlist[`memory]!enlist[20]] tab
```
with default values of the ```memory``` being set to ```10```. The function calculates the standard deviation for the last ```memory``` ticks.

The function updates the table as:

| c	| t	|
|  ----- |  ---- |
| price | f |
| priceMSD20	| f |

### Bollinger Bands

The Bollinger bands are defined as a moving average with a confidence interval:

```
Middle Band = memory mavg price
```

and confidence bands

```
Upper Band = Middle Band + (2 * Moving Standard Deviation)
```
```
Lower Band = Middle Band - (2 * Moving Standard Deviation)
```
with ```MSD``` being the simple moving standard deviation.

The function is provided as:

```
.quantQ.ta.bollingerBands[`price;enlist[`memory]!enlist[20]] tab
```
with default values of the ```memory``` being set to ```10```.

The function updates the table as:

| c	| t	|
|  ----- |  ---- |
| price	| f	|		
| priceMiddleBound10	| f	|
| priceLowerBound10	| f	|		
| priceUpperBound10	| f	|		


### Momentum Indicator

The momentum indicator is provided as:

```
.quantQ.ta.momentum[`price;enlist[`memory]!enlist[20]] tab
```
with default values of the ```memory``` being set to ```10```. The momentum is defined as:

```
MOM =  price - price[i-memory]
```

The function updates the table as:

| c	| t	|
|  ----- |  ---- |
| price	| f |		
| priceMom20	| f |		

### Acceleration Bands

The accelerator bands are provided as:

```
.quantQ.ta.accBands[`open`high`low`clos;enlist[`memory]!enlist[20]] tabOHLC
```
with default values of the ```memory``` being set to ```10```.

The accelerator bands is a moving average with error bands, where the error is defined with standard deviation estimated using OHLC table:

```
Upper Band = Middle Band + (mHigh * ( 1 + 4 * (mHigh - mLow) / (mHigh + mLow)))
```

```
Middle Band = memory mavg close
```

```
Lower Band = Simple Moving Average (Low * (1 - 4 * (High - Low)/ (High + Low)))
```

where
```
mHigh = memory mmax high
```
and
```
mLow = memory mmin low
```

The function updates the table as:

| c	| t	|
|  ----- |  ---- |
| open	| f |		
| high	| f |		
| low	| f |		
| close	| f |		
| closeUpperBand10	| f |		
| closeMiddleBand10	| f |		
| closeLowerBand10	| f |		


### Ichimoku indicator

The Ichimoku indicator can be used as:

```
 .quantQ.ta.Ichimoku[`open`high`low`close;()!()] tabOHLC
 ```
with default values of the ```memoryShort``` being 9, the ```memoryMedium``` being 26, and the ```memoryLong``` being 52.

The Ichimoku indicator is a combination of several moving averages using the extreme prices, lead forward by a specified number of ticks. The Ichimoku indicator is using the OHLC table. The Ichimoku indicator is composed of two defining lines:


```
Turning Line = ((memoryShort mmax High) + (memoryShort mmin Low)) / 2
```

```
Standard Line = ((memoryMedium mmax High) + (memoryMedium mmin Low)) / 2
```

which are then being used in the spanning line:

```
Leading Span 1 = (prev/)[memoryMedium;] ( Standard Line + Turning Line ) / 2
```

and analogously to the first spanning line, we define the second spanning line

```
Leading Span 2 = (prev/)[memoryMedium;] ((memoryLong mmax High) + (memoryLong mmin Low)) / 2
```

The function updates the table as:

| c	| t	|
|  ----- |  ---- |
| open	| f |		
| high	| f |		
| low	| f |		
| close	| f |		
| closeTurningLine9x26x52	|  f |		
| closeStandardLine9x26x52	|  f |		
| closeLeadingSpan1x9x26x52	|  f |		
| closeLeadingSpan2x9x26x52	|  f |



### MACD

The MACD indicator can be used as:

```
.quantQ.ta.MACD[`price;()!()] tab
 ```
with default values of the ```memoryFast``` being 10, the ```memorySlow``` being 20, and the ```memory``` being 15.

The Moving Average Convergence Divergence, or MACD, is a combination of the moving averages. It is defined as:

```
MACD = (memoryFast mavg price) - (memorySlow mavg price)
```
the signal line is the moving average of the MACD:
```
MACDSignalLine = (memory mavg MACD)
```
and the MACD histogram is then derived from the MACD and the signal line as:
```
MACDhistogram = MACD - MACDSignaLine
```

The function updates the table with price as:

| c	| t	|
| ----- |  ---- |
| price	|  f |
| priceMACDHistogram10x20x15	|  f |
| priceMACDSignalLine10x20x15	|  f |
| priceMACD10x20	|  f |

Remark: The MACD does not have the memory parameter in its name as it does not depend on it.

### Aroon

The Aroon indicator is called as:

```
.quantQ.ta.Aroon[`price`;()!());] tab  
```
with default value of the ```memory``` being 20.

The Aroon indicator is composed of three components: The first two is the number of periods since high/low tick within the past memory ticks:

```
AroonUp = ((memory - lag where price = memory mmax price) / memory)
```
and
```
AroonDown = ((memory - lag where price = memory mmin price) / memory)
```
Lastly, the Aroon oscillator is defined as:
```
AroonOscillator = AroonUp - AroonDown
```
The function updates the table as:

| c	| t	|
| ----- |  ---- |
| price	|  f |
| priceAroonUp20	|  f |
| priceAroonDown20	|  f |
| priceAroonOscillator20 | f |


### Relative Strength Indicator

The Relative Strength indicator is invoked as:

```
 .quantQ.ta.RSI[`open`close;()!()] tabOHLC
 ```
with default values of the ```memory``` being 20. The RSI needs the ```open``` and ```close``` columns to assess the price change per every time step. The two strings pointing to the respective columns need to be specified.

The RSI is defined as:

```
RSI = 100 - (100 / (1 + gainAVG/lossAVG))
```
with
```
gainAVG = memory mavg (close-open) where (close-open)>0
```
```
lossAVG = memory mavg abs[close-open] where (close-open)<0
```

The ratio ```gainAVG/lossAVG``` is between 0 and infinity such that the RSI is between 0 and 100.

The function updates the table as:

| c	| t	|
|  ----- |  ---- |
| open	| f |		
| high	| f |		
| low	| f |		
| close	| f |		
| closeRSI20	|  f |

### Average True Range

The Average True Range indicator is applied as:

```
 .quantQ.ta.ATR[`open`high`low`close;()!()] tabOHLC
 ```
with default values of the ```memory``` being 20. The ATR expects the open/high/low/close table.

The ATR is defined as:
```
ATR = memory mavg TR
```
where True Range is defined as maximum of the three ranges:
```
TR = max[(high-low;abs[high-close];abs[low-close])]
```

The function updates the table as:

| c	| t	|
|  ----- |  ---- |
| open	| f |		
| high	| f |		
| low	| f |		
| close	| f |		
| closeATRI20	|  f |
