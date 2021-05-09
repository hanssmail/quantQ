# Technical Analysis

The library ```.quantQ.ta``` comprises the implementation of many standard technical indicators used in the Technical Analysis. The provided set of functions aims to help users with feature engineering and to analyse information in the data. For the sake of simplicity, it is the user's responsibility to have data sorted by the correct time index. Time index is used only when it is necessary, e.g., the time distance between observations is needed.

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

The function adds one or more columns, depending on the indicator. Whenever possible, the parameters are part of the column name, .i.e., the simple moving average technical indicator with 10 step memory applied on a column named ```ret``` update the table with a new column called ```retMA10```. This allows users to add a number of the same technical indicators with different parameters, and the process is down the line (either by visualising the new columns, or, for instance, using a more quantitative approach with the Lasso regression, as implemented in ```.quantQ.reg.LAR```).

## Implemented Functions and Examples

We present an essential use of the provided functions.

### Synthetic Data Generation

We use several tables with synthetically generated data to illustrate functionality.

The simple table with a single dummy column called return

```
tab:([] ret:100?1.0);
```

### Simple Moving Average

The simple moving average is provided as

```
.quantQ.ta.ma[`ret;enlist[`memory]!enlist[20]] tab
```
with default values of the ```memory``` being set to ```10```.

The table is updated as:

| c	| t	|
|  ----- |  ---- |
| ret | f |
| retMA20	| f |

### Exponential Moving Average

The exponential moving average is provided as

```
.quantQ.ta.ema[`ret;enlist[`memory]!enlist[20]] tab
```
with default values of the ```memory``` being set to ```10```. The exponential moving average is updated as:

```
emaNew = column * (2.0%(memory+1))+ emaOld * (1-2.0%(memory+1))
```

The function updates the table as:

| c	| t	|
|  ----- |  ---- |
| ret | f |
| retEMA20	| f |
