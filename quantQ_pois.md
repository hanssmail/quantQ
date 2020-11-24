# Poisson distribution and Poisson regression in quantQ

The Poisson process is a random point process where the probability of an event to occur (point) is driven by the intensity function. The process is characterised by the Poisson distribution, which specifies the probability of ```N``` events to occur given the provided arrival intensity.

The extension of the Poisson process is the non-homogenous Poisson process which allows for the intensity to be driven by the exogenous factors. Then, the Poisson regression is a non-linear regression method to estimate the intensity of the non-homogenous Poisson process given the realisation of the process and the exogenous features.

The ```.quantQ.pois``` namespace contains the routines to obtain the Poisson distribution (probability, and the PDF and the CDF in the form of a table), and function to estimate the Poisson regression. Further, the model can be estimated with ```L2``` regularisation and cross-validated.

The Poisson package depends on the ```.quantQ.stats``` and ```.quantQ.opt``` packages to calculate factorial and optimise the functions.

## Examples: Poisson Distribution

First, we explore the Poisson distribution. The library provides three functions. Let us start with calculating the probability to observe ```n``` events given the provided intensity

```
.quantQ.pois.poissonPDF[0.3;3]
```

Further, it is convenient to get the Probability Distribution Function and the Cumulative Distribution Function (or Mass Distribution Function). For that purpose, we can run the following two queries:


```
.quantQ.pois.poissonPDFTab[0.3;()!()]
```

and

```
.quantQ.pois.poissonCDFTab[0.3;()!()]
```

We can proceed with an example of the Poisson regressions.


## Examples: Poisson regression

Let us start with a synthetic data set:

```
nTab:500;
betaTrue: (-1.0;0.5;2.0);
```

where we randomly generate three factors and create intensity as a linear combination of them. We then draw the number of events using the Poisson CDF function:

```
tab: update xVar: (x1,'x2,'x3) from ([] x1:neg[1]+nTab?2.0; x2:nTab?2.0;x3:nTab?3.0);
tab: update yVar: {[x] exec first n from .quantQ.pois.poissonCDFTab[x;()!()] where prob>=first 1?1.0 } each sum betaTrue*flip xVar from tab;
```

The Poisson regression can be performed in a straightforward way as:

```
bucket:()!();
model: .quantQ.pois.regPoisson[tab;bucket];
```
where we have used empty dictionary in order to utilise the default values. The outcome of the model can be accessed as:

```
model
```

and contains the estimated coefficients as well as the maximum likelihood, which is the loss function in the regression.

In order to use the estimated model and predict the values for any dataset, we can proceed as follows:

```
.quantQ.pois.predPoisson[tab;model]
```

The prediction contains both the most likely outcome as well as the sorted list of potential outcomes according to their likelihood (i.e., the 2-nd most likely outcome, the 3-rd most likely outcome...).


### Change the default values

We can override the default values by setting our parameters. In this case, we decrease the precision of the optimisation routine:

```
bucket: enlist[`precision]!enlist[0.01];
```

and re-estimate the model:

```
model: .quantQ.pois.regPoisson[tab;bucket];
model
```

This concludes the basic Poisson regression.


## Regression with regularisation

In addition to the direct Poisson regression, the library contains the regularised version. We have implemented the ```L2``` regularisation of the estimated terms. In order to account for the constant and treat the estimated coefficient of the constant term appropriately, we need to specify in the input if we want to add constant or not (default option is to add the constant).

First, let us run the regularised model with fully pre-set default setup:

```
bucket:()!();
lambda:0.2;
modelReg:.quantQ.pois.regPoissonReg[tab;lambda;bucket];
```

and the estimated model can be accessed with:

```
modelReg
```

In order to apply the model on any dataset and make a prediction, we can use the following routine:

```
.quantQ.pois.predPoissonReg[tab;modelReg]
```

In order to run the regression without adding the constant, we need to override the default behaviour:

```
bucket:enlist[`addConstant]!enlist[0b];
lambda:0.2;
```

and then call the regression

```
modelRegNC:.quantQ.pois.regPoissonReg[tab;lambda;bucket];
```

with the model

```
modelRegNC
```

In the same way, as we did above, we can predict the model without added constant (the information if the constant is added or not being stored in the model):

```
.quantQ.pois.predPoissonReg[tab;modelRegNC]
```

Finally, we can use the ```n```-fold cross-validation technique to find the optimal regularisation parameter (with log-likelihood being optimised)

```
lambdaArray:0.0 0.001 0.01 0.1 0.2 0.3;
bucket:()!();
nCV:6;
([] lambda: lambdaArray; logLNeg: avg each .quantQ.pois.regPoissonCV[tab;bucket;nCV;] each lambdaArray)
```

where the outcome of the cross-validation being organised in a table for better visualisation.

We encourage readers to try the algorithm by yourselves.
