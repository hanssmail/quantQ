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














































## Examples

We start with generating simple data set, which comprises of two-dimensional independent variable living in ```(0;1)x(0;1)``` support and single binary dependent variable, defined by cutting the support by a linear fragment:

```
data:update y:1+2*{-1*x=0} (x1+x2)>1 from ([] x1:100?1.0; x2:100?1.0);
data: select yVar:y,xVar:(x1,'x2) from data;
```

Let us run the SVM straight ahead, using the default setup of the algorithm. First, we define an empty dictionary such that the algorithm can apply the default choices:

```
bucket:()!();
```

and estimate the model:

```
bucketFit:.quantQ.svm.fitSoftMargin[data;bucket];
```

In order to apply the model on the data set, we need to run:

```
dataModel:.quantQ.svm.evalModel[bucketFit;data];
```

The function applies the estimated model on any data set with ```xVar``` of the same dimensionality. We can use this function to assess the estimate in-sample as well as out-of-sample.

In order to estimate the quality of the fit, we have provided a function which gives a number of common metrics used for evaluating binary fits, which can be accessed as:


```
.quantQ.svm.evalStats[dataModel]
```

The function is general enough and can be used for any table where we need to compare a table with binary variables (-1/1) ```yVar``` and ```yVarPred```.

### Cross-validation

Further, the Soft Margin specification introduces a regularisation parameter, which penalises the ```L2``` norm of the loading coefficients. The value of the regularisation parameter can be found using the cross-validation technique.

First, we set the range of regularisation parameter to be explored and the number of folds to be used:

```
lambdaArray:0.0001 0.01 0.1 1.0;
nCV:6;
bucket:()!();
```

and then we can cross-validated the range of provided values and find the corresponding cross-validated soft margin:

```
{([] lambda:x; softMarginCV:y)}[lambdaArray;] avg each .quantQ.svm.cvFit[data;bucket;nCV;] each lambdaArray
```

### Using Defaults

The default setup of the SVM algorithm can be overridden, and the performance of the algorithm be fine-tuned. Let us list the parameters, which can be changed. We start with an empty dictionary with parameters:

```
bucket:()!();
```

First, we change the function, which governs the learning rate. The default setup is a constant learning rate set to ```0``` throughout the simulation. We set the power-law function defined in the ```.quantQ.nn``` namespace when we have been introducing the deep neural networks:

```
bucket:bucket,(`learningRateFunc`learningRateParams)!(
        .quantQ.nn.learningFuncs;(`method`learningRate`power!(`powerLaw;0.01;-0.5)));
```

The soft margin formulation of the SVM allows us to introduce regularisation into the formula. We can set the regularisation parameter in the ```L2``` norm:

```
bucket:bucket,enlist[`lambda]!enlist[0.2];     
```

We can override the maximum number of iterations the algorithm can perform as well:

```
bucket:bucket,enlist[`maxIters]!enlist[30000];  
```

The algorithm contains a second stopping rule, which is based on the convergence of soft margin. In particular, we use two exponentially moving averages of the soft margin to smooth out the variance in the soft margin. The memory of the moving average---slow and fast---are derived from the sample size, which can be configured as well. Altogether, we can fine-tune this stopping rule as:

```
bucket:bucket,(`sMAlpSlow`sMAlpFast)!(0.9;0.1);
```
and then specify the threshold for the convergence between two moving averages (which is the parameter to be fine-tuned rather than the memory):

```
bucket:bucket,enlist[`softMarginDecayError]!enlist[-0wf];
```

By default, the algorithm returns the solution; however, in many cases, it can be desirable to inspect the convergence path, which can help us to decide about the stopping rule. The default behaviour can be overridden as:

```
bucket:bucket,enlist[`monitorConvergence]!(enlist[1]);
```

Finally, we can run the estimation with the customised dictionary:

```
bucketFit:.quantQ.svm.fitSoftMargin[data;bucket];
```

and visualise the convergence path---the soft margin as a function of the iteration step---where due to a large number of observations, we choose only every 20-th data point:

```
select i, softMargin from bucketFit where (i%20)=floor (i%20), i>0
```
