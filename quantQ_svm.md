# Support Vector Machine in quantQ

The Support Vector Machine is a popular method used for the binary classification. We have implemented the Soft Margin version of the SVM, which is calibrated using the Stochastic Gradient Descent. In particular, we update the parameters of the model using a single observation at the time, where observations are chosen randomly (without a guarantee that every observation will be chosen).

The model is implemented using an interface with defaults parameters, which allows a quick start. The default settings can be overridden to find the algorithm for the specific problem at hand. We present below a simple tutorial on how to use the Support Vector Machine within ```quantQ``` framework.

The SVM depends on the ```.quantQ.nn``` library with learning functions.

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

We encourage you to try the algorithm by yourselves.
