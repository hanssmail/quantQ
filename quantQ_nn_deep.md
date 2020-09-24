# Deep Neural Networks in quantQ

The neural network namespace contains the set of functions, which introduce the deep neural networks into q. The single layer neural network can be implemented by the deep framework and thus this supersedes the existing neural network function.

The number of features implemented in the deep neural network is inherited from the single layer functions (classifier, multivariate classifier, and regressions, for example) and thus we recommend users to go through this section in the book.

The deep neural network estimation can
a be invoked by

```
.quantQ.nn.modelNNDeep[input;output;bucket]
```
where ```input``` is the array of the inout data, ```output``` is the array of the outputs we want to model, and ```bucket``` is a dictionary with parameters (see more below).

The function sets the architecture and iterates over the forward-backward algorithm implemented in ```.quantQ.nn.forwardBackwardNNDeepwErrorDrop```.

Once we will estimate the model, we can obtain predicted values (in or out-of-sample) using ```.quantQ.nn.predictNNDeep``` function, which takes as an input the features and the estimated model.

From user perspective, it is worth to mention ```.quantQ.nn.learningFuncs```, which defines number of common ways how we can parametrise the time learning/regularisation throughout the iteration life-cycle: constant, power-law decay of the parameters, step-wise decay, and power-law annealing (the parameter decays, while it regularly jumps up, which in the case of learning parameter allows escaping from the local optima).

## Examples

We illustrate the deep neural network using the iris data set discussed in the book.

First, load the data:

```
irisData:("FFFFS"; enlist ",") 0:`$"irisData.csv"
```

Then, we split the data set into in-sample and out-of-sample parts:

```
randomDraw100:neg[100]?count irisData;
irisDataIN: select from irisData where i in randomDraw100;
irisDataOUT: select from irisData where not i in randomDraw100;
```

In addition, the problem is multivariate classification and thus we need to define the distinct classes we want to estimate. This can be achieved by

```
distinctClass:distinct irisDataIN`class;
```
where we will use the defined classes in the function ```.quantQ.nn.encoderMulti```, which specifies the output in the suitable form.

The input and output arrays are given as:

```
input: flip (irisDataIN`sepalLength;irisDataIN`sepalWidth;irisDataIN`petalLength;irisDataIN`petalWidth);
output: .quantQ.nn.encoderMulti[distinctClass] irisDataIN`class;
```

Further, we specify three functions: First, the learning rate used to update the coefficients after every iteration. We will use step-wise decaying learning rate starting at value ```0.01``` and decreasing by factor ```0.8``` every ```30``` steps:

```
learningRateFunc:.quantQ.nn.learningFuncs[(`learningRate`method`frac`step)!(0.01;`stepwise;0.8;30);];
```

The two remaining functions characterise the dropout rate at the input layer, and in the hidden layers. For both, we use constant functions with probability to keep the input feature equal to ```0.9```, and the hidden node equal to ```0.8```:

```
pKeepInputFunc:.quantQ.nn.learningFuncs[(`learningRate`method)!(0.9;`const);];
pKeepHiddenFunc:.quantQ.nn.learningFuncs[(`learningRate`method)!(0.8;`const);];
```

Finally, we specify the dictionary with all the parameters:

bucket:(`nNeurons`typeTerminal`argTerminal`learningRateFunc`pKeepInputFunc`pKeepHiddenFunc`monitorConvergence`func)!(
(3 3) ;`count ;1000 ;learningRateFunc ;pKeepInputFunc;pKeepHiddenFunc;1;`multiClassifier);

where, among others, the parameter ```nNeurons``` is an array which specifies the number of hidden nodes layer by layer, starting at the input. Changing this array thus opens the user the spectrum of different architectures and user can choose the one which fits the problem at hand the best (it is not surprising that iris data set can be described well even with single layer).

The estimation itself is run as:

```
model:.quantQ.nn.modelNNDeep[input;output;bucket]
```
The error of the estimators during the iterations can be accessed as:

```
select Counter, Error from model
```

and the predicted values are available with (this function have been missing in the book)

```
.quantQ.nn.predictNNDeep[input;last model]
```

### Using defaults

The implementation allows us to rely on the default values. The minimum set of parameters for the iris data set reads:

```
bucket:(`nNeurons`typeTerminal`argTerminal`func)!(
(6 6);`count;100000;`multiClassifier);
```

and estimation is then run as usual:

```
model:.quantQ.nn.modelNNDeep[input;output;bucket]
```

### Using Batches

The implementation allows us to train network in batches. The default value is one batch, i.e., iterating through the data set at once.

We can easily specify the deep neural network, where data set is split into three batches and iteration through the entire data set is done in three steps. Further, we want to run ```100``` epochs:

```
bucket:(`nNeurons`typeTerminal`argTerminal`learningRateFunc`pKeepInputFunc`pKeepHiddenFunc`monitorConvergence`func`nEpochs`nBatches
  )!(
(enlist 6);`count;1000;learningRateFunc;pKeepInputFunc;pKeepHiddenFunc;1;`multiClassifier;100;3);
```

The estimation is run as usual:

````
model:.quantQ.nn.modelNNDeep[input;output;bucket]
```
