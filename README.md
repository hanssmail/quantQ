# quantQ

The repository for the Machine Learning and Big Data with kdb+/q book by Novotny et al.

Order the book at https://www.wiley.com/en-us/Machine+Learning+and+Big+Data+with+KDB%2B+Q-p-9781119404750.

## Getting Started

```
$ cd quantQ/lib
$ q quantQ.q -p 5000
```

The naming convention for each .q file reflects the corresponding book chapter and the functions it defines reside under a homonymous namespace, as outlined in Section 3.1.0.1.

## Errata

### Chapter 7: Joins

| Section       | Note                
| ------------- |--------------------
| 7.1           | The example of comma join corresponds to ```t1,t6``` and not ```t5,t6```.


| Section       | Note                
| ------------- |--------------------
| 7.2.10        | In the text we state that *we aim to aggregate the data from the table dataSet2 over a window starting 1 minute prior to the trade and ending at the time of the trade.*; however, the ```window``` is defined as starting at the time of the trade and ending 1 minute after the trade. The example should read ```window: (-00:01:00;0) +\: exec time from dataSet1;```.

### Chapter 14: Time Series Econometrics

| Section       | Note                
| ------------- |--------------------
| 14.1.6.1      | Ordering of ```phi``` vector inside implementation of  ```.quantQ.ts.simAR``` should be reversed to be in line with definition 14.1 and example on page 276. Addressed in the repo. Variation of the function using adverbs also provided, under ```.quantQ.tse.simAR```


### Chapter 15: Fourier Transform

| Section       | Note                
| ------------- |--------------------
| 15.3          | Example following implementation of the Hamilton product (15.29-15.32) should read as ```.quantQ.quat.mult[quat1;quat3]```. Definition of ```.quantQ.quat.mult``` has a typo which is fixed in the repo.    


### Chapter 22: Neural Networks

| Section       | Note                
| ------------- |--------------------
| 22.2          | Missing functions ```.quantQ.nn.funcNN``` and ```.quantQ.nn.funcErrNN``` have been added to the repo.

## Extension beyond the book

### Mathematical functions

We have added ```.quantQ.math``` namespace with various mathematical functions, constants and identities. Currently, there are constants, hyperbolic functions, number of special functions and polynomials (defined in the real domain), and the most frequently used PDF and CDF. In order to obtain the multivariate normal distribution, we have included tetrachoric expansion.

### Biostatistics

We have added into the ```.quantQ.stats``` namespace functions to work with contingency tables, namely Exact Fisher test and Barnard test.

### Optimization

We have added a native implementation of Nelder-Mead (amoeba) optimisation method. Functions can be found in the ```.quantQ.amoeba``` namespace. The illustration shows a solution to Rosenbrock's function.  

### Bioinformatics

We have implemented the Needleman-Wunsch algorithm developed in bioinformatics which can be used to align two sequences (nucleotide sequences or general finite-length sequences) using principles of dynamic programming. Functions are within the ```.quantQ.stats``` namespace. The local matching using the Smith-Waterman algorithm is available as well. Levenshtein distance calculated with Wagner-Fischer algorithm has been added.

### Dynamic Time Warp

We have added the algorithms to perform the Dynamic Time Warp calculations in the ```.quantQ.dtw``` namespace. Details of the algorithm can be found [here](data/dtw_intro.pdf).

### Deep Neural Networks

We have added deep neural networks. The features implemented include the dropout, batches, and different functional forms of annealing in learning and regularisation. The architecture is specified in the input dictionary; the rest is fully automated. More examples can be found [here](quantQ_nn_deep.md). The implementation is part of the ```.quantQ.nn``` namespace.

### Stochastic Optimisation

We have implemented the stochastic optimisation, which can be used to minimise the provided n-dimensional function using random search. The provided method is an iterative procedure which explores at every step a neighbourhood comprising of random points lying on the n-sphere, where the radius of the sphere is shrinking with occasional attempts to investigate further points. The implementation is part of the ```.quantQ.so```  namespace.

### Support Vector Machine

We have added the Support Vector Machine for binary classification using the Soft Margin and Stochastic Gradient Descent (using one observation per step). Functions are specified with a set of the default setup, which is customised based on the provided dataset. The Soft Margin specification is provided with regularisation parameter, which can be optimised using the built-in n-fold cross-validation method.

The implementation includes a function which calculates more than 20 statistics used to evaluate the binary classification, including specificity, accuracy, precision, or F1. More examples can be found [here](quantQ_svm.md). The implementation is part of the ```.quantQ.svm``` namespace.

### Poisson Regression

We have added the library with the Poisson distribution and the Poisson regression to estimate the integer-valued Poisson process. The functions are based on the maximum likelihood optimised using routines from ```.quantQ.opt``` library. The library also contains the ```L2```-regularised version with ```n```-fold cross-validation.

Examples and basic usage of the library can be found in [here](quantQ_pois.md). The implementation is part of the ```.quantQ.pois``` namespace.


### Quantum Computing

We have added name space ```.quantQ.quantum``` which contains a set of routines to set up and perform quantum computing using quibits. The library is not connected to any actual quantum computer and is for demonstration purposes only.

Examples and basic library usage can be found in [here](quantQ_quantum.md).

### Topological Data Analysis

We have added the library ```.quantQ.tda``` which allows us to perform the topological data analysis for a cloud on ```n```-dimensional points. The routines provided include calculation of the distance between points, identification of the all Vietoris-Rips complexes given the provided threshold, the routine to calculate all unique loops, and several analytical functions to get insight into data using the TDA. 

Examples and basic library usage can be found in [here](quantQ_tda.md).

### Utility Library

We have added the library ```.quantQ.util``` with various utility functions.
