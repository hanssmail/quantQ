# Topological Data Analysis

The library ```.quantQ.tda``` contains basic building blocks to perform the topological data analysis (TDA). The approach is based on calculating properties of Vietoris-Rips complexes (all Rips complexes, as they are known in the literature). The metric used is ```L_p``` norm, where ```p``` is a choice parameter.

Waring: The TDA involves working with graphs, and many problems in graph theory are computationally very expensive. We have implemented all the algorithms from scratch, and some implementation may not be the most optimal. It is, therefore, necessary to start analyzing the data at a small scale and increase the scope (and threshold for the distance between points to be considered connected). We leave with readers to calculate how many different single loops (with one eye) can create from, say, 20 fully connected points. This is the dimensionality we are working with.  

Last but not least, we encourage readers to either contribute and provide a more efficient implementation.

We introduce the concepts in the TDA implemented in our library using particular toy examples. Besides, we add a few routines to perform birth-death analysis, which is getting traction in the literature.  

## Concepts and Examples

Let us start with a set of ```n```-dimensional points. We simulate a set of random points using:

```
dataSet: update index: i from ([] point:(20 3)#neg[0.5]+60?1.0);  
```

We work with three-dimensional points. In order to define the topology, we need to provide a metric (distance between points). We use the ```L_2``` norm, which we have wrapped as follows:

metric: .quantQ.tda.distance[.quantQ.tda.Lp[2]];

where ```metric``` is a function expecting two arrays of the same size and returns their distance. We can provide any metric (even distance between two words based on alphabetical distance, if we wish to analyze vocabulary).

Having the metric ```metric```, we create a table with all edges within our dataset and assign length to every edge using metric:

```
tabEdges:.quantQ.tda.getAllEdges[dataSet;metric];
```

Our implementation calculates all the combinations among the points in the cloud. We work with all points as we later analyze the topological objects' properties on the maximum distance between points to be considered connected.

The basic building block for further analysis is a Vietoris-Rips complex with parameter ```thr```, which is defined as a set of points from the cloud of points, such that the distance between every pair of points is smaller or equal to ```thr```.

The smallest VR complex is an edge of the length smaller or equal to ```thr```. Such a complex is trivial. The first interesting complex is a triangle, an object composed of three points, connected with an edge having the length smaller or equal to ```thr```.

We can recover all the triangles as:

```
thr:0.4;  
arrayTriangles:.quantQ.tda.getTriangles[tabEdges;thr];    
arrayTriangles
```

We have stored the triangles into an array as we will use it below to derive complexes of larger cardinality.

From triangles, it is natural to ask a question about larger complexes, i.e., a set of points where every mutual pair-wise distance is below the threshold ```thr```.

First, we calculate the most extensive possible complexes. The algorithm works as follows: We start with a first triangle (seed) and search among all points and form the most extensive possible complex. Once the complex is formed, we store the complex into the final output and remove all the triangles that have contributed to the new complex from the array of triangles considered seeds.

```
.quantQ.tda.getAllVR[arrayTriangles]   
```

Possible outcome can look as follows:

```
0 3 5 14
0 3 6
1 4 11 18
1 4 16
2 3 5 14
2 5 8 14
6 12 19
10 15 16
```

We see that among the list of complexes, there are two with an overlapping triangle. To get a (one possible) list of largest constructed complexes, we add another routine, which omits any overlapping complexes. This is achieved in a greedy fashion where once we construct a large complex, we discard all the nodes contributing to the complex in the consequent iterations   

```
.quantQ.tda.getAllVRUnique[arrayTriangles]
```

with a corresponding outcome:

```
0 3 5 14
0 3 6
1 4 11 18
1 4 16
2 3 5
2 3 14
2 5 8 14
6 12 19
10 15 16
```

As the algorithm is greedy, the outcome depends on the order of the triangles going as an input. As opposed to the greedy approach, we can build a utility-based algorithm, resulting in a (more) unique outcome.

The complexes are connected to the generators of the 0-homology group, and their properties are the basis for the TDA analysis. We often study the cardinality of all complexes, their histograms, or their stability as a function of the threshold ```thr```.

One particular class are persistent diagrams. The persistent diagram is a visual way to understand the threshold (scale) at which complexes emerge and how they die. This is called the born-death diagram. These diagrams are stable to noise and can serve as a tool to classify underlying datasets.

The simplest example of the persistent diagram is a case of isolated nodes, i.e., nodes not belonging to any edge given a threshold ```thr```. The nodes (points) are created at ```thr=0``` and with increasing ```thr```, they can become part of an edge. We can thus easily calculate how many nodes are part of at least one edge for a given ```thr```, for which we have prepared a function:

```
.quantQ.tda.persistentGraphPoint[tabEdges;] 0.1+0.02*til 20
```

The outcome is an increasing function. To get the increments (in terms of single dead nodes), we apply ```deltas```.

Further, we can analyse the persistence of cardinality for complexes with various number of nodes. We provide a function, which returns for every threshold a histogram. We can use the function and focus on a complex of certain size:

```
{select thr, numberVR from x where orderVR=3}
.quantQ.tda.nComplexesHistThr[tabEdges;] (0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6)
```

giving

| thr	| numberVR |
| --- | -------- |
| 0.3	| 2 |
| 0.35	| 1 |
| 0.4	| 4 |
| 0.45	| 6 |
| 0.5	| 3 |
| 0.6	| 2 |

We see that number of complexes of size 3 is not ever-increasing with an increasing threshold, as with a large threshold, the larger objects absorb the smaller complexes. We leave on readers to verify that we ultimately end up with one large complex, including all points.

## Loops

Another interesting object within the TDA domain is loops. We define a loop as ordered list nodes, where every two consecutive nodes are connected by an edge (and the last node is connected to the first node -- we include the first node as a last node in the calculations).

The first non-trivial loop is a 3-node complex. This loop is not considered a loop in the literature. The literature considers a pure loop to be a connected path of nodes with a (single) hole. We will get to this point in a few paragraphs. First, let us start with calculating all the loops:

```
arrayEdges: exec edge from tabEdges where distance<=thr;
.quantQ.tda.getAllLoops[arrayEdges]  
```

where we have created a function to take as input all the edges and work with edges only. Let us stress that we denote edge as ```(0 1)``` which is equivalent to ```(1 0)```, i.e., edges are not ordered, while loops are (```(0 5 3 1 0)``` is different to ```(0 3 5 1 0)```).

The function that calculates all the loops is computationally costly, and we encourage users to scale the problem carefully. The computational complexity behind working with graphs is challenging.

Further, we consider a subset of loops, which are defined as:
Being a loop as defined above
There exists at least one pair of points on the loop, which are not connected by an edge
There are no three points on the loop which would create a complex

Condition 2) excludes the loops with multiple eyes, while condition 3) ensures that the loops have minimum width (we can consider threshold ```thr``` as a resolution, below which we are not able to distinguish two points. Thus, a doughnut-like object can be a ring, which we want to exclude).

We can calculate the `pure` loops using following function:

```
thr:0.4;
.quantQ.tda.getAllPureLoops[tabEdges;thr]  
```

where the function returns a list of loops. We have another implementation for the same purpose, ``` .quantQ.tda.getAllPureLoops_NaiveFilter``` which accepts the list of all loops and all edges and then prune the list of loops. This implementation is less efficient relative to the previously mentioned function.

The (pure) loops are linked to the generators of the 1-homology group and are frequently used in the TDA to classify the data. The persistent graphs can be defined for loops as well.

First, we can focus on all the loops (non-pure ones) and count number of different loops as a function of threshold:

```
.quantQ.tda.nLoopsThr[tabEdges;]  0.2+0.005*til 160
```

The more important information about the data set is encoded in the pure loops. There is also less number of pure loops relative to all possible loops. Thus, we provide a function, which accepts the table with all edges and list of threshold to analyse and it returns for every threshold information about the threshold at which the loop was born and also the value if its death:

```
thrList: (0.25 0.26 0.27 0.28 0.29 0.30 0.31 0.32 0.33 0.34 0.35 0.36 0.37 0.38 0.39 0.40 0.41 0.42 0.43 0.44 0.45 0.46 0.47 0.48 0.49 0.50 0.51 0.52);
loopAnalytics:.quantQ.tda.bornAndDeathLoops[tabEdges;thrList]
```

The table produced by the function looks as follows:

| loop	| born	| death	| sizeLoop	| isBornBefore	| isDeathAfter |
|  ----- |  ---- | ----- | ----- | ----- | ----- |
| 0 3 2 5 0	| 0.25	| 0.28	| 4	| true	| false |
| 0 3 2 5 8 14 0	| 0.25	| 0.26	| 6	| true	| false |
| 0 3 2 5 14 0	| 0.25	| 0.26	| 5	| true	| false |
| 0 3 2 8 5 0	| 0.25	| 0.28	| 5	| true	| false |

We can further work with this table and get various marginal persistent diagrams: A simple persistent scatter graph can be obtained with:

```
select born, death from loopAnalytics
```

or the born process of loops of length 4 can be analyzed as:

```
update sums nBorned from
select nBorned: count sizeLoop by born from loopAnalytics where sizeLoop=4
```

There are much more possibilities to employ the provided function and obtain new insight into your data set. This is a tester to illustrate the library usage with a synthetic dataset. We encourage readers to try the algorithm by yourself using your real datasets.
