# Quantum Computation in quantQ

The library ```.quantQ.quantum``` is derived from the qiskit implementation of the quantum computing. This library aims to show the versatility of q and learn elements of quantum computing. For readers seriously interested in quantum computing, we recommend qiskit as a nice platform to start with as well as to learn about quantum computing, see  
[here](https://qiskit.org/).

In a nutshell, we define two objects: the quantum state composed of ```n``` quibits, and the bench, which is composed of a set of quantum gates. The quantum states propagate through the bench---it interacts with the individual elements---and changes its state. At the end, we take a measurement of the propagated state.

We have implemented most of the standard single quibit gates, see ```key .quantQ.quantum.singleDict```:

```
`PauliX`PauliY`PauliZ`One`Null`U2`U1`H`S`Sconj`T`Tconj`Rx`Ry`Rz
```

and one double quibit gate, namely controlled not, or CNOT.

## Example

Let us illustrate the library on a simple example: The system is composed of three quibits (the system is composed of three entangled spin-1/2 particles)

```
Nquibits:3;
```

We create the initial state:

```
state: .quantQ.quantum.setSystem[Nquibits];
```

A quibit represents each individual particle: ``` \alpha*|0>+\beta*|1>```, where ```\alpha``` and ```\beta``` are complex numbers (we utilise the ```.quantQ.complex``` library to work with complex numbers). The initial state is a tensor product of three individual states such that the three-particle state is a linear composition of ```|0 0 0>```, ```|0 0 1>```, ```|0 1 0>``` etc. The dimension of the three-particle state is ```2^3```.

It is worth to clarify the indexing of quibits. In particular, the indices are sorted as follows ```|quibit_2 quibit_1 quibit_0>```. The indices are running from right to left. If we want to create our gate, we need to consider this.  

The initial state is composed of a state, where each particle is in pure ```|0>```. Thus, the only non-zero coefficient is with ```|0 0 0>``` element.

Second, we intialize the bench:

```
bench0: .quantQ.quantum.setBench[Nquibits];
```

A bench is an object with two elements: The first one---```bench```---is a list of matrices, which corresponds to elements of the gate. Each matrix is ```2^{Nquibits}x2^{Nquibits}``` unitary matrix (preserving probability such that particles do not escape nor are created within our system). The second object---```view```---is a table, which helps us to visualise how the bench looks, what elements have been added and in what order.


The bench itself can be accessed as:

```
bench0[`bench]
```

while the table with visualised output can be obtained as:

```
bench0[`view]
```

The initialised bench corresponds to the identity matrix, which does not change the quantum states.

Further, we add two operations into the bench. First, we apply Hadamard gate on 0-th quibit:

```
bench1:  .quantQ.quantum.appendSingleQuibitMat[bench0;`H;0;()!()];
```

then, we add the controlled not gate, where 0-th quibit is a control quibit, while 1-st quibit is the target one:

```
bench2: .quantQ.quantum.appendCNOTMat[bench1;0;1];
```

Every time we add a gate, we add a unitary matrix into the list of matrices. It is possible to add more than one single quibit gate simultaneously, but we have not implemented the formalism in that way (the outcome would be equivalent).

We can visualise the bench as:

```
bench2[`view]
```

The outcome looks like follows (using qStudio and its unappealling representation of the list of complex numbers):

| index	| quibit0	| quibit1	| quibit2|
| ------- | -------- | -------- | ------- |
| 0	| (kx.c$Flip@ba9ea38;kx.c$Flip@1b4591d4)	| (kx.c$Flip@3f1a2e4b;kx.c$Flip@198f00c5)	| (kx.c$Flip@5d33d1e6;kx.c$Flip@6a4f20a5) |
| 1	| H	|  	|  |
| 2	| cCNOT	| tCNOT	|  |

Besides the intuitive single quibit gate, the symbols distinguish the control and target quibits of the CNOT gate, see c/t.

Once we have set the bench, we take the initial state and propagate it through the bench as follows:

```
state2: .quantQ.quantum.evolveState[state;bench2];
```

The state itself is ```2^3``` vector of complex numbers. What is important are the squares of the coefficients, which correspond to probabilities to find a system in a given configuration. We have for this purpose a function. In addition, we have arranged the outcome into a nice table, where we can see initial and final states next to each other:

```
(1!`state`probs0 xcol .quantQ.quantum.measurement[state])
uj
(1!`state`probs2 xcol .quantQ.quantum.measurement[state2])
```

The result reads:

| state	| probs0 | probs2 |
| ------- | ------ | ------ |
| 0 0 0	| 1.0	| 0.5 |
| 0 0 1	| 0.0	| 0.0 |
| 0 1 0	| 0.0	| 0.0 |
| 0 1 1	| 0.0	| 0.5 |
| 1 0 0	| 0.0	| 0.0 |
| 1 0 1	| 0.0	| 0.0 |
| 1 1 0	| 0.0	| 0.0 |
| 1 1 1	| 0.0	| 0.0 |

This is the expected outcome as the quibits of 0-th and 1-st particle are aligned, while the quibit of the 2-nd particle (the one with the first index in the state column) is 0 as it was initialised in this state, and there was no gate affecting 2-nd particle.

We encourage readers to try the algorithm by yourselves.
