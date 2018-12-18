# Code for local multi-qubit Clifford equivalence of graph state

## SAGE

The file LMQC.sage contain classes for sage that can be used to study graph states under local multi-qubit Clifford operations.

Using the class SimpleGraph one can efficiently test whether two graphs are equivalent under local complementations and therefore if their corresponding graph states are equivalent under local Clifford operations.

Using the class SimpleGraphLMQC one can test whether two graphs states are equivalent under local multi-qubit Clifford operations (LMQC equivalence), given a partition of the shared vertex set of the two graphs. The partition specifies certain nodes with (multiple) vertices, and in every node with multiple vertices (qubits) multi-qubit Clifford operations are allowed. Local multi-qubit Clifford operations correspond to a sequence of edge addition/deletion and local complementations. One can test LMQC equivalence with G.is_LMQC_eq(Gp). Note that this method is not efficient.

To use the classes, start sage and load the script by writing `load("LMQC.sage")`.

#### Some useful methods:

* `SimpleGraphLMQC(data)`: Creates an instance of the class `SimpleGraphLMQC`. `data` can be a an instance of the `Graph`-class already in SAGE, a dictionary describing the neighbors of vertices, etc.
* `G1.is_LC_eq(G2)`: Checks if the graph `G1` is LC-equvialent to `G2`, where both graphs are instances of `SimpleGraph`. The method used the algorithm descibed by Bouchet in https://link.springer.com/article/10.1007/BF01275668
* `G.set_partition(T)`: Initializes the partition of `G` as `T`. This is needed for `G.is_LMQC_eq(Gp)`.
* `G.is_LMQC_eq(Gp,method='brute')`: Checks if `G` is local 'T'-multi-qubit Clifford equivalent to `Gp`. Note that there are two methods with an improved runtime, 'gate_tele' and 'conj'. 


#### Example (method=`brute`):
```
G=SimpleGraphLMQC([0,1,2,3])        #Creates a complete graph on the vertices 0,1,2,3
G.set_partition([[0,1],[2],[3]])    #Sets partition of G as [[0,1],[2],[3]]
Gp = SimpleGraph({0:[2,3],1:[]})    #Creates a (disconnected) graph on the vertices 0,1,2,3
G.is_LMQC_eq(Gp)
```
#### Example (using all three methods)
```
G=SimpleGraphLMQC({0:[2,3],1:[2]})        #Creates a graph on the vertices 0,1,2,3
G.set_partition([[0,1],[2],[3]])    #Sets partition of G as [[0,1],[2],[3]]
Gp = SimpleGraph({2:[0,1,3],1:[3]})    #Creates a  graph on the vertices 0,1,2,3
print G.is_LMQC_eq(Gp,method='brute')
print G.is_LMQC_eq(Gp,method='conj')
print G.is_LMQC_eq(Gp,method='gate_tele')
```
#### Example (two random graph states)
```
G = SimpleGraphLMQC(graphs.RandomGNP(4,0.7))
Gp = SimpleGraph(graphs.RandomGNP(4,0.9))
G.set_partition([[0],[1],[2,3]])
print G.is_LMQC_eq(Gp,method='brute')
print G.is_LMQC_eq(Gp,method='conj')
print G.is_LMQC_eq(Gp,method='gate_tele')
```
#### Example (one three-qubit node) (!THIS COULD TAKE LONG!)
```
G = SimpleGraphLMQC(Graph({0:[3],1:[4],2:[5]}))
Gp = SimpleGraph(graphs.RandomGNP(6,0.9))
G.set_partition([[0,1,2],[3],[4],[5]])
G.is_LMQC_eq(Gp,method='brute')
```
