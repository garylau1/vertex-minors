# Code for local multi-qubit Clifford equivalence of graph state

## SAGE

The file LMQC.sage contain classes for sage that can be used to study graph states under local multi-qubit Clifford operations.

Using the class SimpleGraph one can efficiently test whether two graphs are equivalent under local complementations and therefore if their corresponding graph states are equivalent under local Clifford operations.

Using the class SimpleGraphLMQC one can test whether two graphs states are equivalent under local multi-qubit Clifford operations (LMQC equivalence), given a partition of the shared vertex set of the two graphs. The partition specifies certain nodes with (multiple) vertices, and in every node with multiple vertices (qubits) multi-qubit Clifford operations are allowed. Local multi-qubit Clifford operations correspond to a sequence of edge addition/deletion and local complementations. One can test LMQC equivalence with G.is_LMQC_eq(Gp). Note that this method is not efficient.

To use the classes, start sage and load the script by writing `load("LMQC.sage")`.

#### Some useful methods:

* `SimpleGraph(data)`: Creates an instance of the class `SimpleGraph`. `data` can be a an instance of the `Graph`-class already in SAGE, a dictionary describing the neighbors of vertices, etc.
* `G1.is_LC_eq(G2)`: Checks if the graph `G1` is LC-equvialent to `G2`, where both graphs are instances of `SimpleGraph`. The method used the algorithm descibed by Bouchet in https://link.springer.com/article/10.1007/BF01275668
* `G.set_partition(T)`: Initializes the partition of `G` as `T`. This is needed for `G.is_LMQC_eq(Gp)`.
* `G.is_LMQC_eq(Gp)`: Checks if `G` is local multi-qubit Clifford equivalent to `Gp` given a partition `T`.

#### Example:
```
G=SimpleGraphLMQC([0,1,2,3])        #Creates a complete graph on the vertices 0,1,2,3
G.set_partition([[0,1],[2],[3]])    #Sets partition of G as [[0,1],[2],[3]]
Gp = SimpleGraph({0:[2,3],1:[]})    #Creates a (disconnected) graph on the vertices 0,1,2,3
G.is_LMQC_eq(Gp)
```
