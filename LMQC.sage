load("Github/vertex-minors/LC.sage")
load("Github/vertex-minors/GraphClasses.sage")
import itertools
import numpy as np
from scipy.linalg import block_diag

class SimpleGraphLMQC(SimpleGraph):
    def __init__(self,*args,**kwargs):
        """Can be initialized in the same way as Graph(). An additional option is to give a list, which by default will give a complete graph on the vertices provided in the list, if format='empty' is given the graph is instead a graph on the vertices in the list with no edges."""

        try:
            super(SimpleGraph,self).__init__(*args,**kwargs)
            try:
                self._order=list(args[0].get_order())
            except:
                self._order=self.vertices()
        except StandardError as e:
            try:
                data=list(args[0])
                n=len(data)
                if kwargs.has_key('format'):
                    if kwargs['format']=='empty':
                        A=zero_matrix(GF(2),n)
                    else:
                        A=ones_matrix(GF(2),n)+identity_matrix(GF(2),n)
                else:
                    A=ones_matrix(GF(2),n)+identity_matrix(GF(2),n)
                super(SimpleGraph,self).__init__(A)
                self.relabel(data)
                self._order=self.vertices()
            except:
                raise e
        #CM: This is to add partition T of V to self
        try:
            self.partition = kwargs["T"]
        except:
            self.partition = [[i] for i in range(self.order())]

    def calc_Q_local(self):
        """Uses the partition of the vertex set to calculate which entries of Q_local are possibly non-zero. Returns list of columns."""
        nonzero_pos_A = []
        for i in self.partition:
            for j in list(itertools.product(i, i)):
                x,y = j
                nonzero_pos_A.append(y+x*self.order())
        self.nonzero_positions = [nonzero_pos_A
                                +[i+self.order()**2 for i in nonzero_pos_A]
                                +[i+2*self.order()**2 for i in nonzero_pos_A]
                                +[i+3*self.order()**2 for i in nonzero_pos_A]][0]

    def from_Qi_to_Q(self,Q_i):
        Q = []
        for A in Q_i:
            M = A[0]
            for m in A[1:]:
                M = block_diag(M,m)
            Q.append(matrix(M))
        return tuple(Q)

    def check_symp_constraint(self,vec_list,test=False):
        Q_vec = np.sum(vec_list,axis=0).tolist()
        node_list = [len(i) for i in self.partition]
        A = [];B = [];C = [];D = []
        z = len(Q_vec)/4
        i=0
        M = [MatrixSpace(GF(2),j,j) for j in range(1,1+max(node_list))]
        for size_node in node_list:
            A.append( transpose( M[size_node-1](Q_vec[    i:    i+size_node**2]) ) )
            B.append( transpose( M[size_node-1](Q_vec[  z+i:  z+i+size_node**2]) ) )
            C.append( transpose( M[size_node-1](Q_vec[2*z+i:2*z+i+size_node**2]) ) )
            D.append( transpose( M[size_node-1](Q_vec[3*z+i:3*z+i+size_node**2]) ) )

            i+=size_node**2
            if test:
                continue
            if transpose(A[-1])*D[-1]+transpose(C[-1])*B[-1] == matrix.identity(size_node):
                pass
            else:
                return False,[]

            if transpose(A[-1])*C[-1]==transpose(C[-1])*A[-1]:
                pass
            else:
                return False,[]

            if transpose(B[-1])*D[-1]==transpose(D[-1])*B[-1]:
                pass
            else:
                return False,[]
        return True,self.from_Qi_to_Q((A,B,C,D))

    def powerset(iterable,max_length=None):
        "powerset([1,2,3]) -->  (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
        s = list(iterable)
        if max_length==None:
            max_length=len(s)
        return chain.from_iterable(itertools.combinations(s, r) for r in range(1,max_length+1))

    def is_LMQC_equiv(self,other):

        self.calc_Q_local()

        if len(self.nonzero_positions) == 4*self.order():
            return self.is_LC_eq(other)

        n=self.order()
        M = np.concatenate(
        (np.kron(matrix.identity(n),other.adjacency_matrix()),
                        np.kron(self.adjacency_matrix(),other.adjacency_matrix()),
                        np.kron(matrix.identity(n),matrix.identity(n)),
                        np.kron(self.adjacency_matrix(),matrix.identity(n))),
        axis=0)
        M = Matrix(GF(2),M[self.nonzero_positions])

        V = M.kernel()
        V_basis_array = np.asarray(V.basis())
        for combi in powerset(range(len(V.basis()))):
            bool_result,Q_result = self.check_symp_constraint(V_basis_array[list(combi)])
            if bool_result:
                return bool_result,Q_result

        return False,[]
