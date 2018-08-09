
load("GraphClasses.sage")
import itertools
import numpy as np
from scipy.linalg import block_diag

class SimpleGraphLMQC(SimpleGraph):
    def __init__(self,*args,**kwargs):
        """Can be initialized in the same way as SimpleGraph(). In addition, it is possible to specify a partition of the vertex set."""
        super(SimpleGraphLMQC,self).__init__(*args,**kwargs)

        self.partition = kwargs.get('partition')

    def calc_Q_local(self):
        """Uses the partition of the vertex set to calculate which entries of Q_local are possibly non-zero. Returns list of columns."""
        nonzero_pos_A = []
        for i in self.partition:
            for x,y in list(itertools.product(i, i)):
                nonzero_pos_A.append(y+x*self.order())
        self.nonzero_positions = [
                                     [i+0*self.order()**2 for i in nonzero_pos_A]
                                    +[i+1*self.order()**2 for i in nonzero_pos_A]
                                    +[i+2*self.order()**2 for i in nonzero_pos_A]
                                    +[i+3*self.order()**2 for i in nonzero_pos_A]
                                    ][0]

    def from_Qi_to_Q(self):
        """Form Q = (A,B,C,D) from Qi, where Q_i is a list of lists with matrices which are placed on the diagonal of A, B, C or D."""
        Q = []
        for A in self.Q_i:
            M = A[0]
            for m in A[1:]:
                M = block_diag(M,m)
            Q.append(matrix(M))
        self.Q_i = tuple(Q)

    def check_symp_constraint(self,vec_list,debug=False):
        """from list of basis vectors, form Q_i's and check symplectic constraint. Return Q if symplectic"""
        Q_vec = np.sum(vec_list,axis=0).tolist()
        node_list = [len(i) for i in self.partition]
        A = [];B = [];C = [];D = []
        z = len(Q_vec)/4
        i=0
        M = [MatrixSpace(GF(2),j,j) for j in range(1,1+max(node_list))]
        for size_node in node_list:
            A_new = transpose( M[size_node-1](Q_vec[    i:    i+size_node**2]) )
            C_new = transpose( M[size_node-1](Q_vec[2*z+i:2*z+i+size_node**2]) )

            if not A_new.transpose()*C_new==C_new.transpose()*A_new:
                return False,[]

            B_new = transpose( M[size_node-1](Q_vec[  z+i:  z+i+size_node**2]) )
            D_new = transpose( M[size_node-1](Q_vec[3*z+i:3*z+i+size_node**2]) )

            if not B_new.transpose()*D_new==D_new.transpose()*B_new:
                return False,[]

            if not A_new.transpose()*D_new+C_new.transpose()*B_new == matrix.identity(size_node):
                return False,[]

            A.append(A_new);B.append(B_new);C.append(C_new);D.append(D_new)
            i+=size_node**2
        self.Q_i = (A,B,C,D)
        self.from_Qi_to_Q()

        if debug:
            Gamma = self.adjacency_matrix()
            Gamman = self.other.adjacency_matrix()
            A,B,C,D = self.Q_i
            lin_eq_res = Gamman*B*Gamma + D*Gamma + C + Gamman*A
            print "Gamman B Gamma + D Gamma + C + Gamman A = \n", lin_eq_res
            print "A^TC+C^TA = \n", A.transpose()*C+C.transpose()*A
            print "B^TD+D^TA = \n", B.transpose()*D+D.transpose()*B
            print "A^TD+C^TB = \n", A.transpose()*D + C.transpose()*B

        return True,self.Q_i

    def powerset(self,max_length=None):
        "powerset([1,2,3]) -->  (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
        if max_length==None:
            max_length=len(self.iterable)
        return chain.from_iterable(itertools.combinations(self.iterable, r) for r in range(1,max_length+1))

    def is_LMQC_equiv(self,other):
        """Decide whether self and other are LMQC equivalent."""
        assert type(other) in (SimpleGraphLMQC,SimpleGraph,Graph)
        self.other=other

        assert self.vertices()==self.other.vertices(), "The vertex set of self and other are not equal."

        if self.partition is None:
            print "No partition given, assuming every node has one qubit"
            return self.is_LC_eq(self.other)

        assert type(self.partition) is list, "Partition should be a list (of lists)"
        for node in self.partition:
            assert type(node) is list, "Partition should be a list of lists"

        if len(self.partition) == self.order():
            return self.is_LC_eq(self.other)

        self.calc_Q_local()
        iden = matrix.identity(self.order())
        M = np.concatenate((
                        np.kron(iden,other.adjacency_matrix()),
                        np.kron(self.adjacency_matrix(),other.adjacency_matrix()),
                        np.kron(iden,iden),
                        np.kron(self.adjacency_matrix(),iden)
                            ), axis=0)

        M = Matrix(GF(2),M[self.nonzero_positions])
        V = M.kernel()
        V_basis_array = np.asarray(V.basis())

        self.iterable = range(len(V.basis()))
        for combi in self.powerset():
            bool_result,Q_result = self.check_symp_constraint(V_basis_array[list(combi)])
            if bool_result:
                return bool_result,Q_result

        return False,[]

    def run_tests(self):
        assert type(self) in (SimpleGraphLMQC,SimpleGraph,Graph)

        G = SimpleGraphLMQC(Graph({0:[1,2,3]}),**{'partition':[[0,1],[2],[3]]})
        H = SimpleGraphLMQC(Graph({0:[1,3],2:[3]}))
        dum, _ = G.is_LMQC_equiv(H)
        if dum: print "Succesful tested with two LMQC equivalent graphs"
        else: print "Something went wrong"

        G = SimpleGraphLMQC(Graph({0:[1,2,3]}),**{'partition':[[0,1],[2],[3]]})
        H = SimpleGraphLMQC(Graph({0:[3],2:[1]}))
        dum, _ = G.is_LMQC_equiv(H)
        if not dum: print "Succesful tested with two not-LMQC equivalent graphs"
        else: print "Something went wrong"
