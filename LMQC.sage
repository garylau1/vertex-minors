
load("GraphClasses.sage")
import itertools
import numpy as np
from scipy.linalg import block_diag
from itertools import combinations,product

class SimpleGraphLMQC(SimpleGraph):
    def __init__(self,*args,**kwargs):
        """Can be initialized in the same way as SimpleGraph(). In addition, it is possible to specify a partition of the vertex set."""
        super(SimpleGraphLMQC,self).__init__(*args,**kwargs)
        # self.partition = kwargs.get('partition')

    def set_partition(self,partition):
        self.partition = partition

    def calc_Q_local(self):
        """Uses the partition of the vertex set to calculate which entries of Q_local are possibly non-zero. Returns list of columns."""
        #We start by only looking at A of Q. Then we know that A is a block diagonal matrix. Nonzero_pos_A is a list of numbers which indicate the position of possible nonzero elements,
        #where every number corresponds to an element of list(Nonzero_pos_A). So (0,1) is the second element, and (1,0) is the self.order()+1 th element.
        nonzero_pos_A = []
        for i in self.partition:
            for x,y in list(product(i, i)):
                nonzero_pos_A.append(y+x*self.order())
        #Now we shift the positions for A to the positions for B,C,D. They are {1,2,3}*self.order()**2 further in the list respectively.
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
        return chain.from_iterable(combinations(self.iterable, r) for r in range(1,max_length+1))

    def is_LMQC_eq(self,other,method = 'brute'):
        """Decide whether self and other are LMQC equivalent."""
        assert type(other) in (SimpleGraphLMQC,SimpleGraph,Graph)
        self.other=other
        assert self.vertices()==self.other.vertices(), "The vertex set of self and other are not equal."
        if self.partition is None:
            print "No partition given, assuming every node has one qubit"
            return self.is_LC_eq(self.other,allow_disc=True)
        assert type(self.partition) is list, "Partition should be a list (of lists)"
        for node in self.partition:
            assert type(node) is list, "Partition should be a list of lists"
        if self.is_LC_eq(self.other,allow_disc=True):
            return True
        if len(self.partition) == self.order():
            return self.is_LC_eq(self.other,allow_disc=True)
        multi_qubit_nodes = [i for i in self.partition if len(i)>1]
        if method=='gate_tele':
            for a,b in multi_qubit_nodes:
                try: self.delete_edge(a,b)
                except: pass
            Gp_list = [self]
            Gtmp = Graph()
            Gtmp.add_vertices(['A0','A1','Aa','Ab'])
            F_list = [copy(Gtmp) for _ in range(10)]
            edges_list = [[('A0', 'Aa'), ('A1', 'Ab')],[('A0', 'Aa'), ('A1', 'Ab'), ('Aa', 'Ab')],[('A0', 'Aa'), ('A1', 'Aa'), ('A1', 'Ab')],[('A0', 'Ab' ), ('A1', 'Aa')],[('A0', 'Ab'), ('A1', 'Aa'), ('Aa', 'Ab')],[('A0', 'Aa' ), ('A0', 'Ab'), ('A1', 'Aa')],[('A0', 'Ab'), ('A1', 'Aa'), ('A1', 'Ab')],[('A0', 'Aa' ), ('A0', 'Ab'), ('A1', 'Ab')],[('A0', 'A1' ), ('A0', 'Aa'), ('A1', 'Ab')],[('A0', 'A1' ), ('A0', 'Ab'), ('A1', 'Aa')]]
            for i in range(10): F_list[i].add_edges(edges_list[i])
            for a,b in multi_qubit_nodes:
                tmp_list = []
                meas = ([a,b,'A0','A1'],['Z','Z','Z','Z'])
                for Gp in Gp_list:
                    if not set(Gp.neighbors(a)).union(Gp.neighbors(b)) - set([a,b]):
                        #if a,b are disconnected from the other vertices in Gp
                        Gp.flip_edge((a,b))
                        tmp_list.append(Gp)
                    else:
                        for F in F_list:
                            G_F = SimpleGraph(G.union(F))
                            G_F.add_edges([(a,'A0'),(b,'A1')])
                            G_F.tau_seq([a,'A0',a,b,'A1',b],inplace=True)
                            for corr_a,corr_b in [("I","I"),("S","I"),("I","S"),("S","S")]:
                                G_oper = copy(G_F)
                                if corr_a == "S":
                                    G_oper.tau('A0',inplace=True)
                                if corr_b == "S":
                                    G_oper.tau('A1',inplace=True)
                                G_oper.meas_seq(meas[0],meas[1],inplace=True)
                                G_oper.relabel({"Aa":a,"Ab":b})
                                tmp_list.append(G_oper)
                Gp_list = tmp_list
                tmp_list = []
            eq = False
            for G_op in Gp_list:
                if G_op.is_LC_eq(self.other,allow_disc=True):
                    eq = True
                    break
            return eq
        elif method == 'conj':
            assert self.is_connected() and other.is_connected(), "Self or other is not a locally connected graph"
            for a,b in multi_qubit_nodes:
                try: self.delete_edge(a,b)
                except: pass
            self.calc_Q_local()
            iden = matrix.identity(self.order())
            #Vectorization of the linear equation (1|adj')PQ(1|adj)^T = 0 with Q = (A,B,C,D)
            M = np.concatenate((
                            np.kron(iden,other.adjacency_matrix()),
                            np.kron(self.adjacency_matrix(),other.adjacency_matrix()),
                            np.kron(iden,iden),
                            np.kron(self.adjacency_matrix(),iden)
                                ), axis=0)
            #We can disregard part of this matrix as they correspond to neccessary zero elements of the solution vector.
            M = Matrix(GF(2),M[self.nonzero_positions])
            #and we solve for the vector (A,B,C,D).
            V = M.kernel()
            V_basis_array = np.asarray(V.basis())
            #Now we check every vector in the linear solution space to see if it satisfies the symplectic constraint
            self.iterable = range(len(V.basis()))
            num_multi_qubit_nodes = len([i for i in self.partition if len(i)>1])
            for combi in self.powerset(2+2*num_multi_qubit_nodes):
                bool_result,Q_result = self.check_symp_constraint(V_basis_array[list(combi)])
                if bool_result:
                    # return bool_result,Q_result
                    return True
            elif method == 'brute':
                self.calc_Q_local()
                iden = matrix.identity(self.order())
                #Vectorization of the linear equation (1|adj')PQ(1|adj)^T = 0 with Q = (A,B,C,D)
                M = np.concatenate((
                                np.kron(iden,other.adjacency_matrix()),
                                np.kron(self.adjacency_matrix(),other.adjacency_matrix()),
                                np.kron(iden,iden),
                                np.kron(self.adjacency_matrix(),iden)
                                    ), axis=0)
                #We can disregard part of this matrix as they correspond to neccessary zero elements of the solution vector.
                M = Matrix(GF(2),M[self.nonzero_positions])
                #and we solve for the vector (A,B,C,D).
                V = M.kernel()
                V_basis_array = np.asarray(V.basis())
                #Now we check every vector in the linear solution space to see if it satisfies the symplectic constraint
                self.iterable = range(len(V.basis()))
                num_multi_qubit_nodes = len([i for i in self.partition if len(i)>1])
                for combi in self.powerset():
                    bool_result,Q_result = self.check_symp_constraint(V_basis_array[list(combi)])
                    if bool_result:
                        return True
                        # return bool_result,Q_result
            #If none of the vectors satisfy the symplectic constraint, return False
            # return False,[]
            return False

    def run_tests(self):
        assert type(self) in (SimpleGraphLMQC,SimpleGraph,Graph)

        G = SimpleGraphLMQC(Graph({0:[1,2,3]}),**{'partition':[[0,1],[2],[3]]})
        H = SimpleGraphLMQC(Graph({0:[1,3],2:[3]}))
        dum = G.is_LMQC_eq(H)
        if dum: print "Succesful tested with two LMQC equivalent graphs"
        else: print "Something went wrong"

        G = SimpleGraphLMQC(Graph({0:[1,2,3]}),**{'partition':[[0,1],[2],[3]]})
        H = SimpleGraphLMQC(Graph({0:[3],2:[1]}))
        dum = G.is_LMQC_eq(H)
        if not dum: print "Succesful tested with two not-LMQC equivalent graphs"
        else: print "Something went wrong"
