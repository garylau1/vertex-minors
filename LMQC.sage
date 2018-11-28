
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

    def is_LMQC_eq(self,other):
        """Decide whether self and other are LMQC equivalent."""
        assert type(other) in (SimpleGraphLMQC,SimpleGraph,Graph)
        self.other=other

        assert self.vertices()==self.other.vertices(), "The vertex set of self and other are not equal."

        if self.partition is None:
            print "No partition given, assuming every node has one qubit"
            return self.is_LC_eq(self.other,allow_disc=True),[]

        assert type(self.partition) is list, "Partition should be a list (of lists)"
        for node in self.partition:
            assert type(node) is list, "Partition should be a list of lists"

        if len(self.partition) == self.order():
            return self.is_LC_eq(self.other,allow_disc=True),[]

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
        for combi in self.powerset():
            bool_result,Q_result = self.check_symp_constraint(V_basis_array[list(combi)])
            if bool_result:
                return bool_result,Q_result

        #If none of the vectors satisfy the symplectic constraint, return False
        return False,[]

    @staticmethod
    def generateTwoEquivRandomGraphs(V,nlist = []):
        #Generate node list
        # nlist =[]
        two_qubit_oper_list = []
        two_qubit_qubits = []
        i = 0
        if nlist == []:
            while i < V:
                coin = randint(0,3)
                if coin != 1:
                    nlist.append([i])
                    i+=1
                else:
                    nlist.append([i,i+1])
                    two_qubit_oper_list.append([i,i+1])
                    two_qubit_qubits.extend([i,i+1])
                    i+=2
                if i == V-1:
                    nlist.append([i])
                    i+=1
        else:
            two_qubit_oper_list = [i for i in nlist if len(i)>1]
            for i,j in two_qubit_oper_list:
                two_qubit_qubits.extend([i,j])

        #Start with a random tree
        G = SimpleGraphLMQC(graphs.RandomTree(V))

        #Randomly choose a number of edges of the graph, we will add those later
        max_number_edges = choice(range(V-1,V*(V-1)))

        #Initialize the complement of G for later use.
        Gc = G.complement()

        #While G has less number of edges than the allowed maximum number of edges,
        #add a random new edge from the edges not yet in G (so they are in Gc)
        while G.size()<max_number_edges and Gc.size()>0:
            i,j,_ = Gc.random_edge()
            G.add_edge(i,j)
            Gc.delete_edge(i,j)

        #Now we have G, we want a local equivalant graph H
        H = copy(G)

        #We make a list of all possible operations, the first part is a list of
        #all possible vertices to do local complementations,
        #the second part is a list of all possible nodes with multiple qubits
        #to do multiqubit operations.

        #We value local complementations a bit more than edge swaps.
        all_operations = 3*[(0,v) for v in range(V)]+[(1,l) for l in two_qubit_oper_list]

        U = []
        #We do random number of random operations on the graph
        for _ in range(randint(V,V^2)):
            #oper is a tuple where the first element is 0 or 1 and the second element
            #is a vertice or a pair, depending on the first.
            oper = choice(all_operations)

            if oper[0]:
                H.flip_edge(oper[1])
                U.append(("CZ",oper[1]))
            else:
                U.append(("LC",oper[1]))
                for i in two_qubit_qubits:
                    if i in H.neighbors(oper[1]):
                        U.append(("LC_N",i))
                H.tau(oper[1],inplace=True)
        return G,H,nlist,U

    @staticmethod
    def find_F_list():
        V = 4
        H = SimpleGraphLMQC({0:[3],1:[2]})
        F_list = []
        for i in powerset(combinations(range(V),2)):
            G = SimpleGraphLMQC(V)
            G.set_partition([[0,1],[2],[3]])
            G.add_edges(i)
            if G.is_LMQC_eq(H)[0]:
                G.relabel({0:"Aa",1:"Ab",2:"A1",3:"A0"})
                F_list.append(G)
        return F_list

    @staticmethod
    def do_GT(G,F,a,b):
        relabel = False
        if not [i for i in G.neighbors(a) if i!=b]:
            if [i for i in G.neighbors(b) if i!=a]:
                a,b=b,a
                relabel = True
            else:
                raise ValueError('a,b are disconnected from other parts of G')
        FG = SimpleGraph(G.union(F))
        FG.add_edges([(a,'A0'),(b,'A1')])
        if G.has_edge(a,b) and F.has_edge('A0','A1'):
            # print "This part is used"
            FG.tau_seq([a,'A0',a],inplace=True)
            #Are we allowed to pick a qubit in another multi-qubit node?
            n_G_of_A1 = [i for i in FG.neighbors('A1') if i in list(set(range(G.order()))-set([a,b]))]
            u = n_G_of_A1[0]
            FG.tau_seq([u,'A1',u],inplace=True)
            n_G_of_b = [i for i in FG.neighbors(b) if i in ['Aa','Ab']]
            v = n_G_of_b[0]
            FG.tau_seq([v,b,v],inplace=True)
            if FG.has_edge('A0','A1'):
                raise ValueError('Two Y measurements will not commute!')
        else:
            FG.tau_seq([a,'A0',a,b,'A1',b],inplace=True)
        if relabel:
            a,b=b,a
        return FG

    def is_LMQC_eq_GT(self,H,F_list = None,debug=False):
        multi_qubit_nodes = [i for i in self.partition if len(i)>1]
        if self.is_LC_eq(H,allow_disc=True):
            return True
        eq = True
        Gp_list = []
        Gp_list.append(self)
        if F_list == None:
            F_list = SimpleGraphLMQC.find_F_list()
        for a,b in multi_qubit_nodes:
            meas = ([a,b,'A0','A1'],['Z','Z','Z','Z'])
            tmp_list = []
            for Gp in Gp_list:
                if debug:
                    print "Set of neighbors connected to a,b",set(Gp.neighbors(a)).union(Gp.neighbors(b)) - set([a,b])
                if not set(Gp.neighbors(a)).union(Gp.neighbors(b)) - set([a,b]):
                    #if a,b are disconnected from the other vertices in Gp
                    Gp.flip_edge((a,b))
                    tmp_list.append(Gp)
                    if debug: print "Gp with flipped edge has been added"
                    continue
                else:
                    if debug: print "Gp with flipped edge has not been added"
                    pass

                if debug:
                    print "(a,b) in G and Gp",G.has_edge(a,b),Gp.has_edge(a,b)
                for i,F in enumerate(F_list):
                    #This comes from the assumption that HSZ=HSZsqrt(-iZ)
                    if i in [0, 1, 2, 4, 5, 6, 8, 10, 12, 16]:
                        pass
                    else:
                        continue
                    G_F = SimpleGraphLMQC.do_GT(Gp,F,a,b)
                    if not Gp.has_edge(a,b) and G_F.has_edge('A0','A1'):
                        raise ValueError('One of the used lemmas is not satisfied')
                    if Gp.has_edge(a,b) and not G_F.has_edge('A0','A1') and not F.has_edge('A0','A1'):
                        raise ValueError('One of the used lemmas is not satisfied')

                    if F.has_edge('A0','A1') and Gp.has_edge(a,b):
                        #Case (I,I)
                        G_oper = copy(G_F)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case II for (0,1) in F and (a,b) in Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case 2 (S,I)
                        G_oper = copy(G_F)
                        G_oper.tau('A0',inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case SI for (0,1) in F and (a,b) in Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case 3 (I,S)
                        G_oper = copy(G_F)
                        G_oper.tau('A1',inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case IS for (0,1) in F and (a,b) in Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case 4 (S,S)
                        G_oper = copy(G_F)
                        G_oper.tau_seq(['A0','A1'],inplace=True)
                        if debug:
                            if G_oper.has_edge('A0','A1'):
                                print "(A0,A1) are connected"
                                print "Edges of F",F.edges()
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case SS for (0,1) in F and (a,b) in Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                    elif F.has_edge('A0','A1') and not Gp.has_edge(a,b):
                        #Case (I,I)
                        G_oper = copy(G_F)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case II for (0,1) in F and (a,b) notin Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case 2 (S,I)
                        G_oper = copy(G_F)
                        G_oper.tau('A0',inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case SI for (0,1) in F and (a,b) notin Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case 3 (I,S)
                        G_oper = copy(G_F)
                        G_oper.tau('A1',inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case IS for (0,1) in F and (a,b) notin Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case 4 (S,S)
                        G_oper = copy(G_F)
                        G_oper.tau_seq(['A0','A1'],inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case SS for (0,1) in F and (a,b) notin Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()
                    elif not F.has_edge('A0','A1') and Gp.has_edge(a,b):
                        #Case (I,I)
                        G_oper = copy(G_F)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case II for (0,1) notin F and (a,b) in Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case 2 (S,I)
                        G_oper = copy(G_F)
                        G_oper.tau('A0',inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case SI for (0,1) notin F and (a,b) in Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case 3 (I,S)
                        G_oper = copy(G_F)
                        G_oper.tau('A1',inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case IS for (0,1) notin F and (a,b) in Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case 4 (S,S)
                        G_oper = copy(G_F)
                        G_oper.tau('A1',inplace=True)
                        u = [i for i in G_oper.neighbors('A0') if i not in [a,b,'A0','A1']][0]
                        G_oper.tau_seq(['A0',u,'A0'],inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case SS for (0,1) notin F and (a,b) in Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()
                    elif not F.has_edge('A0','A1') and not Gp.has_edge(a,b):
                        #Case (I,I)
                        G_oper = copy(G_F)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case II for (0,1) notin F and (a,b) notin Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case (S,I)
                        G_oper = copy(G_F)
                        G_oper.tau('A0',inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case SI for (0,1) notin F and (a,b) notin Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case (I,S)
                        G_oper = copy(G_F)
                        G_oper.tau('A1',inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case IS for (0,1) notin F and (a,b) notin Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()

                        #Case 4 (S,S)
                        G_oper = copy(G_F)
                        G_oper.tau_seq(['A0','A1'],inplace=True)
                        G_oper.meas_seq(meas[0],meas[1],inplace=True)
                        G_oper.relabel({"Aa":a,"Ab":b})
                        tmp_list.append(G_oper)
                        if debug:
                            print "Case SS for (0,1) notin F and (a,b) notin Gp"
                            if G_oper.is_LC_eq(H,allow_disc=True):
                                print F.edges()
            Gp_list = tmp_list
            tmp_list = []

        eq = False
        for G_op in Gp_list:
            if G_op.is_LC_eq(H,allow_disc=True):
                if debug: show(G_op,H)
                eq = True
                break
        # return result.val
        return eq

    def run_tests(self):
        assert type(self) in (SimpleGraphLMQC,SimpleGraph,Graph)

        G = SimpleGraphLMQC(Graph({0:[1,2,3]}),**{'partition':[[0,1],[2],[3]]})
        H = SimpleGraphLMQC(Graph({0:[1,3],2:[3]}))
        dum, _ = G.is_LMQC_eq(H)
        if dum: print "Succesful tested with two LMQC equivalent graphs"
        else: print "Something went wrong"

        G = SimpleGraphLMQC(Graph({0:[1,2,3]}),**{'partition':[[0,1],[2],[3]]})
        H = SimpleGraphLMQC(Graph({0:[3],2:[1]}))
        dum, _ = G.is_LMQC_eq(H)
        if not dum: print "Succesful tested with two not-LMQC equivalent graphs"
        else: print "Something went wrong"
