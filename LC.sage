def M2vec(M):
    """Returns a vectorization of M, column-stacking."""
    return vector(GF(2),sum([list(row) for row in M.transpose()],[]))

def E(i,n):
    """Returns the matrix E_i which has a single nonzero value 1 at pos (i,i)"""
    tmp=zero_matrix(GF(2),n)
    tmp[i,i]=1
    return tmp

def E_tim_M(M,i):
    """Returns E_ii*M where E_i has a single nonzero value 1 at pos (i,i)"""
    n=M.nrows()
    tmp=zero_matrix(GF(2),n)
    tmp.set_row(i,M.row(i))
    return tmp

def M_tim_E(M,i):
    """Returns M*E_ii where E_i has a single nonzero value 1 at pos (i,i)"""
    n=M.nrows()
    tmp=zero_matrix(GF(2),n)
    tmp.set_column(i,M.column(i))
    return tmp

def constr_mat(A1,A2):
    """Returns the matrix corresponding to the equation A1*B*A2+D*A2+A1*A+C=0, but in terms of the vector with the diagonals."""
    if not (A1.is_square() and A2.is_square()):
        raise TypeError("Both matrices not square.")
    n=A1.nrows()
    if not n==A2.nrows():
        raise TypeError("Matrices not of the same size.")
    mat1=matrix(GF(2),[M2vec(M_tim_E(A1,i)) for i in range(n)]).transpose()
    mat2=matrix(GF(2),[M2vec(M_tim_E(A1,i)*A2) for i in range(n)]).transpose()
    mat3=matrix(GF(2),[M2vec(E(i,n)) for i in range(n)]).transpose()
    mat4=matrix(GF(2),[M2vec(E_tim_M(A2,i)) for i in range(n)]).transpose()
    return block_matrix(GF(2),[[mat1,mat2,mat3,mat4]],subdivide=False)

def cond(vec):
    """Tests the cond for ABCD to be symplectic transformation AD+BC=I, but in terms of the vector consisting of the matrices diagonals and for the equivalent statement that 2x2 matrices with A_aa and so on if nonsingular"""
    n=vec.length()/4
    return all([not matrix(GF(2),[[vec[a],vec[n+a]],[vec[2*n+a],vec[3*n+a]]]).is_singular() for a in range(n)])

def test_cond(A1,A2,ret_sol=False):
    """Tests if the graphs with adj.matrices A1 and A2 are LC-equivalent, efficiently O(n^4)"""
    V=constr_mat(A1,A2).right_kernel()
    if V.dimension()<=4:
        for vec in V:
            if cond(vec):
                return True
        return False
    else:
        B=V.basis()
        for b1 in B:
            if cond(b1):
                return True
            for b2 in B:
                if cond(b1+b2):
                    return True
        return False
