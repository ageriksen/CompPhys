"""
Callable class to compute the unknown function "u" of a scaled and Taylor expanded
differential expansion, given the matrix representation Au = d.
Here, A is taken as a tri-diagonal matrix, with the diagonal and it's upper and 
lower neighbours are the only nonzero elements. these 3 are represented by the vectors
a, b and c. This forms a linear equation set, which we can solve by a forward and backwards
reduction.
"""


class TriSubstitution:
    """Substitution of a Tri-diagonal matrix to solve the
    matrix equation Au = d, with a, b, c and d a known."""

    def __init__(self, u):
        self.u = u
    
    def __call__(self, a, b, c, d):
        import numpy as np
        self.a = a
        self.b = b
        self.c = c
        self.d = d
        
        #creating empty arrays
        #for substitution
#        bT = np.zeros(len(b))
#        dT = np.zeros(len(d)) 
#        bT[0] = b[0]
#        dT[0] = d[0]
#        self.bT = bT
#        self.dT = dT

    def Substitute(self, n):
        """ Forward substitution of the tridiagonal"""
        #forward substitution
        n = int(n)
        #TIME BEGIN
        for i in range(2, n+1):
            #precalculate repeated operations
#            print(
#                'i = ', i, '\n',
#                'a[',i-1,'] = ', self.a[i-1], '\n',
#                'b[',i-1,'] = ', self.b[i-1], '\n',
#                    )
            e = self.a[i-1]/float(self.b[i-1])
            self.b[i] = self.b[i] - e*self.c[i-1]
            self.d[i] = self.d[i] - e*self.d[i-1]
        #backward substitution
        self.u[n] = self.d[n]/self.b[n]
        for i in range(n-1, 0, -1):
            self.u[i] = (self.d[i] - self.c[i]*self.u[i+1])/self.b[i]
        #TIME END
        # the 1st loop has 5(n-2) FLOPS and the 2nd has 3(n-1) FLOPS.
        # this makes a total 8n-13 FLOPS, which for large n approximates to 
        # FLOPS: 8n
        FLOPS = 8*n
        return self.u, FLOPS
