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
        
    def general(self, n):
        """ forward & backward subst. general tri diagonal matrix. 
        The 1st loop contains 5(n-2) FLOPS and the 2nd 3(n-1). 
        All told, they make 8n-13, or roughly 8n FLOPS"""
        #forward substitution
        n = int(n)
        #TIME BEGIN
        for i in range(2, n+1):
            #precalculate repeated operations
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
        return self.u, FLOPS #, time

    def special(self, n):
        """ forward & backward substitution of special system, where a[i]=c[i]=-1 & b[i]=2.
        By precalculating vector b, we free up to roughly 4n FLOPS total."""
        n = int(n)
        self.b -= 0.5
        #TIME BEGIN
        #forward
        for i in range(2, n+1):
            self.d[i] = self.d[i] + 0.5*self.d[i-1]
        #backwards
        self.u[n] = self.d[n]/self.b[n]
        for i in range(n-1, 0, -1):
            self.u[i] = (self.d[i] - u[i+1])/b[i]
        FLOPS = 4*n
        return self.u, FLOPS #, time
