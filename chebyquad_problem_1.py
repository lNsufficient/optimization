# -*- coding: utf-8 -*-
"""
Chebyquad Testproblem

Course Material for the course FMNN25

Created on Wed Nov 23 22:52:35 2011

@author: Claus Führer
"""
from  __future__  import division
from  scipy       import dot,linspace
import scipy.optimize as so
from numpy import array
import optimization as op
import numpy 


def T(x, n):
    """
    Recursive evaluation of the Chebychev Polynomials of the first kind
    x evaluation point (scalar)
    n degree 
    """
    if n == 0:
        return 1.0
    if n == 1:
        return x
    return 2. * x * T(x, n - 1) - T(x, n - 2)

def U(x, n):
    """
    Recursive evaluation of the Chebychev Polynomials of the second kind
    x evaluation point (scalar)
    n degree 
    Note d/dx T(x,n)= n*U(x,n-1)  
    """
    if n == 0:
        return 1.0
    if n == 1:
        return 2. * x
    return 2. * x * U(x, n - 1) - U(x, n - 2) 
    
def chebyquad_fcn(x):
    """
    Nonlinear function: R^n -> R^n
    """    
    n = len(x)
    def exact_integral(n):
        """
        Generator object to compute the exact integral of
        the transformed Chebychev function T(2x-1,i), i=0...n
        """
        for i in range(n):
            if i % 2 == 0: 
                yield -1./(i**2 - 1.)
            else:
                yield 0.
    exint = exact_integral(n)
    
    def approx_integral(i):
        """
        Approximates the integral by taking the mean value
        of n sample points
        """
        return sum(T(2. * xj - 1., i) for xj in x) / n

    exintlist = [i for i in exint]
        
#    return array([approx_integral(i) - exint.next() for i in range(n)]) 
    return array([approx_integral(i) - exintlist[i] for i in range(n)]) 
def chebyquad(x):
    """            
    norm(chebyquad_fcn)**2                
    """
    chq = chebyquad_fcn(x)
    return dot(chq, chq)

def gradchebyquad(x):
    """
    Evaluation of the gradient function of chebyquad
    """
    chq = chebyquad_fcn(x)
    UM = 4. / len(x) * array([[(i+1) * U(2. * xj - 1., i) 
                             for xj in x] for i in range(len(x) - 1)])
    return dot(chq[1:].reshape((1, -1)), UM).reshape((-1, ))


    
if __name__ == '__main__':
    x=linspace(0,1,4)
    xmin= so.fmin_bfgs(chebyquad,x,gradchebyquad)  # should converge after 18 iterations  
    columns = numpy.eye(numpy.size(x,0))
    g = [lambda x, i=i: gradchebyquad(columns[:,i]*x)[i] for i in range(numpy.size(x))]
    print(xmin)
    optProb = op.OptimizationProblem(chebyquad, g) 
    optim = op.BFGS(optProb, True)
    xmin2 = optim(x)
    print("Function value from our method", chebyquad(xmin2))
        
