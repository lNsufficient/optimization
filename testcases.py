#AUTHORS
#Edvard Johansson
#David Petersson
#Linnéa Stövring-Nielsen
#John Hellborg



import unittest
import numpy as N
import optimization as O

class testOptimization(unittest.TestCase):

    def testSecondDegreePolynomialNewtonExact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.Newton(op,True)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        print(xmin)
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.001)

    def testSecondDegreePolynomialNewtonInexact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.Newton(op,False)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        print(xmin)
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.001)


    def testSecondDegreePolynomialGoodBroydenExact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.GoodBroyden(op,True)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        print(xmin)
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.001)

    def testSecondDegreePolynomialGoodBroydenInexact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.Newton(op,False)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        print(xmin)
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.001)



if __name__== '__main__':
    unittest.main()
