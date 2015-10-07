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
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.00001)

    def testSecondDegreePolynomialNewtonInexact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.Newton(op,False)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.00001)


    def testSecondDegreePolynomialGoodBroydenExact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.GoodBroyden(op,True)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.00001)

    def testSecondDegreePolynomialGoodBroydenInexact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.GoodBroyden(op,False)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.00001)

    def testSecondDegreePolynomialBadBroydenExact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.BadBroyden(op,True)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.00001)

    def testSecondDegreePolynomialBadBroydenInexact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.BadBroyden(op,False)
        xmin = minimize(N.array([1.,2.5,0]))
        expected = N.array([1,3,-1])
        print('FAIL'+str(xmin))
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.1)

    def testSecondDegreePolynomialDFPExact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.DFP(op,True)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.00001)

    def testSecondDegreePolynomialDFPInexact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.DFP(op,False)
        xmin = minimize(N.array([1.,2.5,0]))
        expected = N.array([1,3,-1])
        print('FAIL'+str(xmin))
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.00001)

    def testSecondDegreePolynomialBFGSExact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.BFGS(op,True)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.00001)

    def testSecondDegreePolynomialBFGSInexact(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.BFGS(op,False)
        xmin = minimize(N.array([1.,2.5,0]))
        expected = N.array([1,3,-1])
        print('FAIL'+str(xmin))
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.00001)




    def testRosenbrockNewtonExact(self):
        f = lambda x: 100*(x[1]-x[0]**2)**2+(1-x[0])**2
        op = O.OptimizationProblem(f)
        minimize = O.Newton(op,True)
        xmin = minimize(N.array([0.,0]))
        expected = N.array([1.,1])
        for i in range(2):
            self.assertTrue(abs(xmin[i]-expected[i])<0.001)

    def testRosenbrockNewtonInexact(self):
        f = lambda x: 100*(x[1]-x[0]**2)**2+(1-x[0])**2
        op = O.OptimizationProblem(f)
        minimize = O.Newton(op,False)
        xmin = minimize(N.array([0.,0]))
        expected = N.array([1.,1])
        for i in range(2):
            self.assertTrue(abs(xmin[i]-expected[i])<0.0000001)





if __name__== '__main__':
    unittest.main()
