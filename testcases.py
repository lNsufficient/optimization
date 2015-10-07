#AUTHORS
#Edvard Johansson
#David Petersson
#Linnéa Stövring-Nielsen
#John Hellborg



import unittest
import numpy as N
import optimization as O

class testOptimization(unittest.TestCase):
    def testSecondDegreePolynomialNewton(self):
        f = lambda x: (1-x[0])**2+(x[1]-3)**2+(x[2]+1)**2
        op = O.OptimizationProblem(f)
        minimize = O.Newton(op,True)
        xmin = minimize(N.array([0.,0,0]))
        expected = N.array([1,3,-1])
        print(xmin)
        for i in range(3):
            self.assertTrue(abs(xmin[i]-expected[i])<0.001)

#,msg="Newton can not minimize second degree polynomial")


if __name__== '__main__':
    unittest.main()
