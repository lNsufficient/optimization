#AUTHORS
#Edvard Johansson
#David Petersson
#Linnéa Stövring-Nielsen
#John Hellborg

#naming convention a^(k+1) is called a_kp1 ("kPlus1"), and a^(k-1) a_km1(Minus)

import numpy as N
import scipy
import scipy.linalg





class OptimizationProblem(object):

    def __init__(self, function, gradient=None): #This should be okay even for
        self.function = function #(...) subclasses, if they don't send any 
        self.h = 1e-9
        self.gradient = gradient

class Optimization(object):
    def __init__(self, problemObject):
        self.function = problemObject.function
        self.gradient = problemObject.gradient
        self.h = problemObject.h
        self.defineConstants()

    def defineConstants(self):
        self.rho = 0.1
        self.xi = 9
        self.sigma = 0.7
        self.tau = 0.1

    def __call__(self, x0):
        x_min = self.minimize(x0)
        print("x_min: " + str(self.x_k))
        return x_min 

    
    def computeGradient(self, x_k):
        h = self.h
        xsize = N.size(x_k)
        d = N.eye(xsize)*h
        g = N.array(list([lambda x, i=i: (self.function(x+d[i,:])-self.function(x-d[i,:]))/(2*h) for i in range(xsize)]))
        return N.array(g)

    def minimize(self, x0):
        #This is needed to get the while loop running
        if (self.gradient == None): #Could not find any other way to know its 
            self.gradient = self.computeGradient(x0) #dim, else than this.
        if (self.quasiNewton == False):
            H_k = self.setupHessian(x0) 
        else: 
            H_k = N.eye(N.size(x0))
        self.g_k = [self.gradient[i](x0) for i in range(N.size(x0))]
        self.x_k = x0
        tol = self.h

        #This "The loop"
        while (True):
            #print("g_k" + str(self.g_k))
            print(self.x_k)
            #Compute s^k = -H^k*g^k
            s_k = self.newtonDirection(H_k, self.g_k)
            if s_k is None:
                print("s_k is nonetype")

#This was commented 13/10 14:38
#            if (self.quasiNewton == False):
 #               s_k = self.newtonDirection(H_k, self.g_k)
  #          else:
   #             s_k = -N.dot(H_k,self.g_k)
            #Line search for alpha^k

            alpha_k = self.lineSearch(self.x_k, s_k, tol)
            if (alpha_k == 0):
                print("Alpha == 0!")    
                print(x_k)
                alpha = 0.1
                
            #if abs(alpha_k)<=0.0000001:
             #   print(N.linalg.norm(alpha_k*s_k))
                #alpha_k=0.1
            if alpha_k is None:
                print("alpha_k is None")

                print("H_k" + str(H_k))
                print("s_k" + str(s_k))
                print("x_k" + str(self.x_k))
                print("g_k" + str(self.g_k))

            (self.x_k, self.x_km1)= (self.x_k+alpha_k*s_k, self.x_k)
            (self.g_k, self.g_km1) = (N.array([self.gradient[i](self.x_k) for i in range(N.size(self.x_k))]), self.g_k)
            if (self.quasiNewton == False):
                H_k = self.setupHessian(self.x_k)
            else:
                H_k = self.hessian(H_k, self.gamma_k(), self.delta_k()) 
            if (N.linalg.norm(alpha_k*s_k) < tol and ((self.function(self.x_k)-self.function(self.x_k+alpha_k*s_k))<tol)):
                break
            print(self.x_k)
            print(self.function(self.x_k))
        return self.x_k

    def gamma_k(self): 
        return self.g_k - self.g_km1

    def delta_k(self):
        return self.x_k - self.x_km1

class Newton(Optimization): 
    def __init__(self, oP, isExact=False):
        super().__init__(oP)
        self.isExact = isExact
        self.quasiNewton = False
    
    def setupHessian(self, x_k):
        G = self.numericHessian(x_k)
        G = N.array([[self.numericHessian(x_k)[i,j](x_k) for i in range(N.size(x_k))] for j in range(N.size(x_k))]) 
        G = 1/2*(G.conj()+ G.T.conj())
        return G

    def numericHessian(self, x):
        h = self.h
        d = N.eye(N.size(x))*h
        H  = N.array([[lambda x, i=i, j=j: (self.gradient[i](x+d[:,j])-self.gradient[i](x-d[:,j]))/(2*h) for j in range(N.size(x))] for i in range(N.size(x))])
        return H

    def newtonDirection(self, G, g_k):
        #VERY IMPORTANT: Make sure that setup hessian has been run before running this, otherwise this will be super incorrect.        
         
        try:
            s_k = -scipy.linalg.cho_solve(scipy.linalg.cho_factor(G),g_k)
        except N.linalg.linalg.LinAlgError:
            LU = scipy.linalg.lu_factor(G)
            s_k = scipy.linalg.lu_solve(LU, g_k)
        return s_k 

    def lineSearch(self,x_k,s_k,f_bar): 
        if (self.isExact == True):
            return self._exactLineSearch_(x_k,s_k,f_bar)
        else:
            return self._inexactLineSearch_(x_k,s_k,f_bar)

    def _inexactLineSearch_(self, x_k,s_k,f_bar,alphai=1):
        '''
        Calculate the alpha which minimizes the function supplied to the class in init

        Args:
            x_k: A vector containing the x values for which alpha is to be calculated.
            s_k: A vector containing the direction to look in.
            f_bar: The tolerance level, meaning all alpha that return a function value lower or equal to f_bar are ok.
            alphai: The starting value of alpha, with default value 1, used to calculate the correct alpha.
        Returns:
            The alpha which minimizes the function supplied to the class
        '''
        self.fderive = self._derive_(x_k,s_k) #Self????
        SomethingLarge=100
        alpha_prev = 0
        f_prev = self.function(x_k)
        f0=self.function(x_k)
        self.fderive0 = self.fderive(0) #Self???
        if self.fderive0==0:
            mu = 10*(f_bar-f0)/self.rho
        else:
            mu=(f_bar-f0)/(self.rho*self.fderive0) 
        for i in range(0,SomethingLarge):
            f=self.function(x_k+alphai*s_k)
            if f<=f_bar:
                return alphai
            if (f>f0+alphai*self.rho*self.fderive0) or f>=f_prev: #could be wrong, fletcher...f'(0)
                ai=alpha_prev
                bi=alphai
                #terminate block end
                return self._NextIteration_(ai,bi,i,x_k,s_k,f0)
            f_deriv = self.fderive(alphai)
            if abs(f_deriv)<=-self.sigma*self.fderive0:
                return alphai
            if f_deriv>=0:
                ai = alphai
                bi = alpha_prev
                #terminate block end
                return self._NextIteration_(ai,bi,i,x_k,s_k,f0)
            if mu<=2*alphai-alpha_prev:
                (alphai,alpha_prev) = (mu,alphai)
            else:
                tau1=9
                (alphai,alpha_prev) = (self._choose_(2*alphai-alpha_prev,min(mu,alphai+tau1*(alphai-alpha_prev))),alphai)
            f_prev=f
            
    def _NextIteration_(self, aj,bj,i,x_k,s_k,f0):   
        somethinglarge=100
        tau2=0.1
        tau3=1/2
        epsilon=10**(-10)
        for j in range(i, somethinglarge):
            alphaj=self._choose_(aj+tau2*(bj-aj),bj-tau3*(bj-aj))
            f=self.function(x_k+alphaj*s_k)
            fderivaj=self.fderive(aj)
            if ((aj-alphaj)*fderivaj<=epsilon):
                return alphaj
            if f>f0+self.rho*alphaj*self.fderive0 or f>=self.function(x_k+aj*s_k):
                bj=alphaj
            else:
                f_deriv = self.fderive(alphaj) 
                if abs(f_deriv)<=-self.sigma*self.fderive0:
                    return alphaj
                if (bj-aj)*f_deriv>=0:
                    bj=aj
                aj=alphaj
        return alphaj

    def _choose_(self, minaj,maxaj):
        return (maxaj+minaj)/2

    def _derive_(self,x,s): 
        ''' 
        Returns the derivative of (self.)function with respect to alpha
        '''
        h=self.h
        return lambda alpha: (self.function(x+(alpha+h)*s)-self.function(x+(alpha-h)*s))/(2*h)
        

    def _exactLineSearch_(self, x_k,s_k,f_bar):

        mu=100 #Maybe, mu should be chosen in a smarter way than this, but for now, we let it stay like this. 
        test = lambda alpha: self.function(x_k+alpha*s_k) 
        tmin = 0
        tmax = mu
        points = 100
        for j in range(9):
            fpoints = [test(i) for i in N.linspace(tmin, tmax, points)]
            deltat = (tmax-tmin)/points
            (tmin, tmax) = (tmin + deltat*(N.argmin(fpoints)-1), tmin + deltat*(N.argmin(fpoints)+1))
        return  tmin + deltat*N.argmin(fpoints)

class QuasiNewton(Newton):
    def __init__(self, oP, isExact=True):
        super().__init__(oP, isExact)
        self.quasiNewton = True

    def newtonDirection(self, H_k, g_k):
        return -N.dot(H_k, g_k)        

class GoodBroyden(QuasiNewton):
    #Rank 1 update -- see wikipedia broyden's method

    def hessian(self, H, gamma, delta):
        #Sherman Morrison formula, see 3.17 - do this by hand
        u = delta - N.dot(H,gamma)
        return H + 1/N.dot(u,gamma)*N.outer(u,u)
    
class BadBroyden(QuasiNewton):
    #Every step, caclulate Q(k), then use G^-1 = H, where H = Q(k)^-1 -- see wikipedia broyden's method 

    def hessian(self, H, gamma, delta):
        return H + N.dot((delta - N.dot(H,gamma))/(N.dot(gamma, gamma)),gamma)

class DFP(QuasiNewton): 

    def hessian(self, H, gamma, delta):
        return H + N.outer(delta, delta)/N.dot(delta, gamma) - N.outer(N.dot(H,gamma),N.dot(gamma,H))/N.dot(gamma,N.dot(H,gamma))

class BFGS(QuasiNewton):

    def hessian(self, H, gamma, delta):
        return H + (1 + N.dot(N.dot(gamma,H),gamma)/N.dot(delta,gamma))*N.outer(delta,delta)/N.dot(delta,gamma)-(N.outer(delta,N.dot(gamma,H)) + N.outer(N.dot(H,gamma),delta))/N.dot(delta,gamma)

"""
f = lambda x: (1-x[0])**4+x[1]**2-x[2]**2
f = lambda x: 100*(x[1] - x[0]**2)**2+(1-x[0])**2
f = lambda x: 100*(x[0] - x[1]**2)**2+(1-x[1])**2
g = lambda x: N.array([2*x[0], 2*x[1], 2*x[2]])
op = OptimizationProblem(f)
minimize = Newton(op, True)
minimize = GoodBroyden(op, True)
x_min = minimize(N.array([53,-260]))
print(x_min)
"""

