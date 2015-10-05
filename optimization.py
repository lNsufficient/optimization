#AUTHORS
#Edvard Johansson
#David Petersson
#Linnéa Stövring-Nielsen
#John Hellborg



#Notes:
#g^k = g(x^k)
#slide 3.14 - there is an alpha missing for calculating x^(k+1)?
#naming convention a^(k+1) is called a_kp1 (as in Plus), and a^(k-1) a_km1 (as in 
#Minus)

import numpy as N

class OptimizationProblem(Object):

    def __init__(self, function, gradient=none): #This should be okay even for
        self.function = function #(...) subclasses, if they don't send any 
        if (gradient == none):
            gradient = computeGradient(function);
        self.gradient = gradient #(...) gradient?

    def computeGradient(function, x):
        h = 1e-7
        d = N.eye(N.size(x))*h
        g = array([lambda x: (f(x+d[:,i])-f(x-d[:,i]))/(2*h) for i in range(N.size(x))])

class Optimization(Object):
    def __init__(self, problemObject):
        self.function = problemObject.function
        self.gradient = problemObject.gradient

    def __call__(self, x0):
        minimize(x0)

    
    def minimize(self, x0):
        H_k = setupHessian(x0) #what will be done the first run?!
                #what is g_km1? x_k = x0, g_k = gradient(x_k)
                #what will H_k be the first run?
        self.g_k = self.gradient(x0)
        x_k = x0

        while (true):
            #Compute s^k = -H^k*g^k
            if (self.quasiNewton == false):
                s_k = self.newtonDirection(H_k, self.g_k)
            else:
                s_k = -H_k*self.g_k
            #Line search for alpha^k
            alpha_k = self.lineSearch(x_k, s_k)
            (self.x_k, self.x_km1)= (x_k+alpha_k*s_k, x_k)
            (self.g_k, self.g_km1) = (self.gradient(x_k), self.g_k) #Was very  convenient to make g_k attribute

            if (self.quasiNewton == false):
                H_k = setupHessian(x_k)
            else:
                H_k = self.hessian(H_k, gamma_k(), delta_k())                 
            if (finished = true):
                break

    def gamma_k(self): #This will only be needed once every run, so it was 
                        #decided not to make this an attribute
        return self.g_k - self.g_km1

    def delta_k(self):
        return self.x_k - self.x_km1

class Newton(Optimization): #This should probably inherit from QuasiNewton instead, even though it feels strange.
    def __init__(self, function, gradient = none, isExact = false):
        super().__init__(function, gradient, isExact)
        self.isExact = isExact
    
    def setupHessian(self, x_k):
        G = numericHessian(x_k) #Visst sa claus att vi skulle hitta fkn?
        self.cholesky = scipy.linalg.cho_factor(G)
        
            #raise ValueError('Hessian is not positive definite')
            #This will never happen - the decomposition will fail instead,
            #then python will raise an error itself.
        return 1/2*(G.conj()+ G.T.conj())

    def numericHessian(self, x)
        h = 1e-7
        d = N.eye(size(x))*h
        H  = [[lambda x: (self.gradient(x+d[:,j])[i]-self.gradient(x-d[:,j])[i])/(2*h) for j in range(size(x))] for i in range(size(x))] 
        return H

    def newtonDirection(self, G):
        #VERY IMPORTANT: Make sure that setup hessian has been run before running this, otherwise this will be super incorrect.        


        #Ginv*g => G*a = g, a = Ginv*g = cho_solve(cho_factor(G), g)
        return scipy.linalg.cho_solve(self.cholesky,g_k) 

    def lineSearch(self): #Since this is just an optimization class
        if (self.isExact == true):
            return exactLineSearch()
        else
            return inexactLineSearch()
        #Here goes the line search algorithm! :) 
        #It is okay for this method to call for hessians and things like that,
        #since the classes that will be created themself contain the hessian
        #method (even though this does not contain it)

    def inexactLineSearch()
    
    def exactLineSearch()

class GoodBroyden(Newton):
    #Rank 1 update -- see wikipedia broyden's method

    def hessian(self, H, gamma, delta):
        #Sherman Morrison formula, see 3.17 - do this by hand
        u = delta - N.dot(H,gamma)
        return H + 1/N.dot(u,gamma)*N.outer(u,u)
    
class BadBroyden(Newton):
    #Every step, caclulate Q(k), then use G^-1 = H, where H = Q(k)^-1 -- see wikipedia broyden's method 

    def hessian(self, H, gamma, delta):
        #Find Q as described in 3.15
        #H = inv(Q)
        return H + N.dot((delta - N.dot(H,gamma))/(N.dot(gamma, gamma)),gamma)

class DFP(Newton) #Uses 3.18 for H

    def hessian(self, H, gamma, delta):
#        H = self.H #This increases readability a lot, and since they are just objects, it will probably not waste a lot of computational power
 #       gamma = self.gamma
  #      delta = self.delta
        return H + N.outer(delta, delta)/N.dot(delta, gamma) - N.dot(N.dot(H,gamma),N.dot(gamma,H))/N.dot(gamma,N.dot(H,gamma))

class BFGS(Newton)

    def hessian(self, H, gamma, delta):
    #    H = self.H
     #   gamma = self.gamma
      #  delta = self.delta
        return H + (1 + N.dot(N.dot(gamma,H),gamma)/N.dot(delta,gamma))*N.outer(delta,delta)/N.dot(delta,gamma)-(N.outer(delta,N.dot(gamma,H)) + N.outer(N.dot(H,gamma),delta))/dot(delta,gamma)
