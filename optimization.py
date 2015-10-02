#Notes:
#g^k = g(x^k)
#slide 3.14 - there is an alpha missing for calculating x^(k+1)?


class Optimization(Object):
    def __init__(self, function, gradient=none): #This should be okay even for
        self.function = function #(...) subclasses, if they don't send any 
        self.gradient = gradient #(...) gradient?

    def __call__(self, x0):
        minimize(x0)
    
    def minimize(self, x0):
        #Compute s^k = -H^k*g^k
        #Line search for alpha^k
        alpha = lineSearch()
        #Compute x^(k+1) = x^k+alpha^k*s^k
        #Find next H^k

    def lineSearch(self): #Since this is just an optimization class
        if (self.isExact == true):
            return exactLineSearch()
        else
            return inexactLineSearch()
        #Here goes the line search algorithm! :) 
        #It is okay for this method to call for hessians and things like that,
        #since the classes that will be created themself contain the hessian
        #method (even though this does not contain it)



class QuasiNewton(Optimization):
    def __init__(self, function, gradient=none, isExact=false):
        super().__init__(function, gradient, isExact)
        #Sedan behöver inte klasserna som ärver av denna ha en egen init!

    def __init__(self, function, gradient = none, isExact = false):
        super().__init__(function, gradient, isExact)
        self.isExact = isExact

class Newton(QuasiNewton): #This should probably inherit from QuasiNewton instead, even though it feels strange.
    
    def hessian(self, ):
        #This is done in different ways depending on the class...inheritance
        #decided to call it hessian even here because after all G^-1 will be calculated... 
    

class GoodBroyden(QuasiNewton):
    #Rank 1 update -- see wikipedia broyden's method

    def hessian(self):
        #Sherman Morrison formula, see 3.17 - do this by hand
        self.H = H + /(N.dot(delta,gamma)-N.dot(N.dot(H,gamma),gamma))
    
class BadBroyden(QuasiNewton):
    #Every step, caclulate Q(k), then use G^-1 = H, where H = Q(k)^-1 -- see wikipedia broyden's method 

    def hessian(self):
        #Find Q as described in 3.15
        #H = inv(Q)

class DFP(QuasiNewton) #Uses 3.18 for H

    def hessian(self):
        H = self.H #This increases readability a lot, and since they are just objects, it will probably not waste a lot of computational power
        gamma = self.gamma
        delta = self.delta
        self.H = H + N.outer(delta, delta)/N.dot(delta, gamma) - N.dot(N.dot(H,gamma),N.dot(gamma,H))/N.dot(gamma,N.dot(H,gamma))
        return self.H

class BFGS(QuasiNewton)

    def hessian(self):
        H = self.H
        gamma = self.gamma
        delta = self.delta
        self.H = H + (1 + N.dot(N.dot(gamma,H),gamma)/N.dot(delta,gamma))*N.outer(delta,delta)/N.dot(delta,gamma)-(N.outer(delta,N.dot(gamma,H)) + N.outer(N.dot(H,gamma),delta))/dot(delta,gamma)
