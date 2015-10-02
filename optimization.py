#Notes:
#g^k = g(x^k)
#slide 3.14 - there is an alpha missing for calculating x^(k+1)?
#naming convention a^(k+1) is called a_kp1 (as in Plus), and a^(k-1) a_km1 (as in 
#Minus)

class Optimization(Object):
    def __init__(self, function, gradient=none): #This should be okay even for
        self.function = function #(...) subclasses, if they don't send any 
        self.gradient = gradient #(...) gradient?

    def __call__(self, x0):
        minimize(x0)
    
    def minimize(self, x0):
        setup() #what will be done the first run?!
                #what is g_km1 and x_km1? x_k = x0, g_k = gradient(x_k)
                #what will H_k be the first run?

        #Compute s^k = -H^k*g^k
        (self.g_k, self.g_km1) = (self.gradient(x_k), self.g_k) #Was very 
        H_k = self.hessian(H_k, gamma_k(), delta_k())             #convenient to make g_k attribute
        s_k = -H_k*self.g_k
        #Line search for alpha^k
        alpha_k = self.lineSearch()
        (self.x_k, self.x_km1)= (x_k+alpha_k*s_k, x_k)
        #Find next H^k
        self.H = self.hessian()

    def gamma_k(self): #This will only be needed once every run, so it was 
                        #decided not to make this an attribute
        return self.g_k - self.g_km1

    def delta_k(self):
        return self.x_k - self.x_km1

class Newton(Optimization): #This should probably inherit from QuasiNewton instead, even though it feels strange.
    def __init__(self, function, gradient = none, isExact = false):
        super().__init__(function, gradient, isExact)
        self.isExact = isExact
    
    def hessian(self, H_k, gamma, delta):
        H = some.function.hessian(H_k) #Visst sa claus att vi skulle hitta fkn?
        G = scipy.linalg.cho_factor(H)
        
            #raise ValueError('Hessian is not positive definite')
            #This will never happen - the decomposition will fail instead,
            #then python will raise an error itself.
        return 1/2*(G.conj()+ G.T.conj())

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

    def hessian(self):
        #Sherman Morrison formula, see 3.17 - do this by hand
        return H + /(N.dot(delta,gamma)-N.dot(N.dot(H,gamma),gamma))
    
class BadBroyden(Newton):
    #Every step, caclulate Q(k), then use G^-1 = H, where H = Q(k)^-1 -- see wikipedia broyden's method 

    def hessian(self):
        #Find Q as described in 3.15
        #H = inv(Q)

class DFP(Newton) #Uses 3.18 for H

    def hessian(self):
        H = self.H #This increases readability a lot, and since they are just objects, it will probably not waste a lot of computational power
        gamma = self.gamma
        delta = self.delta
        return H + N.outer(delta, delta)/N.dot(delta, gamma) - N.dot(N.dot(H,gamma),N.dot(gamma,H))/N.dot(gamma,N.dot(H,gamma))

class BFGS(Newton)

    def hessian(self):
        H = self.H
        gamma = self.gamma
        delta = self.delta
        return H + (1 + N.dot(N.dot(gamma,H),gamma)/N.dot(delta,gamma))*N.outer(delta,delta)/N.dot(delta,gamma)-(N.outer(delta,N.dot(gamma,H)) + N.outer(N.dot(H,gamma),delta))/dot(delta,gamma)
