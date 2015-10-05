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

    def lineSearch(self,x_k,s_k,f_bar): #Since this is just an optimization class
        if (self.isExact == true):
            return exactLineSearch(x_k,s_k,f_bar)
        else
            return inexactLineSearch(x_k,s_k,f_bar)
        #Here goes the line search algorithm! :) 
        #It is okay for this method to call for hessians and things like that,
        #since the classes that will be created themself contain the hessian
        #method (even though this does not contain it)

    def _inexactLineSearch_(x_k,s_k,f_bar,alphai=1):
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
                (alphai,alpha_prev) = (self._choose_(2*alphai-alpha_prev,min(mu,alphai+tau1*(alphai-alpha_prev)),alphai)
            f_prev=f
            
    def _NextIteration_(aj,bj,i,x_k,s_k,f0):   
        somethinglarge=100
        tau2=0.1
        tau3=1/2
        epsilon=10**(-10)
        for j in range(i, somethinglarge):
            alphaj=self._choose_(aj+tau2*(bj-aj),bj-tau3*(bj-aj))
            f=self.function(x_k+alphaj*s_k)
            fderivaj=self.fderive(aj)
            if (aj-alphaj)*fderivaj<=epsilon):
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

    def _choose_(minaj,maxaj):
        return (maxaj+minaj)/2

    def _derive_(x,s): 
        ''' 
        Returns the derivative of (self.)function with respect to alpha
        '''
        h=10**(-10)
        snorm =N.linalg.norm(s)
        return lambda alpha: (self.function(x+(alpha+h)*s)-self.function(x+(alpha-h)*s))/(2*h*snorm)
        

    def _exactLineSearch_(x_k,s_k,f_bar):
        fderive=_derive_(x_k,s_k)
        mu=(f_bar-self.function(x_k)/(self.rho*fderive(x_k))
        alpha = N.linspace(0,mu,1000)
        test = lambda alpha: self.function(x_k+alpha*s_k)
        for i in range(0,1000):
            alpha[i]=test(alpha[i])
        return mu*N.argmin(alpha)/1000

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
