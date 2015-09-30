class QuasiNewton:
    #Description:
    
    __init__(self, function, isExact=false, sigma=0.7, xi=9, thau=0.1, rho=0.1):
    #default values from page 113-114 in Practical Optimization
        self.function = function
        self.isExact = isExact
        self.sigma = sigma
        self.xi = xi
        self.thau = thau
        self.rho = rho

    __call__(self, alpha):

        alpha_k = alpha #?
        if isExact:
            self._exactLineSearch_(alpha_k)
        else:
            self._inexactLineSearch_(alpha_k)

    
    setMode(self, isExact):
        self.isExact = isExact

    _exactLineSearch_(x_k,s_k,f_bar,alphai=1):
        SomethingLarge=100
        alpha_prev = 0
        f_prev = self.function(x_k)
        f0=self.function(x_k)
        mu=(f_bar-f0)/(self.rho**2) #could be wrong, fletcher...f'(0)
        for i in range(0,SomethingLarge):
            f=self.function(x_k+alphai*s_k)
            if f<=f_bar:
                return (x_k,alphai)
            if (f>f0+alphai*self.rho) or f>=f_prev:
                ai=alpha_prev
                bi=alphai
                #terminate block end
                return self._NextIteration_(ai,bi,i,x_k,s_k,f0)
            f_deriv = 5 #This must be changed to f'(alphai)
            if abs(f_deriv)<=-self.sigma*self.rho:
                return (x_k,alphai)
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
            
    _NextIteration_(aj,bj,i,x_k,s_k,f0):   
        somethinglarge=100
        tau2=0.1
        tau3=1/2
        epsilon=10**(-10)
        for j in range(i, somethinglarge):
            alphaj=self._choose_(aj+tau2*(bj-aj),bj-tau3*(bj-aj))
            f=self.function(x_k+alphaj*s_k)
            fderivaj=1 #this must be changed to f'(aj)
            if (aj-alphaj)*fderivaj<=epsilon):
                return (x,aj)
            if f>f0+self.rho**2*alphaj or f>=self.function(x_k+aj*s_k):
                bj=alphaj
            else:
                f_deriv = 5 #this must be changed to f'(alphaj)
                if abs(f_deriv)<=-self.sigma*self.rho:
                    return (x_k,alphaj)
                if (bj-aj)*f_deriv>=0:
                    bj=aj
                aj=alphaj
        return (x,aj)

    _choose_(minaj,maxaj):
        return (maxaj+minaj)/2

    _inexactLineSearch_():    
        pass

    approximatedHessian():
        pass

    gradient():

    #setParameters():
