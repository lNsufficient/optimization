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

    _exactLineSearch_(alphai,x_k,s_k,f_bar):
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
                return self._NewtonIteration_(ai,bi)
            f_deriv = 5 #This must be changed to f'(alphai)
            if abs(f_deriv)<=-self.sigma*self.rho:
                return (x_k,alphai)
            if f_deriv>=0:
                ai = alphai
                bi = alpha_prev
                #terminate block end
                return self._NewtonIteration_(ai,bi)
            if mu<=2*alphai-alpha_prev:
                (alphai,alpha_prev) = (mu,alphai)
            else:
                (alphai,alpha_prev) = (2*alphai-alpha_prev,alphai)
                #alphai belongs to [2*alphai-alpha_prev, min(mu,alphai+tau1*(alphai-alpha_prev))], tau1>1, tau1=9
            f_prev=f
            
    _NewtonIteration_(ai,bi,x_prev,s_k):   
        x=x_prev-s_k 
        return (x,ai)

    _inexactLineSearch_():    
        pass

    approximatedHessian():
        pass

    gradient():

    #setParameters():
