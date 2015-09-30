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

    _exactLineSearch_():
        pass

    _inexactLineSearch_():    
        pass

    approximatedHessian():
        pass

    gradient():

    #setParameters():
