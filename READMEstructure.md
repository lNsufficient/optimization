class optimization
    instructions:
    tell the user the attributes available, so that the user understands how
     they can be changed and what difference that change will make


    ATTRIBUTES:
    function
    interval
    starting values alpha1
    xi 
    thau
    rho
    isExact=false (by default)
    accuracy=something small

    FUNCTIONS:
    init(self, function, isExact=false, xi=default, thau=default, rho=default, accuracy=default) :
        set attributes
        
    call: 
        inparameter:
        accuracy
        optimization(alpha-guess)
        
        evaluate function depending on the value of isInexact
            if isExact
                exactLineSearch
            else
                inexactLineSearch 
    
    setMode:
        isInexact = input
   
    _exactLineSearch_:

    _inExactLineSearch_:

    gradient:

    appriximatedHessian:

    setParameters:
        in case we would realize that parameters somehow depends on eachother, this is how we would recommend that the user changes them

    
