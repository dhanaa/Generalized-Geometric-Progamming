# Generalized-Geometric-Progamming
Generalized Geometric Programming (GGP) Solver for Matlab users

For original paper, please refer to 
Global Optimization in Generalized Geometric Programming, Maranas and Floudas, 1997
http://maranas.che.psu.edu/pub/1997/Maranas_and_Floudas,_Computers_and_Chem._Eng.,_1997.pdf


GGP here means the coefficient before monomials can be negative, which makes the GP non-convex. The definition of GGP here is different
from that in other sources like Boyd's toolbox.

Take a simple example. The following problem is non-convex in terms of x,y:

    min   x

    subject to
  
    0.25 * x + 0.5 * y - (1/16) * x^2 - (1/16)*y^2 -1 <= 0

    (1/14) * x^2 + (1/14) * y^2 + 1 -(3/7) * x - (3/7) * y <=0
    
    1 <= x <= 5.5
    
    1 <= y <= 5.5
    
To use the solver, we have to transform the problem:

    % use vector to represent monomial
    % use struct to represent each line (including objective function)

We express objective function as: 

    obj.positive = [1 1 0];
    obj.negative = 0;

We express constraint 1 as:

    con1.positive = [0.25 1 0; 0.5 0 1];
    con1.negative = [1/16 2 0; 1/16 0 2; 1 0 0];

We express constraint 2 as:

    con2.positive = [1/14 2 0; 1/14 0 2; 1 0 0];
    con2.negative = [3/7 1 0; 3/7 0 1];

We stack all constraint inequalities: 

    constraint = [con1, con2]; % stack all constraints

We need some global upper and lower bounds for variables:

    lowerbound = [1, 1]; % global bounds for variables
    upperbound = [5.5, 5.5];
    
We also tell the Solver the number of unknowns

    N = 2; % number of variables

Lastly, we can call the solver:

    [obj, var] = GGPSolver(N, obj,constraint,lowerbound,upperbound);
