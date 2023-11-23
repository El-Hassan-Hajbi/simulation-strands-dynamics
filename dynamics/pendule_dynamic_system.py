import numpy as np
from scipy.optimize import fsolve
from .abstract_dynamic_system import AbstractDynamicSystem

g, l = 9.8, np.sqrt(0.5**2 + 1**2)

def F(X, t):
	return np.array([X[1], -(g/l)*np.sin(X[0])],np.float64)

def implicit_step(X, delta, t):
    # Define the function to find the new state implicitly
    def equation(Y):
        return X - Y + delta * F(Y, t)
    
    # Use a numerical solver to find the new state
    X_new = fsolve(equation, X)
    return X_new

## Dummy dynamic system just to test
class PenduleDynamicSystem(AbstractDynamicSystem):

    def __init__(self, mesh, theta, scheme):
        ## Constructor
        # @param self
        # @param mesh  
        # @param theta0  
        super().__init__()
        self.mesh = mesh
        self.scheme = scheme
        # Animations parameters
        self.it = 60.
        self.delta = 0.01
        self.period = 120.
        self.theta0 = np.arctan(self.mesh.positions[1]/self.mesh.positions[2])
        self.X = np.array([self.theta0,0],np.float64) # [theta, theta']
        

    def step(self):
        if self.scheme=="explicit":
            self.X += self.delta * F(self.X, self.it)
        elif self.scheme=="implicit":
            self.X = implicit_step(self.X, self.delta, self.it)
        theta = self.X[0]
        self.mesh.positions[2], self.mesh.positions[3] = l*np.sin(theta), l*(1-np.cos(theta))
		
        self.it += self.delta